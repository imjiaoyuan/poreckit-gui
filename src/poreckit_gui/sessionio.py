"""Shared session persistence primitives.

The four views that save session JSON (main browser, sequence browser,
assembly editor and plot workbench) used to hand-roll file writing,
JSON (de)serialization and relative-path handling. This module centralizes
those primitives so every session file is written atomically and paths are
(de)relativized with identical semantics.

Schema, file naming and dialogs stay the responsibility of each view.
"""
import json
import os
import tempfile
from pathlib import Path

from .remote import is_remote


def save_json(path, data):
    """Serialize ``data`` to ``path`` atomically.

    The JSON text is written to a temporary file in the same directory,
    flushed to disk, then moved onto the target with :func:`os.replace`.
    An interrupted save (crash, full disk, exception mid-write) therefore
    leaves the previous file at ``path`` untouched and at worst leaves a
    stray temporary file behind. Any ``TypeError`` from serialization or
    ``OSError`` from the filesystem propagates to the caller.
    """
    path = Path(path)
    text = json.dumps(data, ensure_ascii=False, indent=2)
    directory = path.parent if str(path.parent) else Path(".")
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    try:  # Make the directory entry durable too; best effort on exotic filesystems.
        dirfd = os.open(directory, os.O_RDONLY)
        try:
            os.fsync(dirfd)
        finally:
            os.close(dirfd)
    except OSError:
        pass


def load_json(path):
    """Read and parse the JSON file at ``path``.

    Raises ``OSError`` if the file cannot be read and ``ValueError`` with a
    clear, path-bearing message if the content is not valid JSON. JSON
    decoding errors never surface as opaque ``json.JSONDecodeError`` text
    without file context.
    """
    text = Path(path).read_text(encoding="utf-8")
    try:
        return json.loads(text)
    except ValueError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc


def to_relative(path, base):
    """Return ``path`` relative to ``base`` for sessions.

    Remote URLs pass through unchanged and empty paths stay empty; this
    matches the semantics previously duplicated in
    ``sequenceview.session_data``. Cross-filesystem cases where
    :func:`os.path.relpath` fails keep the original absolute path.
    """
    if not path or is_remote(path):
        return path
    try:
        return os.path.relpath(path, base)
    except ValueError:
        return path


def to_absolute(path, base):
    """Return the absolute, resolved form of a saved ``path`` under ``base``.

    Remote URLs pass through unchanged; relative saved paths are resolved
    against ``base`` (the directory containing the session file). This
    mirrors the previous ``sequenceview.restore`` behavior; callers keep
    their own existence checks and validation.
    """
    if is_remote(path):
        return path
    p = Path(path)
    p = (Path(base) / p).resolve() if not p.is_absolute() else p
    return str(p)


def cleanup_temp(directory):
    """Best-effort removal of an application-owned temporary directory.

    On Windows a file can stay briefly locked after the worker that owned it
    exits (pending handle release, antivirus scanning), so a single
    TemporaryDirectory.cleanup() may raise PermissionError. Retry with short
    backoff; leave the directory in place if it still will not go away rather
    than crashing at shutdown."""
    import shutil
    import time
    from pathlib import Path
    for delay in (0, .5, 1, 2, 4):
        if delay:
            time.sleep(delay)
        try:
            shutil.rmtree(directory, ignore_errors=False)
            return
        except FileNotFoundError:
            return
        except (PermissionError, OSError):
            continue
    try:
        shutil.rmtree(directory, ignore_errors=True)
    except OSError:
        pass
