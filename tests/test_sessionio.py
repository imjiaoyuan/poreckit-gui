"""Regression tests for the shared session persistence helpers."""

import pytest

from poreckit_gui.sessionio import load_json, save_json, to_absolute, to_relative


def test_save_then_load_roundtrip(tmp_path):
    data = {"version": 1, "tracks": [{"path": "a.cool", "nested": {"é": "ü"}}], "counts": [1, 2, 3]}
    path = tmp_path / "session.json"
    save_json(path, data)
    assert load_json(path) == data


def test_load_json_reports_invalid_json_with_path(tmp_path):
    path = tmp_path / "broken.json"
    path.write_text('{"version": 1,', encoding="utf-8")
    with pytest.raises(ValueError) as excinfo:
        load_json(path)
    assert "broken.json" in str(excinfo.value)


def test_to_relative_and_absolute_roundtrip(tmp_path):
    target = tmp_path / "data" / "map.cool"
    saved = to_relative(str(target), tmp_path)
    assert saved == os_relpath(target, tmp_path)
    assert to_absolute(saved, tmp_path) == str(target.resolve())
    # Remote URLs and empty paths pass through untouched.
    assert to_relative("https://example.com/map.cool", tmp_path) == "https://example.com/map.cool"
    assert to_absolute("https://example.com/map.cool", tmp_path) == "https://example.com/map.cool"
    assert to_relative("", tmp_path) == ""


def os_relpath(path, base):
    import os
    return os.path.relpath(str(path), str(base))


def test_interrupted_save_keeps_previous_session(tmp_path, monkeypatch):
    path = tmp_path / "session.json"
    save_json(path, {"version": 1, "note": "previous"})
    previous = path.read_text()

    from poreckit_gui import sessionio

    def boom(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(sessionio.os, "fsync", boom)
    with pytest.raises(OSError):
        save_json(path, {"version": 1, "note": "interrupted"})
    # The previous session survives byte for byte and no temp files remain.
    assert path.read_text() == previous
    assert list(tmp_path.iterdir()) == [path]


def test_failing_serialization_keeps_previous_session(tmp_path, monkeypatch):
    path = tmp_path / "session.json"
    save_json(path, {"version": 1})
    previous = path.read_text()

    import json as jsonlib

    def bad_dumps(*args, **kwargs):
        raise TypeError("not serializable")

    monkeypatch.setattr(jsonlib, "dumps", bad_dumps)
    with pytest.raises(TypeError):
        save_json(path, {"bad": object()})
    assert path.read_text() == previous
    assert list(tmp_path.iterdir()) == [path]
    assert load_json(path) == {"version": 1}


def test_sequence_view_save_reports_error_and_keeps_old_file(tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui import sessionio
    from poreckit_gui.app import Browser

    _app = QApplication.instance() or QApplication([])
    browser = Browser()
    seq = browser.sequence_browser
    try:
        path = tmp_path / "tracks.json"
        save_json(path, {"type": "poreckit-sequence", "version": 1})
        previous = path.read_text()
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(path), "JSON (*.json)"))

        def boom(target, data):
            raise OSError("disk full")

        monkeypatch.setattr(sessionio, "save_json", boom)
        seq.save_session()
        assert seq.status.text() == "disk full"
        assert path.read_text() == previous
        assert not browser.sequence_browser.current_session
    finally:
        browser.close()
