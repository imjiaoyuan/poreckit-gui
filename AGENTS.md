# AGENTS.md

## What this is

Native PySide6 (Qt Widgets) desktop genome browser and plotting workbench for the sibling [`poreckit`](../poreckit) library (Python + Rust). One Python package (`src/poreckit_gui`, flat ~38 modules) with a single entry point `poreckit-gui = poreckit_gui.app:main`; runs on Windows/Linux/macOS and is frozen with PyInstaller into per-OS installers. Feature scope and user-facing behavior are owned by `README.md`; packaging, CI and release steps by `docs/packaging.md` (written in Chinese).

## Commands

Run from the repo root with `uv`; `--no-sync` skips re-resolving the lockfile.

```sh
uv sync                                            # create .venv; builds the sibling ../poreckit (needs Rust toolchain, Python 3.11+)
uv sync --reinstall-package poreckit               # rebuild poreckit after editing the sibling source
uv run poreckit-gui                                # run the app
uv run pytest -q                                   # full suite; headless — tests/conftest.py forces QT_QPA_PLATFORM=offscreen
uv run pytest -q tests/test_browser.py -k locus    # single file / single test
uv run --no-sync python scripts/write_catalog.py   # regenerate docs/plot-catalog.md from the installed poreckit
```

CI (`.github/workflows/release.yml`) runs `uv run --no-sync pytest -q -n auto --timeout=300` after installing `pytest-xdist` and `pytest-timeout`; those plugins are not in the local dev group. No linter or formatter is configured.

Build installers (per-OS only, no cross-compilation; details in `docs/packaging.md`):

```sh
uv sync --locked --group build
uv run --no-sync python scripts/build_release.py                     # freeze + package -> dist/releases/
uv run --no-sync python scripts/build_release.py --skip-freeze       # repackage only, reuse last freeze
uv run --no-sync python scripts/build_release.py --format deb        # deb | arch | dmg | "exe zip" (current OS only)
```

Frozen self-test: `QT_QPA_PLATFORM=offscreen dist/poreckit-gui/poreckit-gui --self-test report.json` (macOS: the `.app` binary path); check `ok` and `frozen` in the JSON.

## Architecture

- `app.py` — `Browser` main window and entry point.
- `model.py`, `viewport.py`, `states.py`, `sessionio.py` — portable session state, zero-based half-open locus navigation, history, persistence.
- `tracks.py`, `contactview.py`, `dataview.py`, `mapdata.py`, `tracksettings.py` — browser contact/annotation panels and bounded data queries.
- `axistracks.py`, `annotations.py`, `annotationui.py`, `genesearch.py`, `searchui.py`, `expected.py` — axis-side tracks, BEDPE layers, gene search, expected curves.
- `genome.py`/`genomeui.py`, `assembly.py`/`assemblyview.py`, `sequenceio.py`/`sequenceview.py` — whole-genome maps, scaffold assembly editing with coordinate lifting, signal/alignment/variant/reference viewer.
- `bigwig.py`, `remote.py` — local and HTTP byte-range signal readers with bounded caching.
- `workbench.py` — catalog UI, editable track stacks, shared figure viewer; backed by `catalog.py` (discovers poreckit commands/functions/track classes), `forms.py` (typed parameter forms), `plotdata.py` (declarative, non-executable data references), `execution.py` (runs plots through poreckit's original API).
- `jobs.py`, `worker.py`, `cache.py` — reusable cancellable subprocess transport and file-versioned reader/matrix caches; `demo.py` — synthetic demo data.
- `widgets.py`, `windows.py` — shared widgets and linked windows; `selftest.py` — packaged-app self-test.

Supporting trees: `tests/` (pytest against real poreckit rendering), `scripts/` (`build_release.py`, `check_release_tag.py`, `write_catalog.py`), `packaging/` (PyInstaller spec, Inno Setup `.iss`, `arch/PKGBUILD`, icons), `examples/differential-matrix.json` (executable preset), `.github/poreckit-ref` (pinned backend commit).

Storage: no persistent data inside the repo — runtime scratch lives in temp dirs, window settings in `QSettings("poreckit", "GenomeBrowser")`, sessions/presets/annotations are JSON files that reference data paths without copying them, and build outputs land in `dist/releases/`.

## Conventions

- Coordinates are zero-based, half-open everywhere, matching poreckit.
- Call poreckit's public reader/Canvas/track APIs and Click command objects directly; never shell out to poreckit CLI strings. Parsing and computation stay in poreckit's Rust extension.
- Preserve data contracts: workbench functions expect matrices or result dicts — never substitute a contact file for an analysis result; never load pickles or execute code; data references stay declarative (`$matrix`, `$reader`, `$contacts`, `$file`, `$array`).
- Respect the bounded-resource model: dense matrix queries reject more than 2500 bins per side; caches are file-versioned and bounded; plotting runs in cancellable workers; preview renders SVG only, PNG/PDF on demand.
- All application text in English; no custom stylesheet — follow the desktop Qt palette and style.
- `docs/plot-catalog.md` is generated by `scripts/write_catalog.py`; never hand-edit it.
- Backend coupling: `poreckit` is a regular (non-editable) install from `../poreckit`. When the backend changes, write its full 40-hex commit into `.github/poreckit-ref`, refresh `uv.lock`, and review `pyproject.toml`, `uv.lock`, `PKGBUILD` and the ref together.
- Releases: the tag must equal `v<version>` from `pyproject.toml` (enforced by `scripts/check_release_tag.py`); publishing requires all five CI targets to pass. Signing/notarization and verification details live in `docs/packaging.md`.
