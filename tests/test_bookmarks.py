"""Bookmark list, JSON roundtrip, QSettings persistence and slideshow state machine."""
import json
import time

import pytest

from poreckit_gui.demo import create_demo


@pytest.fixture
def browser_env(tmp_path, monkeypatch):
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication, QMessageBox

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    # Never let a modal dialog hang the offscreen test run.
    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.StandardButton.No))
    settings = QSettings(str(tmp_path / "gui.ini"), QSettings.Format.IniFormat)
    browser = Browser(settings=settings)
    browser.show()
    browser.auto_update.setChecked(False)
    errors = []
    browser.error = lambda value: errors.append(str(value))

    def wait():
        until = time.monotonic()+30
        while browser.busy() and time.monotonic() < until:
            app.processEvents()
            time.sleep(.01)
        assert not browser.busy(), browser.log_view.toPlainText()
        app.processEvents()

    browser.load_file(create_demo(tmp_path / "data"))
    wait()
    yield browser, errors, wait
    browser.annotations.model.mark_saved()
    browser.bookmarks.mark_saved()
    browser.close()
    settings.sync()


def _add(browser, monkeypatch, name):
    from PySide6.QtWidgets import QInputDialog
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: (name, True))
    browser.bookmarks.add()


def _set_file_dialog(monkeypatch, path):
    from PySide6.QtWidgets import QFileDialog
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(path), ""))


def test_bookmark_menu_wiring():
    from PySide6.QtWidgets import QApplication

    from poreckit_gui.app import Browser
    QApplication.instance() or QApplication([])
    browser = Browser()
    try:
        menus = {action.text(): action.menu() for action in browser.menuBar().actions()}
        menu = [item.text() for item in menus["Bookmarks"].actions()]
        for title in ("Bookmark current view…", "Import bookmarks…", "Export bookmarks…"):
            assert title in menu
        assert "Bookmarks" in menu
        for text in ("Add…", "Update", "Rename…", "Remove", "Up", "Down",
                     "Import…", "Export…", "Next"):
            assert any(button.text() == text for button in browser.bookmarks.dock.findChildren(type(browser.bookmarks.play)))
    finally:
        browser.close()


def test_bookmark_add_jump_json_roundtrip_and_missing_skip(browser_env, tmp_path, monkeypatch):
    browser, errors, wait = browser_env
    controller = browser.bookmarks
    _add(browser, monkeypatch, "Full")
    wait()
    browser.locus.setText("chr1:0-1000000")
    browser.render()
    wait()
    _add(browser, monkeypatch, "Zoom")
    wait()
    assert [entry["name"] for entry in controller.entries] == ["Full", "Zoom"]
    # Jump roundtrip through the state object.
    controller.list.setCurrentRow(1)
    controller.open_current()
    wait()
    assert browser.current_viewport().x.end == 1000000
    controller.list.setCurrentRow(0)
    controller.open_current()
    wait()
    assert browser.current_viewport().x.end == 2000000
    # Rename, reorder, delete.
    monkeypatch.setattr("PySide6.QtWidgets.QInputDialog.getText",
                        lambda *a, **k: ("Renamed", True))
    controller.list.setCurrentRow(1)
    controller.rename()
    assert controller.entries[1]["name"] == "Renamed"
    controller.move(-1)
    assert [entry["name"] for entry in controller.entries] == ["Renamed", "Full"]
    controller.list.setCurrentRow(0)
    controller.remove()
    assert [entry["name"] for entry in controller.entries] == ["Full"]
    # Rebuild two bookmarks, then export to the session JSON dialect.
    _add(browser, monkeypatch, "Second")
    wait()
    path = tmp_path / "bookmarks.json"
    _set_file_dialog(monkeypatch, path)
    controller.export()
    data = json.loads(path.read_text())
    assert data["version"] == 1 and len(data["bookmarks"]) == 2
    assert Path_is_relative(data["bookmarks"][0]["view"]["input"], tmp_path)
    # Missing-file bookmarks are skipped with a clear error, valid ones import.
    data["bookmarks"][0]["view"]["input"] = str(tmp_path / "gone.hic")
    path.write_text(json.dumps(data))
    controller.entries.clear()
    controller.refresh()
    controller.mark_saved()
    controller.load()
    wait()
    assert [entry["name"] for entry in controller.entries] == ["Second"]
    assert errors and "Second" not in errors[0] and "skipped" in errors[0]
    # A clean import appends and marks the controller saved.
    controller.entries.clear()
    controller.refresh()
    data["bookmarks"][0]["view"]["input"] = data["bookmarks"][1]["view"]["input"]
    path.write_text(json.dumps(data))
    errors.clear()
    controller.load()
    wait()
    assert len(controller.entries) == 2 and not errors


def Path_is_relative(path, base):
    from pathlib import Path
    path = Path(path)
    if not path.is_absolute():
        path = Path(base) / path
    return path.resolve().is_relative_to(Path(base).resolve())


def test_bookmark_qsettings_persistence(browser_env, tmp_path, monkeypatch):
    from PySide6.QtCore import QSettings

    from poreckit_gui.app import Browser
    browser, errors, wait = browser_env
    _add(browser, monkeypatch, "Kept")
    wait()
    browser.locus.setText("chr1:0-500000")
    browser.render()
    wait()
    _add(browser, monkeypatch, "Also kept")
    wait()
    assert browser.bookmarks.entries
    browser.bookmarks.persist()
    browser.settings.sync()
    settings = QSettings(str(tmp_path / "gui.ini"), QSettings.Format.IniFormat)
    stored = json.loads(settings.value("bookmarks"))
    assert [entry["name"] for entry in stored] == [entry["name"] for entry in browser.bookmarks.entries]
    other = Browser(settings=settings)
    try:
        assert [entry["name"] for entry in other.bookmarks.entries] == \
               [entry["name"] for entry in browser.bookmarks.entries]
        assert not other.bookmarks.dirty
    finally:
        other.close()


def test_slideshow_state_machine_no_wraparound(browser_env, monkeypatch):
    browser, errors, wait = browser_env
    controller = browser.bookmarks
    for name in ("One", "Two", "Three"):
        _add(browser, monkeypatch, name)
        wait()
    controller.list.setCurrentRow(0)
    controller.play.setChecked(True)
    wait()
    assert controller.play.isChecked()
    assert controller.indicator.text() == "1/3 — One"
    assert controller.timer.isActive()
    controller.advance()
    wait()
    assert controller.list.currentRow() == 1
    assert controller.indicator.text() == "2/3 — Two"
    controller.advance()
    wait()
    assert controller.list.currentRow() == 2
    assert controller.indicator.text() == "3/3 — Three"
    # No wraparound: the last advance stops the slideshow at the end.
    controller.advance()
    wait()
    assert not controller.play.isChecked()
    assert not controller.timer.isActive()
    assert not controller.indicator.isVisible()
    assert controller.play.text() == "Start slideshow"
    # Escape exits a running slideshow.
    controller.play.setChecked(True)
    wait()
    controller.escape.activated.emit()
    assert not controller.play.isChecked() and not controller.timer.isActive()
