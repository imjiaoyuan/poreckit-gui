import json
import time

import pytest

from poreckit_gui.annotations import AnnotationModel, read_bedpe, write_bedpe, feature_payload
from poreckit_gui.model import Region
from poreckit_gui.demo import create_demo


def test_annotation_roundtrip_edits_and_mirrored_hits(tmp_path):
    path = tmp_path / "loops.bedpe"
    path.write_text("chr1\t100\t200\tchr2\t300\t400\tloop\t7\t+\t-\n")
    model = AnnotationModel()
    key = model.add_layer("Loops", read_bedpe(path))
    model.mark_saved()
    value = model.layer(key)["features"][0]
    fid = value["id"]
    model.update(key, fid, x0=120, x1=220)
    assert model.dirty
    model.undo()
    assert value["x0"] == 100 and not model.dirty
    model.redo()
    assert value["x0"] == 120
    with pytest.raises(ValueError):
        model.update(key, fid, x1=0)
    assert value["x1"] == 220
    hit = feature_payload(model.layers, Region("2", 0, 1000), Region("1", 0, 1000))[0]
    assert hit["reflected"] and hit["coords"] == (300, 400, 120, 220)
    model.set_layer(key, visible=False)
    assert not feature_payload(model.layers, Region("2", 0, 1000), Region("1", 0, 1000))
    export = tmp_path / "saved.bedpe"
    write_bedpe(export, model.layer(key))
    restored = read_bedpe(export)[0]
    assert restored["extra"] == ["7", "+", "-"]
    assert restored["x0"] == 120
    model.remove_layer(key)
    model.undo()
    assert model.layer(key)["features"][0]["id"] == fid


def test_gui_annotation_gestures_and_bookmark_roundtrip(tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication, QFileDialog, QInputDialog
    from PySide6.QtCore import Qt, QPoint
    from PySide6.QtTest import QTest
    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.resize(1400, 1000)
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
        assert not errors, errors
        app.processEvents()

    try:
        browser.load_file(create_demo(tmp_path / "data"))
        wait()
        controller = browser.annotations
        key = controller.model.add_layer("Manual")
        fid = controller.model.add(key, dict(xchrom="chr1", x0=400000, x1=600000,
                                            ychrom="chr1", y0=800000, y1=1000000, name="Loop"))
        browser.render()
        wait()
        image = browser.image
        point = image.mapFromScene(image.data_rect.x()+image.data_rect.width()*.25,
                                   image.data_rect.y()+image.data_rect.height()*.45)
        QTest.mouseClick(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.ControlModifier, point)
        app.processEvents()
        assert controller.dock.isVisible()
        assert controller.table.item(0, 6).text() == "Loop"
        # Dock geometry has settled before dragging the selected rectangle.
        point = image.mapFromScene(image.data_rect.x()+image.data_rect.width()*.25,
                                   image.data_rect.y()+image.data_rect.height()*.45)
        target = point+QPoint(20, 10)
        QTest.mousePress(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.ControlModifier, point)
        QTest.mouseMove(image.viewport(), target, delay=10)
        QTest.mouseRelease(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.ControlModifier, target)
        value = controller.model.layer(key)["features"][0]
        assert value["x0"] > 400000 and value["y0"] > 800000
        controller.model.undo()
        assert value["x0"] == 400000
        modifiers = Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier
        QTest.mousePress(image.viewport(), Qt.MouseButton.LeftButton, modifiers, point)
        QTest.mouseMove(image.viewport(), target, delay=10)
        QTest.mouseRelease(image.viewport(), Qt.MouseButton.LeftButton, modifiers, target)
        assert value["x0"] == 400000 and value["x1"] > 600000 and value["y1"] > 1000000
        controller.model.undo()
        assert value["x1"] == 600000
        monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Annotated", True))
        browser.bookmarks.add()
        browser.locus.setText("chr1:0-1000000")
        browser.render()
        wait()
        monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Detail", True))
        browser.bookmarks.add()
        browser.bookmarks.play.setChecked(True)
        wait()
        browser.bookmarks.advance()
        wait()
        browser.bookmarks.play.setChecked(False)
        assert browser.current_viewport().x.end == 2000000
        browser.bookmarks.list.setCurrentRow(0)
        browser.bookmarks.open_current()
        wait()
        assert browser.current_viewport().x.end == 2000000
        assert controller.model.layer(key)["features"][0]["id"] == fid
        path = tmp_path / "session.json"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(path), ""))
        browser.save_session()
        data = json.loads(path.read_text())
        assert data["bookmarks"][0]["name"] == "Annotated"
        assert not browser.bookmarks.dirty
        browser.bookmarks.entries.clear()
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(path), ""))
        browser.load_session()
        wait()
        assert len(browser.bookmarks.entries) == 2
        assert len(browser.image.map_data.spec["features"]) == 2
    finally:
        controller.model.mark_saved()
        browser.bookmarks.mark_saved()
        browser.close()
