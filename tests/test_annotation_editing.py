import time

from poreckit_gui.annotations import AnnotationModel, draw_layers, hit_test, read_bedpe, write_bedpe
from poreckit_gui.demo import create_demo
from poreckit_gui.model import Region


def model_with_feature():
    model = AnnotationModel()
    key = model.add_layer("Manual", [dict(xchrom="chr1", x0=400000, x1=600000,
                                          ychrom="chr1", y0=800000, y1=1000000, name="Loop", score=12.5)])
    return model, key


def test_editable_layer_state_and_feature_defaults():
    model, key = model_with_feature()
    assert model.layer(key)["unsaved"]
    model.mark_saved(key)
    assert not model.layer(key)["unsaved"] and not model.dirty
    fid = model.layer(key)["features"][0]["id"]
    model.update(key, fid, name="Renamed")
    assert model.layer(key)["unsaved"] and model.dirty
    model.undo()
    assert not model.layer(key)["unsaved"]
    model.redo()
    assert model.layer(key)["unsaved"]


def test_hit_test_prefers_nearest_anchor_and_interior():
    model, key = model_with_feature()
    layers = model.layers
    x = Region("chr1", 0, 2000000)
    inside = hit_test(layers, x, x, 500000, 900000, 4000, 4000)
    assert inside and inside["id"] == model.layer(key)["features"][0]["id"]
    # Loop anchor centre on the diagonal, just outside the rectangle.
    anchor = hit_test(layers, x, x, 500000, 500000, 4000, 4000)
    assert anchor and anchor["coords"] == (400000, 600000, 800000, 1000000)
    # A few pixels further the anchor is out of reach.
    assert hit_test(layers, x, x, 500000, 503000, 1000, 1000) is None
    # A second closer feature wins.
    model.add(key, dict(xchrom="chr1", x0=700000, x1=720000, ychrom="chr1", y0=700000, y1=720000))
    near = hit_test(layers, x, x, 710000, 710000, 4000, 4000)
    assert near["coords"] == (700000, 720000, 700000, 720000)


def test_delete_feature_and_undo_redo(tmp_path):
    model, key = model_with_feature()
    fid = model.layer(key)["features"][0]["id"]
    model.remove(key, fid)
    assert not model.layer(key)["features"]
    model.undo()
    assert model.layer(key)["features"][0]["id"] == fid
    model.redo()
    assert not model.layer(key)["features"]
    path = tmp_path / "out.bedpe"
    write_bedpe(path, model.layer(key))
    assert read_bedpe(path) == []


def test_bedpe_roundtrip_preserves_score_and_edits(tmp_path):
    path = tmp_path / "loops.bedpe"
    path.write_text("chr1\t100\t200\tchr1\t300\t400\tloopA\t7\t+\n"
                    "chr2\t10\t20\tchr2\t30\t40\ttadB\n")
    model = AnnotationModel()
    key = model.add_layer("Loaded", read_bedpe(path))
    fid = model.layer(key)["features"][0]["id"]
    model.update(key, fid, score=3, x1=250, name="loopA-edit")
    export = tmp_path / "saved.bedpe"
    write_bedpe(export, model.layer(key))
    rows = read_bedpe(export)
    assert rows[0]["name"] == "loopA-edit" and rows[0]["score"] == "3"
    assert (rows[0]["xchrom"], rows[0]["x0"], rows[0]["x1"]) == ("chr1", 100, 250)
    assert rows[0]["extra"] == ["+"]
    assert rows[1]["name"] == "tadB" and rows[1]["score"] == ""
    assert (rows[1]["xchrom"], rows[1]["x0"], rows[1]["x1"], rows[1]["y0"], rows[1]["y1"]) == ("chr2", 10, 20, 30, 40)


def test_unsaved_layers_draw_dashed_and_saved_solid():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    model, key = model_with_feature()
    x = Region("chr1", 0, 2000000)
    figure = plt.figure()
    ax = figure.add_subplot()
    draw_layers(ax, model.layers, x, x)
    assert ax.patches[0].get_linestyle() == (0, (4, 3))
    model.mark_saved(key)
    ax.clear()
    draw_layers(ax, model.layers, x, x)
    assert ax.patches[0].get_linestyle() == "solid"
    plt.close(figure)


def test_draw_mode_toggle_leaves_plain_click_alone(tmp_path, monkeypatch):
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication, QInputDialog

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

    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("", True))
    try:
        browser.load_file(create_demo(tmp_path / "data"))
        wait()
        controller = browser.annotations
        assert controller.model.layers == []
        image = browser.image
        # Draw mode off: a plain drag pans, no annotation is created.
        start = image.mapFromScene(image.data_rect.center())
        QTest.mousePress(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start)
        QTest.mouseMove(image.viewport(), start+QPoint(40, 30), delay=10)
        QTest.mouseRelease(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start+QPoint(40, 30))
        wait()
        assert controller.model.layers == []
        assert not errors
        # Draw mode on: the same plain drag creates a layer with one feature.
        controller.draw.trigger()
        assert image.annotation_mode
        start = image.mapFromScene(image.data_rect.center())
        QTest.mousePress(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start)
        QTest.mouseMove(image.viewport(), start+QPoint(-40, -30), delay=10)
        QTest.mouseRelease(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start+QPoint(-40, -30))
        wait()
        assert len(controller.model.layers) == 1
        features = controller.model.layers[0]["features"]
        assert len(features) == 1 and controller.model.layers[0]["unsaved"]
        # Toggle off restores plain navigation clicks.
        controller.draw.trigger()
        assert not image.annotation_mode
        viewport = browser.current_viewport()
        QTest.mouseClick(image.viewport(), Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, start)
        assert browser.current_viewport() == viewport
        assert len(controller.model.layers[0]["features"]) == 1
    finally:
        controller.model.mark_saved()
        browser.close()
