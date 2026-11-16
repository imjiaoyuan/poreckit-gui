import json
import time

import numpy as np
import pytest

from poreckit_gui.contactview import transform_matrix
from poreckit_gui.demo import create_demo
from poreckit_gui.mapdata import MapData
from poreckit_gui.model import Region
from poreckit_gui.viewport import Viewport
from poreckit_gui.worker import render


def test_independent_axis_navigation():
    view = Viewport(Region("chr1", 100, 500), Region("chr2", 1000, 1800))
    zoomed = view.zoom(.5, .25, .75)
    assert zoomed.x == Region("chr1", 150, 350)
    assert zoomed.y == Region("chr2", 1300, 1700)
    assert view.pan(.5, -.5) == Viewport(Region("chr1", 300, 700), Region("chr2", 600, 1400))
    assert view.select(.25, .75, .5, 1).x == Region("chr1", 200, 400)


def test_control_formulas_and_zero_handling():
    a = np.array([[0, 2], [4, 8]], dtype=np.float32)
    b = np.array([[0, 1], [2, 4]], dtype=np.float32)
    np.testing.assert_allclose(transform_matrix(a, b, "difference"), a-b)
    np.testing.assert_allclose(transform_matrix(a, b, "ratio", 1), (a+1)/(b+1))
    np.testing.assert_allclose(transform_matrix(a, b, "log2ratio", 1), np.log2((a+1)/(b+1)))
    assert np.isnan(transform_matrix(a, b, "ratio", 0)[0, 0])
    split = transform_matrix(a, b, "split")
    assert split[0, 1] == 2 and split[1, 0] == 2


def test_cross_chromosome_payload_and_axis_annotations(tmp_path):
    from poreckit import make_reader
    path = create_demo(tmp_path, [("chr1", 1200000), ("chr2", 800000)])
    signal = tmp_path / "signal.bg"
    signal.write_text("chr1\t0\t1200000\t2\nchr2\t0\t800000\t8\n")
    loops = tmp_path / "loops.bedpe"
    loops.write_text("chr1\t200000\t210000\tchr2\t400000\t410000\ttrans\t5\n")
    view = dict(input=path, region="chr1:100001-899999", region_y="chr2:0-600000", mode="square",
                binsize=10000, cmap="hic_whitered", norm=None, log=False, vmax_percentile=95,
                tracks=[dict(kind="signal", path=str(signal), placement="auto"),
                        dict(kind="loops", path=str(loops), placement="auto")])
    result = render(view, str(tmp_path / "trans.svg"))
    assert result["region"] == "chr1:100,000-900,000"
    assert result["region_y"] == "chr2:0-600,000"
    data = MapData(result["map_data"])
    rd = make_reader(path)
    try:
        expected = rd.matrix("chr2", 0, 600000, 10000, None, chrom2="chr1", start2=100000, end2=900000)
    finally:
        rd.close()
    np.testing.assert_allclose(data.observed, expected)
    assert expected.shape == (60, 80)
    assert result["map_data"]["bounds"][2] > 0
    hover = data.at(.25, .5)
    assert hover["x"] == ("chr1", 300000, 310000)
    assert hover["y"] == ("chr2", 300000, 310000)
    assert hover["observed"] == expected[30, 20]


def test_control_render_and_incompatible_genome(tmp_path):
    import h5py
    observed = create_demo(tmp_path / "a")
    control = create_demo(tmp_path / "b", scale=2)
    view = dict(input=observed, control=control, region="chr1:0-1000000", mode="square", binsize=10000,
                cmap="RdBu_r", norm=None, log=False, vmax_percentile=95, tracks=[], display_mode="difference")
    result = render(view, str(tmp_path / "difference.svg"))
    data = MapData(result["map_data"])
    np.testing.assert_allclose(data.displayed, -data.observed)
    assert result["display"]["color_min"] < 0 < result["display"]["color_max"]
    with h5py.File(control, "r+") as f:
        f.attrs["assembly"] = "another-build"
    with pytest.raises(ValueError, match="different genome"):
        render(view, str(tmp_path / "invalid.svg"))


def test_gui_data_drag_zoom_hover_and_session(tmp_path, monkeypatch):
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.show()
    browser.resize(1400, 1000)
    browser.auto_update.setChecked(False)
    errors = []
    browser.error = lambda message: errors.append(str(message))

    def wait():
        deadline = time.monotonic() + 30
        while browser.busy() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.01)
        assert not browser.busy(), browser.log_view.toPlainText()
        assert not errors, errors

    try:
        path = create_demo(tmp_path, [("chr1", 1200000), ("chr2", 800000)])
        browser.load_file(path)
        wait()
        browser.locus.setText("chr1:0-600000")
        browser.independent_y.setChecked(True)
        browser.locus_y.setText("chr2:0-600000")
        browser.render()
        wait()
        assert browser.image.map_data.observed.shape == (60, 60)
        centre = browser.image.mapFromScene(browser.image.data_rect.center())
        QTest.mouseMove(browser.image.viewport(), centre, delay=10)
        assert browser.hovered_bin["x"][0] == "chr1"
        assert browser.hovered_bin["y"][0] == "chr2"
        QTest.mouseDClick(browser.image.viewport(), Qt.MouseButton.LeftButton, pos=centre)
        wait()
        assert 300000 <= browser.current_viewport().x.end-browser.current_viewport().x.start <= 320000
        before = browser.current_viewport()
        centre = browser.image.mapFromScene(browser.image.data_rect.center())
        QTest.mousePress(browser.image.viewport(), Qt.MouseButton.LeftButton, pos=centre)
        QTest.mouseMove(browser.image.viewport(), centre-QPoint(60, 0), delay=20)
        QTest.mouseRelease(browser.image.viewport(), Qt.MouseButton.LeftButton, pos=centre-QPoint(60, 0))
        wait()
        assert browser.current_viewport().x.start > before.x.start
        assert browser.current_viewport().y == before.y
        session = tmp_path / "xy-session.json"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(session), "JSON (*.json)"))
        browser.save_session()
        data = json.loads(session.read_text())
        assert data["region_y"].startswith("chr2:")
        browser.independent_y.setChecked(False)
        browser.render()
        wait()
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(session), "JSON (*.json)"))
        browser.load_session()
        wait()
        assert browser.independent_y.isChecked()
        assert browser.current_viewport().y.chrom == "chr2"
    finally:
        browser.close()
