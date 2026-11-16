import json
import time
from pathlib import Path

import pytest

from poreckit_gui.demo import create_demo
from poreckit_gui.model import Region, validate_view
from poreckit_gui.worker import inspect_file, render


def test_navigation_bounds_and_coordinate_convention():
    region = Region.parse("chr1:1,000-5,000")
    assert region.start == 1000 and region.end == 5000
    assert region.pan(-0.5, 10000) == Region("chr1", 0, 4000)
    assert region.pan(10, 10000) == Region("chr1", 6000, 10000)
    assert region.zoom(100, 10000) == Region("chr1", 0, 10000)
    assert region.zoom(0.5, 10000) == Region("chr1", 2000, 4000)
    for text in ("chr1:5-1", "chr1:0-0", "chr1:-1-5", "chr1"):
        with pytest.raises(ValueError):
            Region.parse(text)


@pytest.fixture
def view(tmp_path):
    path = create_demo(tmp_path)
    return {"input": path, "region": "chr1:0-2,000,000", "binsize": None,
            "mode": "triangle", "norm": None, "cmap": "Reds", "log": True, "tracks": []}


def test_metadata_and_large_matrix_guard(view):
    meta = inspect_file(view["input"])
    assert meta["chromsizes"] == [("chr1", 2000000)]
    assert meta["resolutions"] == [10000]
    with pytest.raises(ValueError, match="2500"):
        validate_view({**view, "binsize": 1})


def test_render_real_core_with_tracks_and_vector_export(view, tmp_path):
    genes = tmp_path / "genes.gtf"
    genes.write_text('chr1\tdemo\texon\t200001\t400000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "DEMO1";\n')
    signal = tmp_path / "signal.bedgraph"
    signal.write_text("chr1\t0\t1000000\t2\nchr1\t1000000\t2000000\t8\n")
    tads = tmp_path / "tads.bed"
    tads.write_text("chr1\t0\t400000\tdomain1\n")
    # Overlay order must remain correct even if a TAD follows a signal in the list.
    view["tracks"] = [{"kind": kind, "path": str(path), "visible": True}
                      for kind, path in [("genes", genes), ("signal", signal), ("tads", tads)]]
    result = render(view, str(tmp_path / "triangle.png"))
    assert result["binsize"] == 10000
    assert (tmp_path / "triangle.png").read_bytes().startswith(b"\x89PNG")
    view["mode"] = "square"
    render(view, str(tmp_path / "square.svg"))
    assert "<svg" in (tmp_path / "square.svg").read_text()
    render(view, str(tmp_path / "square.pdf"))
    assert (tmp_path / "square.pdf").read_bytes().startswith(b"%PDF")


def test_gui_worker_navigation_and_failed_load_preserves_view(view, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui.app import Browser

    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.auto_update.setChecked(False)
    errors = []
    browser.error = lambda message: errors.append(str(message))
    browser.show()

    def wait():
        deadline = time.monotonic() + 45
        while browser.busy() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        assert not browser.busy(), browser.log_view.toPlainText()

    try:
        browser.load_file(view["input"])
        wait()
        assert not errors, browser.log_view.toPlainText()
        assert browser.last_view is not None
        assert browser.stack.currentWidget() is browser.image
        browser.zoom(0.5)
        wait()
        assert browser.region == Region("chr1", 500000, 1500000)
        browser.travel(-1)
        wait()
        assert browser.region == Region("chr1", 0, 2000000)
        browser.mode.setCurrentIndex(1)
        browser.render()
        wait()
        assert browser.last_view["mode"] == "square"
        svg = tmp_path / "snapshot.svg"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(svg), "SVG (*.svg)"))
        browser.cmap.setCurrentText("viridis")
        browser.export()
        assert not browser.busy()
        assert svg.read_bytes() == Path(browser.last_files["svg"]).read_bytes()
        browser.cmap.setCurrentText("Reds")
        export_path = tmp_path / "export.pdf"
        browser.render(output=str(export_path))
        wait()
        assert export_path.read_bytes().startswith(b"%PDF")
        session_path = tmp_path / "session.json"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(session_path), "JSON (*.json)"))
        browser.save_session()
        saved = json.loads(session_path.read_text())
        assert saved["input"] == Path(view["input"]).name
        browser.zoom(0.5)
        wait()
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **kw: (str(session_path), "JSON (*.json)"))
        browser.load_session()
        wait()
        assert browser.region == Region("chr1", 0, 2000000)
        assert browser.last_view["mode"] == "square"
        browser.render()
        browser.cancel()
        wait()
        assert browser.statusBar().currentMessage() == "Task cancelled"
        assert browser.last_view["mode"] == "square"
        browser.load_file(str(tmp_path / "missing.cool"))
        wait()
        assert errors
        assert browser.input_path == view["input"]
        assert browser.last_view["mode"] == "square"
    finally:
        browser.close()


def test_signal_panel_preserves_source_scale_and_negative_values(tmp_path):
    import matplotlib.pyplot as plt

    from poreckit_gui.tracks import SignalPanel

    path = tmp_path / "signal.bg"
    path.write_text("chr1\t0\t1000\t-2\nchr1\t1000\t2000\t8\n")
    fig, ax = plt.subplots()
    try:
        SignalPanel(path=str(path), chrom="chr1").draw(ax, "chr1", 0, 2000)
        values = ax.lines[0].get_ydata()
        assert min(values) == -2 and max(values) == 8
        assert ax.get_ylim()[0] < -2 and ax.get_ylim()[1] > 8
    finally:
        plt.close(fig)


def test_square_map_geometry_palette_limits_and_raw_values(view, tmp_path, monkeypatch):
    import numpy as np
    from matplotlib.figure import Figure
    from poreckit import make_reader

    reader = make_reader(view["input"])
    expected = reader.matrix("chr1", 0, 2000000, 10000, None)
    reader.close()
    original_save = Figure.savefig
    observed = []

    def inspect_figure(fig, *args, **kwargs):
        result = original_save(fig, *args, **kwargs)
        ax = fig.axes[0]
        image = ax.images[0]
        observed.append((np.asarray(image.get_array()).copy(), image.get_cmap().name,
                         image.get_clim(), ax.get_window_extent().width / ax.get_window_extent().height,
                         ax.get_ylim()))
        return result

    monkeypatch.setattr(Figure, "savefig", inspect_figure)
    settings = {**view, "mode": "square", "log": False, "cmap": "hic_whitered", "vmax": 50}
    result = render(settings, str(tmp_path / "square.png"))
    pixels, palette, limits, ratio, ylim = observed[-1]
    np.testing.assert_allclose(pixels, expected)
    assert palette == "hic_whitered"
    assert limits == (0, 50)
    assert ratio == pytest.approx(1)
    assert ylim == (2000000, 0)
    assert result["display"]["color_max"] == 50
    settings.update(log=True, cmap="RdBu_r", vmax=None, vmax_percentile=95)
    result = render(settings, str(tmp_path / "log.png"))
    pixels, palette, limits, ratio, ylim = observed[-1]
    np.testing.assert_allclose(pixels, np.log1p(expected))
    positive = expected[np.triu_indices_from(expected, k=1)]
    positive = positive[np.isfinite(positive) & (positive > 0)]
    upper = float(np.percentile(positive, 95))
    assert palette == "RdBu_r"
    assert limits[1] == pytest.approx(np.log1p(upper))
    assert result["display"]["color_max"] == pytest.approx(upper)


def test_display_settings_auto_apply_and_manual_update(view):
    from PySide6.QtWidgets import QApplication

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    errors = []
    browser.error = lambda message: errors.append(str(message))
    browser.show()

    def wait_until(predicate):
        deadline = time.monotonic() + 30
        while not predicate() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.01)
        assert predicate(), browser.log_view.toPlainText()
        assert not errors, errors

    try:
        browser.load_file(view["input"])
        wait_until(lambda: browser.last_view is not None and not browser.busy())
        assert browser.last_view["mode"] == "square"
        assert browser.last_view["log"] is False
        browser.mode.setCurrentIndex(browser.mode.findData("triangle"))
        browser.cmap.setCurrentText("RdBu_r")
        browser.log.setChecked(True)
        browser.color_range.setCurrentIndex(browser.color_range.findData("manual"))
        browser.color_max.setValue(25)
        wait_until(lambda: browser.last_view["mode"] == "triangle" and not browser.busy())
        assert browser.last_view["cmap"] == "RdBu_r"
        assert browser.last_view["log"] is True
        assert browser.last_view["vmax"] == 25
        assert "color max 25" in browser.view_label.text()
        browser.auto_update.setChecked(False)
        browser.mode.setCurrentIndex(browser.mode.findData("square"))
        assert not browser.display_timer.isActive()
        assert browser.last_view["mode"] == "triangle"
        browser.render()
        wait_until(lambda: browser.last_view["mode"] == "square" and not browser.busy())
    finally:
        browser.close()


def test_invalid_session_is_rejected_before_state_changes(view, tmp_path):
    """Regression: a session missing required keys once passed validation and
    then raised KeyError inside loaded(), leaving loading_controls stuck on
    so automatic updates were silently ignored."""
    from PySide6.QtWidgets import QApplication

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.show()
    errors = []
    browser.error = lambda value: errors.append(str(value))

    def wait_until(predicate):
        deadline = time.monotonic() + 30
        while not predicate() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.01)
        assert predicate(), browser.log_view.toPlainText()

    try:
        browser.load_file(view["input"])
        wait_until(lambda: browser.last_view is not None and not browser.busy())
        before = dict(browser.view_data())
        history_before = list(browser.history)
        # 'binsize' is missing: validate_view treats it as auto, but the
        # session loader indexes it directly.
        broken = {**before, "version": 1, "input": view["input"]}
        broken.pop("binsize")
        with pytest.raises(ValueError, match="binsize"):
            browser.validate_session(broken)
        with pytest.raises(ValueError, match="binsize"):
            browser.apply_session(broken)
        assert browser.input_path == view["input"]
        assert browser.history == history_before
        assert browser.loading_controls is False
        wait_until(lambda: not browser.busy())
        assert browser.view_data()["region"] == before["region"]
        # A failure while applying an already-validated session must not
        # leave the flag stuck either.
        genes = tmp_path / "genes.gtf"
        genes.write_text('chr1\tdemo\tgene\t101\t500\t.\t+\t.\tgene_id "g1";\n')
        session = {**before, "version": 1, "input": view["input"],
                   "region": "chr1:0-2,000,000",
                   "tracks": [{"kind": "genes", "path": str(genes), "visible": True}]}
        browser.pending_session = session
        browser.append_track = lambda track: (_ for _ in ()).throw(RuntimeError("boom"))
        try:
            browser.loaded(browser.metadata)
        except RuntimeError:
            pass
        else:
            raise AssertionError("expected the injected failure")
        assert browser.loading_controls is False
        wait_until(lambda: not browser.busy())
    finally:
        browser.close()


def test_session_round_trip_restores_view_and_cleans_temp_dir(tmp_path, monkeypatch):
    """Saving a session, closing the browser and restoring it in a new
    Browser must reinstate region and display settings, render again, and
    clean up each browser instance's temporary directory."""
    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.auto_update.setChecked(False)
    browser.show()
    errors = []
    browser.error = lambda message: errors.append(str(message))

    def wait_until(predicate, label=""):
        deadline = time.monotonic() + 45
        while not predicate() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        assert predicate(), f"{label}: {browser.log_view.toPlainText()}"
        assert not errors, errors

    try:
        browser.load_file(create_demo(tmp_path, [("chr1", 2000000), ("chr2", 1500000)]))
        wait_until(lambda: browser.last_view is not None and not browser.busy(), "load")
        browser.mode.setCurrentIndex(browser.mode.findData("triangle"))
        browser.cmap.setCurrentText("viridis")
        browser.log.setChecked(True)
        browser.navigate(Region.parse("chr2:250,000-1,250,000"))
        wait_until(lambda: browser.region == Region("chr2", 250000, 1250000)
                   and browser.last_view["mode"] == "triangle"
                   and browser.last_view["cmap"] == "viridis"
                   and browser.last_view["log"] is True and not browser.busy(), "settings")
        session_path = tmp_path / "session.json"
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            lambda *a, **kw: (str(session_path), "JSON (*.json)"))
        browser.save_session()
        assert session_path.is_file()
        temp_dir = Path(browser.temp.name)
        browser.close()
        assert not temp_dir.exists()

        restored = Browser()
        restored.auto_update.setChecked(False)
        errors = []
        restored.error = lambda message: errors.append(str(message))
        try:
            monkeypatch.setattr(QFileDialog, "getOpenFileName",
                                lambda *a, **kw: (str(session_path), "JSON (*.json)"))
            restored.load_session()
            wait_until(lambda: restored.last_view is not None and not restored.busy(), "restore")
            assert restored.input_path == browser.input_path
            assert restored.region == Region("chr2", 250000, 1250000)
            assert restored.last_view["mode"] == "triangle"
            assert restored.last_view["cmap"] == "viridis"
            assert restored.last_view["log"] is True
            # A fresh render must succeed in the restored session.
            out = tmp_path / "restored.svg"
            restored.render(output=str(out))
            wait_until(lambda: out.is_file() and out.stat().st_size > 1000
                       and not restored.busy(), "re-render")
            assert "<svg" in out.read_text()
        finally:
            restored_temp = Path(restored.temp.name)
            restored.close()
        assert not restored_temp.exists()
    finally:
        browser.close()
