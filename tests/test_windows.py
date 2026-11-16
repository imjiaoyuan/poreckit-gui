import time
from pathlib import Path

import pytest

from poreckit_gui.demo import create_demo
from poreckit_gui.model import Region
from poreckit_gui.viewport import Viewport
from poreckit_gui.windows import translated_viewport, WindowController


def test_linked_coordinates_require_matching_genomes_and_handle_chr_aliases():
    a = dict(assembly="hg38", chromsizes=[("chr1", 1000), ("chr2", 2000)])
    b = dict(assembly="GRCh38", chromsizes=[("1", 1000), ("2", 2000)])
    view = Viewport(Region("chr1", 100, 500), Region("chr2", 500, 1500))
    mapped = translated_viewport(view, a, b)
    assert mapped.x.chrom == "1" and mapped.y.chrom == "2"
    with pytest.raises(ValueError, match="assemblies differ"):
        translated_viewport(view, a, {**b, "assembly": "mm39"})
    with pytest.raises(ValueError, match="lengths differ"):
        translated_viewport(view, a, {**b, "chromsizes": [("1", 999), ("2", 2000)]})


def test_gui_linked_windows_keep_samples_and_queue_latest_navigation(tmp_path):
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    first, second = Browser(), Browser()
    windows = [first, second]
    errors = []
    for b in windows:
        b.show()
        b.auto_update.setChecked(False)
        b.error = lambda e: errors.append(str(e))
    def wait():
        until = time.monotonic()+30
        while any(b.busy() for b in windows) and time.monotonic() < until:
            app.processEvents()
            time.sleep(.01)
        app.processEvents()
        assert not any(b.busy() for b in windows)
        assert not errors, errors
    try:
        a = create_demo(tmp_path / "a", [("chr1", 1200000), ("chr2", 800000)])
        b = create_demo(tmp_path / "b", [("1", 1200000), ("2", 800000)], scale=2)
        first.load_file(a)
        second.load_file(b)
        wait()
        second.cmap.setCurrentText("viridis")
        WindowController.link_all(True)
        initial = [w.serial for w in windows]
        target = Viewport(Region("chr1", 200000, 600000), Region("chr2", 100000, 500000))
        first.navigate_viewport(target)
        wait()
        assert first.current_viewport() == target
        assert second.current_viewport() == Viewport(Region("1", 200000, 600000), Region("2", 100000, 500000))
        assert second.input_path == b
        assert second.cmap.currentText() == "viridis"
        assert [w.serial for w in windows] == [n+1 for n in initial]
        # A receiver already rendering keeps only the newest pending locus.
        old_first = first.serial
        second.render()
        assert second.busy()
        second.windows.receive(first, Viewport(Region("chr1", 0, 300000), Region("chr2", 0, 300000)))
        final = Viewport(Region("chr1", 300000, 700000), Region("chr2", 200000, 600000))
        second.windows.receive(first, final)
        wait()
        assert second.current_viewport().x == Region("1", 300000, 700000)
        assert second.current_viewport().y == Region("2", 200000, 600000)
        assert first.serial == old_first
        second.metadata["assembly"] = "different-genome"
        serial = second.serial
        first.navigate_viewport(Viewport(Region("chr1", 0, 400000), Region("chr1", 0, 400000)))
        wait()
        assert second.serial == serial
        assert "assemblies differ" in second.statusBar().currentMessage()
    finally:
        for b in windows:
            b.close()
    assert first.windows not in WindowController.instances
    assert second.windows not in WindowController.instances


def test_duplicated_temporary_data_survives_original_window_close():
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    first = Browser()
    first.show()
    first.auto_update.setChecked(False)
    windows = [first]
    errors = []
    first.error = lambda e: errors.append(str(e))
    def wait():
        until = time.monotonic()+30
        while any(w.busy() for w in windows) and time.monotonic() < until:
            app.processEvents()
            time.sleep(.01)
        app.processEvents()
        assert not any(w.busy() for w in windows)
        assert not errors, errors
    try:
        first.open_demo()
        wait()
        key = first.annotations.model.add_layer("Draft loops")
        first.annotations.model.add(key, dict(xchrom="chr1", x0=100000, x1=200000,
                                             ychrom="chr1", y0=400000, y1=500000, name="Loop"))
        duplicate = first.windows.create(duplicate=True)
        assert duplicate is not None
        duplicate.error = lambda e: errors.append(str(e))
        windows.append(duplicate)
        wait()
        assert duplicate.input_path != first.input_path
        assert Path(duplicate.input_path).resolve().is_relative_to(Path(duplicate.temp.name).resolve())
        assert duplicate.annotations.model.layers[0]["name"] == "Draft loops"
        assert duplicate.annotations.model.layers is not first.annotations.model.layers
        first.annotations.model.mark_saved()
        first.close()
        assert Path(duplicate.input_path).is_file()
        duplicate.render()
        wait()
        assert duplicate.image.map_data.observed.shape == (200, 200)
    finally:
        for b in windows:
            b.annotations.model.mark_saved()
            b.close()
