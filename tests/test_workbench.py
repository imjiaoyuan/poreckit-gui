import ast
import inspect
import json
import time
from pathlib import Path

import numpy as np
import pytest

from poreckit_gui.catalog import catalog, fields, cli_arguments, cli_command, resolve
from poreckit_gui.demo import create_demo
from poreckit_gui.execution import execute_plot
from poreckit_gui.plotdata import load_value


def entry(name):
    return next(e for e in catalog().values() if e.title == name)


def assert_assets(result, count=None):
    if count is not None:
        assert len(result["figures"]) == count
    assert result["figures"]
    for figure in result["figures"]:
        assert Path(figure["files"]["png"]).read_bytes().startswith(b"\x89PNG")
        assert Path(figure["files"]["pdf"]).read_bytes().startswith(b"%PDF")
        assert "<svg" in Path(figure["files"]["svg"]).read_text()


def test_catalog_covers_every_public_plot_definition_and_track():
    import poreckit
    import poreckit.tracks as tracks
    from poreckit.tracks.base import BaseTrack
    root = Path(poreckit.__file__).parent
    defined = set()
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and (node.name.startswith("plot_") or node.name == "persp_aggregate"):
                module = "poreckit." + ".".join(path.relative_to(root).with_suffix("").parts)
                defined.add(module + ":" + node.name)
    actual = {e.target for e in catalog().values() if e.kind == "api"}
    assert defined == actual
    expected = {obj.__module__ + ":" + obj.__name__ for name, obj in vars(tracks).items()
                if name.startswith("Track") and inspect.isclass(obj) and issubclass(obj, BaseTrack)}
    assert expected <= {e.target for e in catalog().values() if e.kind == "track"}
    for item in catalog().values():
        described = {field["name"] for field in fields(item)}
        if item.kind != "cli":
            expected = {"__extra__" if p.kind == p.VAR_KEYWORD else p.name
                        for p in inspect.signature(resolve(item.target)).parameters.values()
                        if p.name != "ax" and not p.name.startswith("_") and p.kind != p.VAR_POSITIONAL}
            assert described == expected


def test_cli_arrays_tuples_positional_and_paired_flags(tmp_path):
    path = create_demo(tmp_path)
    command = cli_command("heatmap")
    values = dict(contacts=path, chrom="chr1", start=0, end=2000000, log_=False,
                  signal=["a path.bw", "b.bw"], output=str(tmp_path / "out.pdf"))
    argv = cli_arguments("heatmap", values)
    with command.make_context("heatmap", argv) as ctx:
        assert ctx.params["log_"] is False
        assert ctx.params["signal"] == ("a path.bw", "b.bw")
    argv = cli_arguments("compare", dict(matrix1=path, matrix2=path, chrom="chr1", start=0, end=2000000,
                                        binsize=10000, output=str(tmp_path / "compare.pdf")))
    with cli_command("compare").make_context("compare", argv) as ctx:
        assert ctx.params["matrix1"] == path
    argv = cli_arguments("multiway", dict(source=path, output="out.json", viewpoint=["chr1", 0, 100],
                                          annotation=[["genes", path], ["peaks", path]]))
    with cli_command("multiway").make_context("multiway", argv) as ctx:
        assert ctx.params["viewpoint"] == ("chr1", 0, 100)
        assert ctx.params["annotation"] == (("genes", path), ("peaks", path))


def test_structured_data_loads_arrays_without_executing_code(tmp_path):
    path = tmp_path / "data.npz"
    np.savez(path, a=np.eye(4), b=np.ones(3))
    assert load_value({"$file": str(path), "key": "a"}).shape == (4, 4)
    nested = load_value({"matrix": {"$array": [[1, 2], [3, 4]]}, "names": ["a", "b"]})
    assert isinstance(nested["matrix"], np.ndarray)
    assert nested["names"] == ["a", "b"]
    circular = tmp_path / "circular.json"
    circular.write_text(json.dumps({"$file": str(circular)}))
    with pytest.raises(ValueError, match="Circular"):
        load_value({"$file": str(circular)})
    unsafe = tmp_path / "object.npy"
    np.save(unsafe, np.array([{"code": "no"}], dtype=object))
    with pytest.raises(ValueError):
        load_value({"$file": str(unsafe)})


@pytest.fixture
def contacts(tmp_path):
    from poreckit import ContactDataset, ContactRead, ContactFragment, write_contacts
    reads = []
    for i in range(20):
        fragments = tuple(ContactFragment("chr1", i * 1000 + offset, i * 1000 + offset + 500)
                          for offset in (0, 100000, 300000))
        reads.append(ContactRead(str(i), fragments))
    path = tmp_path / "contacts.tsv"
    write_contacts(ContactDataset(reads), str(path))
    return str(path)


@pytest.mark.parametrize("name", ["heatmap", "tads", "virtual4c", "balance"])
def test_cli_plot_execution_and_all_export_formats(tmp_path, name):
    path = create_demo(tmp_path)
    values = dict(contacts=path, chrom="chr1", start=0, end=1000000, binsize=10000)
    if name == "virtual4c":
        values["viewpoint"] = 200000
    result = execute_plot(dict(kind="cli", entry="cli:" + name, values=values, directory=str(tmp_path / name)))
    assert_assets(result, 2 if name == "tads" else 1)


@pytest.mark.parametrize("name", ["qc", "report", "multiway"])
def test_multi_figure_results_and_reports(contacts, tmp_path, name):
    values = {"source" if name == "multiway" else "pairs": contacts}
    if name == "qc":
        values.update(binsizes="10000,100000", min_contacts=1)
    if name == "report":
        values.update(binsizes="10000,100000", min_contacts=1)
    result = execute_plot(dict(kind="cli", entry="cli:" + name, values=values, directory=str(tmp_path / name)))
    assert_assets(result)
    if name == "qc":
        assert len(result["figures"]) == 3
    if name == "report":
        assert any(p.endswith(".html") for p in result["sidecars"])


def test_api_axes_multi_panel_and_track_stack(tmp_path):
    context = dict(input=create_demo(tmp_path), region="chr1:0-1,000,000", binsize=10000)
    plot = entry("plot_differential_matrix")
    result = execute_plot(dict(kind="api", entry=plot.key,
                              values=dict(m1={"$matrix": True}, m2={"$matrix": True}),
                              context=context, directory=str(tmp_path / "api")))
    assert_assets(result, 1)
    assert_assets(execute_plot(dict(kind="api", entry=entry("plot_aggregate_methods").key,
                                   values={}, directory=str(tmp_path / "multipanel"))), 1)
    stack = [{"entry": entry("TrackHeatmap").key, "values": dict(reader={"$reader": context["input"]},
              chrom="chr1", start=0, end=1000000, bin_size=10000)},
             {"entry": entry("TrackXAxis").key, "values": {}}]
    assert_assets(execute_plot(dict(kind="track", tracks=stack, context=context,
                                   directory=str(tmp_path / "stack"))), 1)


def test_workbench_render_zoom_pan_export_preset_and_cancel(tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication, QFileDialog
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest
    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.show()
    workbench = browser.plot_workbench
    browser.tabs.setCurrentWidget(workbench)
    errors = []
    workbench.error = lambda message: errors.append(str(message))

    def wait():
        end = time.monotonic() + 45
        while workbench.busy() and time.monotonic() < end:
            app.processEvents()
            time.sleep(.01)
        assert not workbench.busy(), workbench.log.toPlainText()

    try:
        key = entry("plot_aggregate_methods").key
        workbench.select_entry(key)
        workbench.run()
        wait()
        assert not errors, workbench.log.toPlainText()
        assert_assets(workbench.result, 1)
        initial = workbench.image.transform().m11()
        workbench.image.zoom_by(2)
        assert workbench.image.transform().m11() == pytest.approx(initial * 2)
        app.processEvents()
        scrollbar = workbench.image.horizontalScrollBar()
        before = scrollbar.value()
        viewport = workbench.image.viewport()
        centre = viewport.rect().center()
        QTest.mousePress(viewport, Qt.MouseButton.LeftButton, pos=centre)
        QTest.mouseMove(viewport, centre + QPoint(-70, 0), delay=20)
        QTest.mouseRelease(viewport, Qt.MouseButton.LeftButton, pos=centre + QPoint(-70, 0))
        assert scrollbar.value() != before
        svg = tmp_path / "export.svg"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(svg), "SVG (*.svg)"))
        workbench.export()
        assert "<svg" in svg.read_text()
        png = tmp_path / "viewport.png"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(png), "PNG (*.png)"))
        workbench.export_viewport()
        assert png.read_bytes().startswith(b"\x89PNG")
        archive_path = tmp_path / "results.zip"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(archive_path), "ZIP (*.zip)"))
        workbench.export_all()
        import zipfile
        with zipfile.ZipFile(archive_path) as archive:
            assert set(archive.namelist()) == {"plot-001.png", "plot-001.svg", "plot-001.pdf"}
        preset = tmp_path / "preset.json"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **kw: (str(preset), "JSON (*.json)"))
        workbench.save_preset()
        workbench.select_entry("cli:heatmap")
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **kw: (str(preset), "JSON (*.json)"))
        workbench.load_preset()
        assert workbench.selected.key == key
        previous = workbench.result
        workbench.select_entry(entry("plot_differential_matrix").key)
        workbench.form.set_values({"m1": {"$array": [1, 2, 3]}, "m2": {"$array": [[1, 2], [3, 4]]}})
        workbench.run()
        wait()
        assert errors
        assert workbench.result is previous
        assert "failed" in workbench.status.text()
        workbench.select_entry(key)
        workbench.run()
        workbench.cancel()
        wait()
        assert "cancelled" in workbench.status.text()
        assert workbench.result["figures"]
    finally:
        browser.close()


def test_matrix_reference_resolves_chromosome_alias_like_the_browser(tmp_path):
    from poreckit import make_reader
    path = create_demo(tmp_path, [("chr1", 40000)])
    reader = make_reader(str(path))
    try:
        expected = reader.matrix("chr1", 0, 40000, 10000, None)
    finally:
        reader.close()
    # The same region must work with either alias, as it does in the browser.
    for region in ("chr1:0-40000", "1:0-40000"):
        result = load_value({"$matrix": True}, context=dict(input=str(path), region=region, binsize=10000))
        np.testing.assert_allclose(result, expected)
        assert result.sum() > 0


def test_workbench_run_directories_are_reaped_during_the_session(tmp_path):
    """Regression: run directories and cancelled/failed runs once accumulated
    until the window closed."""
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    browser = Browser()
    browser.show()
    workbench = browser.plot_workbench
    browser.tabs.setCurrentWidget(workbench)
    errors = []
    workbench.error = lambda message: errors.append(str(message))

    def wait():
        end = time.monotonic() + 45
        while workbench.busy() and time.monotonic() < end:
            app.processEvents()
            time.sleep(.01)
        assert not workbench.busy(), workbench.log.toPlainText()

    def run_dirs():
        return sorted(p.name for p in Path(workbench.temp.name).glob("run-*"))

    try:
        key = entry("plot_aggregate_methods").key
        workbench.select_entry(key)
        for _ in range(3):
            workbench.run()
            wait()
            assert not errors, workbench.log.toPlainText()
        assert len(run_dirs()) == 2, run_dirs()
        # A cancelled run leaves no usable figures; its directory is dropped.
        workbench.run()
        workbench.cancel()
        wait()
        assert "cancelled" in workbench.status.text()
        assert len(run_dirs()) == 2, run_dirs()
        assert workbench.result["figures"]
    finally:
        browser.close()
