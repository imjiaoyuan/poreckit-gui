import math
import time

import h5py
import numpy as np
import pytest

from poreckit_gui.cache import RenderCache
from poreckit_gui.demo import create_demo
from poreckit_gui.genome import draw_genome, genome_data
from poreckit_gui.mapdata import MapData


def test_sparse_coarsening_respects_chromosome_edges_and_storage(tmp_path, monkeypatch):
    import poreckit_gui.genome as genome
    monkeypatch.setattr(genome, "CHUNK", 3)
    chroms = [("chr1", 35000), ("chr2", 25000)]
    path = create_demo(tmp_path, chroms)
    target = [0, 0, 1, 1, 2, 2, 3]
    expected = np.zeros((4, 4))
    with h5py.File(path, "r") as f:
        for b1, b2, count in zip(f["pixels/bin1_id"][:], f["pixels/bin2_id"][:], f["pixels/count"][:], strict=False):
            a, b = target[b1], target[b2]
            expected[a, b] += count
            if a != b:
                expected[b, a] += count
    state = RenderCache()
    try:
        matrix, layout, resolution = genome_data(path, chroms, [10000], 20000, state)
        np.testing.assert_array_equal(matrix, expected)
        again, _, _ = genome_data(path, [list(c) for c in chroms], [10000], 20000, state)
        assert state.hits == 1
        assert again is matrix
        result = draw_genome(matrix, layout, resolution, str(tmp_path / "overview.svg"))
        data = MapData(result["map_data"])
        assert data.at(.875, .375)["x"] == ("chr2", 20000, 25000)
        assert data.at(.875, .375)["y"] == ("chr1", 20000, 35000)
        assert data.at(1, 1)["value"] == expected[-1, -1]
    finally:
        state.close()
    with h5py.File(path, "r+") as f:
        f.attrs["storage-mode"] = "square"
    square, _, _ = genome_data(path, chroms, [10000], 20000)
    assert square[2, 0] == 0 and square[0, 2] > 0
    with h5py.File(path, "r+") as f:
        f["pixels/bin1_id"][0] = 999
    with pytest.raises(ValueError, match="Invalid contact bin"):
        genome_data(path, chroms, [10000], 20000)


def test_square_storage_overview_matches_directed_pixels_across_chunks(tmp_path, monkeypatch):
    """The whole-genome tiled path must honour directed square storage: every
    stored pixel is added exactly once at (row, col) with no mirroring, even
    when the bin and pixel tables are read in several chunks."""
    from test_tracks_and_cache import write_square_storage_cooler

    import poreckit_gui.genome as genome
    # Small CHUNK so both the bin lookup and the pixel loop cross chunk
    # boundaries (65 bins, 2*40*25=2000 stored pixels).
    monkeypatch.setattr(genome, "CHUNK", 17)
    n1, n2, step = 40, 25, 10000
    path = write_square_storage_cooler(tmp_path / "asymmetric.cool", step=step, n1=n1, n2=n2)
    chroms = [("chr1", n1*step), ("chr2", n2*step)]
    # With overview resolution == bin step the genome layout offsets the
    # bins exactly like the file does, so target ids are the stored ones.
    expected = np.zeros((n1+n2, n1+n2), dtype=np.float64)
    with h5py.File(path, "r") as f:
        b1 = f["pixels/bin1_id"][:]
        b2 = f["pixels/bin2_id"][:]
        counts = f["pixels/count"][:]
    for a, b, value in zip(b1, b2, counts, strict=False):
        expected[a, b] += float(value)
    assert expected[3, n1+5] == 3005 and expected[n1+5, 3] == 5003
    state = RenderCache()
    try:
        matrix, layout, resolution = genome_data(path, chroms, [step], step, state)
        assert resolution == step
        assert [c["name"] for c in layout] == ["chr1", "chr2"]
        assert layout[1]["offset"] == n1
        np.testing.assert_array_equal(matrix, expected)
    finally:
        state.close()


def test_native_reader_overview_uses_rectangular_trans_blocks(tmp_path, monkeypatch):
    import poreckit
    closed, calls = [], []
    class Reader:
        def matrix(self, chrom, start, end, resolution, norm, chrom2=None, start2=None, end2=None):
            calls.append((chrom, chrom2))
            rows = math.ceil((end-start)/resolution)
            cols = math.ceil((end2-start2)/resolution) if chrom2 else rows
            return np.full((rows, cols), 7 if chrom2 else 2)
        def close(self):
            closed.append(True)
    monkeypatch.setattr(poreckit, "make_reader", lambda path: Reader())
    matrix, _, _ = genome_data(str(tmp_path / "test.hic"), [("1", 30000), ("2", 15000)], [20000])
    assert len(calls) == 3 and closed
    assert matrix.shape == (3, 3)
    np.testing.assert_array_equal(matrix[:2, 2], [7, 7])
    np.testing.assert_array_equal(matrix[2, :2], [7, 7])
    with pytest.raises(ValueError, match="coarser"):
        genome_data(str(tmp_path / "test.hic"), [("1", 1000000000)], [1], 1)


def test_overview_gui_navigation_and_exports(tmp_path, monkeypatch):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    b = Browser()
    b.show()
    b.auto_update.setChecked(False)
    errors = []
    b.error = lambda e: errors.append(str(e))
    def wait():
        until = time.monotonic()+30
        while b.busy() and time.monotonic() < until:
            app.processEvents()
            time.sleep(.01)
        app.processEvents()
        assert not b.busy(), b.log_view.toPlainText()
        assert not errors, errors
    try:
        b.load_file(create_demo(tmp_path, [("chr1", 1200000), ("chr2", 800000)]))
        wait()
        original = b.current_viewport()
        dialog = b.genome_dialog
        dialog.open()
        wait()
        assert b.current_viewport() == original
        assert dialog.image.map_data.observed.shape == (200, 200)
        for extension in ("svg", "png", "pdf"):
            out = tmp_path / ("whole-genome."+extension)
            monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, _out=out, **k: (str(_out), ""))
            dialog.export()
            wait()
            assert out.stat().st_size > 1000
        out = tmp_path / "genome.npz"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
        dialog.export_matrix()
        with np.load(out, allow_pickle=False) as data:
            assert data["chrom"][150] == "chr2"
            assert data["start"][150] == 300000
        image = dialog.image
        image.zoom_by(1.2)
        assert not image.fit_mode
        image.fit()
        point = image.mapFromScene(image.data_rect.x()+image.data_rect.width()*.25,
                                   image.data_rect.y()+image.data_rect.height()*.8)
        QTest.mouseDClick(image.viewport(), Qt.MouseButton.LeftButton, pos=point)
        wait()
        assert b.current_viewport().x.chrom == "chr1"
        assert b.current_viewport().y.chrom == "chr2"
        assert b.current_viewport().y.end == 800000
        assert b.image.map_data.observed.shape == (80, 120)
    finally:
        b.close()
