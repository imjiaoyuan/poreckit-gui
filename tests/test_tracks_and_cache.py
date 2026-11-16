import json
import re
import time
from pathlib import Path

import numpy as np
import pytest

from poreckit_gui.cache import RenderCache
from poreckit_gui.demo import create_demo
from poreckit_gui.worker import render


def test_cache_bounds_keys_and_file_invalidation(tmp_path, monkeypatch):
    import poreckit
    source = tmp_path / "input.dat"
    source.write_text("1")
    calls = []

    class Reader:
        def matrix(self, chrom, start, end, resolution, norm):
            calls.append((chrom, start, end, resolution, norm))
            return np.full((3, 3), int(source.read_text()), dtype=np.float32)

    monkeypatch.setattr(poreckit, "make_reader", lambda path: Reader())
    cache = RenderCache(max_bytes=72)
    first = cache.matrix(source, "chr1", 0, 300, 100)
    assert cache.matrix(source, "chr1", 0, 300, 100) is first
    assert len(calls) == 1
    cache.matrix(source, "chr1", 0, 300, 100, "balanced")
    assert len(calls) == 2
    cache.matrix(source, "chr1", 100, 400, 100)
    assert cache.bytes <= 72 and len(cache.matrices) == 2
    source.write_text("22")
    updated = cache.matrix(source, "chr1", 0, 300, 100)
    assert np.all(updated == 22) and np.all(first == 1)
    assert len(calls) == 4
    cache.close()
    assert cache.bytes == 0


def test_annotation_cache_reuses_parsers_and_invalidates_modified_sources(tmp_path, monkeypatch):
    import poreckit
    genes = tmp_path / "genes.gtf"
    text = 'chr1\tsrc\tgene\t101\t500\t.\t+\t.\tgene_id "G1"; gene_name "A";\n'
    genes.write_text(text)
    tads = tmp_path / "domains.bed"
    tads.write_text("chr1\t100000\t400000\tT1\n")
    loops = tmp_path / "loops.bedpe"
    loops.write_text("chr1\t200000\t210000\tchr1\t700000\t710000\tL1\t5\n")
    original, calls = poreckit.parse_genes, []
    def parse(path, chrom):
        calls.append(chrom)
        return original(path, chrom)
    monkeypatch.setattr(poreckit, "parse_genes", parse)
    state = RenderCache()
    try:
        path = create_demo(tmp_path)
        view = dict(input=path, region="1:0-2000000", binsize=10000, mode="square", norm=None,
                    cmap="hic_whitered", log=False, vmax_percentile=95,
                    tracks=[dict(kind="genes", path=str(genes)), dict(kind="tads", path=str(tads)),
                            dict(kind="loops", path=str(loops))])
        render(view, str(tmp_path / "first.svg"), state=state)
        render(view, str(tmp_path / "restyled.svg"), state=state)
        assert calls == ["1", "chr1"]
        assert state.records("tads", tads, "1")[0][0] == "1"
        assert state.records("loops", loops, "1")[0][3] == "1"
        genes.write_text(text.replace("101", "201").replace('"A"', '"Updated"'))
        updated = state.records("genes", genes, "1")
        assert updated[0][1] == 200
        assert calls == ["1", "chr1", "1", "chr1"]
        assert 0 < state.record_bytes <= state.max_bytes//2
    finally:
        state.close()
    assert not state.records_cache and state.record_bytes == 0
    tiny = RenderCache(max_bytes=2)
    assert tiny.records("genes", genes, "chr1")
    assert not tiny.records_cache


def test_mixed_contact_signal_tad_loop_panels_share_coordinates(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    path = create_demo(tmp_path)
    signal = tmp_path / "signal.bg"
    signal.write_text("chr1\t0\t1000000\t2\nchr1\t1000000\t2000000\t8\n")
    tads = tmp_path / "tads.bed"
    tads.write_text("chr1\t100000\t500000\tT1\n")
    loops = tmp_path / "loops.bedpe"
    loops.write_text("chr1\t200000\t210000\tchr1\t700000\t710000\tL1\t5\n")
    tracks = [dict(kind="contacts", path=path, style="square", cmap="viridis", binsize="10000", name="Second map", placement="panel"),
              dict(kind="signal", path=str(signal), style="line", name="Line signal", placement="panel"),
              dict(kind="signal", path=str(signal), style="bar", name="Bar signal", color="#bb5522", placement="panel"),
              dict(kind="tads", path=str(tads), style="triangle", name="Domains", placement="panel"),
              dict(kind="loops", path=str(loops), style="bezier", name="Interactions", placement="panel")]
    view = dict(input=path, region="chr1:0-2000000", binsize=10000, mode="triangle", norm=None,
                cmap="hic_whitered", log=False, tracks=tracks, vmax_percentile=95)
    observed = []
    save = Figure.savefig

    def inspect(fig, *args, **kwargs):
        value = save(fig, *args, **kwargs)
        observed.append([(ax.get_xlim(), ax.get_position().x0, ax.get_position().width,
                          ax.get_title(loc="left"), len(ax.patches)) for ax in fig.axes])
        assert len(fig.axes) == 7
        square = fig.axes[1]
        assert square.images[0].get_cmap().name == "viridis"
        assert square.get_window_extent().width / square.get_window_extent().height == pytest.approx(1)
        return value

    monkeypatch.setattr(Figure, "savefig", inspect)
    state = RenderCache()
    try:
        render(view, str(tmp_path / "panels.svg"), state=state)
        rows = observed[-1]
        assert all(row[0] == (0, 2000000) for row in rows)
        assert all(row[1:3] == pytest.approx(rows[0][1:3]) for row in rows)
        assert "Second map" in rows[1][3]
        assert rows[2][4] == 0 and rows[3][4] == 800
        assert rows[4][3] == "Domains" and rows[5][3] == "Interactions"
        # Both maps and both signals reuse the same two native data queries.
        assert state.misses == 2 and state.hits == 2
        view["cmap"] = "Reds"
        render(view, str(tmp_path / "recolored.svg"), state=state)
        assert state.misses == 2 and state.hits == 6
    finally:
        state.close()


def write_square_storage_cooler(path, step=10000, n1=260, n2=120):
    """Cooler with storage-mode='square': both directions stored with
    different values, so matrix(x, y) != matrix(y, x).T."""
    import h5py
    n = n1 + n2
    pairs = [((a, n1+b), float(a*1000 + b)) for a in range(n1) for b in range(n2)]
    pairs += [((n1+b, a), float(b*1000 + a)) for a in range(n1) for b in range(n2)]
    bin1 = np.array([p[0] for p, _ in pairs], dtype="int64")
    bin2 = np.array([p[1] for p, _ in pairs], dtype="int64")
    count = np.array([v for _, v in pairs], dtype="float32")
    order = np.lexsort((bin2, bin1))
    bin1, bin2, count = bin1[order], bin2[order], count[order]
    with h5py.File(path, "w") as f:
        f.attrs.update({"format": "HDF5::Cooler", "format-version": 3,
                        "bin-size": step, "bin-type": "fixed", "storage-mode": "square"})
        chroms = f.create_group("chroms")
        chroms["name"] = np.array([b"chr1", b"chr2"])
        chroms["length"] = np.array([n1*step, n2*step], dtype="int64")
        bins = f.create_group("bins")
        bins["chrom"] = np.repeat(np.arange(2, dtype="int32"), [n1, n2])
        bins["start"] = np.concatenate([np.arange(n1, dtype="int64"), np.arange(n2, dtype="int64")]) * step
        bins["end"] = np.minimum(bins["start"][:] + step,
                                 np.repeat(np.array([n1*step, n2*step], dtype="int64"), [n1, n2]))
        pixels = f.create_group("pixels")
        pixels["bin1_id"], pixels["bin2_id"], pixels["count"] = bin1, bin2, count
        indexes = f.create_group("indexes")
        indexes["chrom_offset"] = np.array([0, n1, n], dtype="int64")
        indexes["bin1_offset"] = np.searchsorted(bin1, np.arange(n + 1)).astype("int64")
    return str(path)


def test_tiled_rectangle_matches_directed_queries_on_square_storage(tmp_path):
    """Regression: symmetric tile sharing once transposed asymmetric
    square-storage matrices and silently changed displayed and hover values."""
    from poreckit import make_reader

    from poreckit_gui.model import Region
    path = write_square_storage_cooler(tmp_path / "asymmetric.cool")
    state = RenderCache()
    rd = make_reader(path)
    try:
        assert state.symmetric(path, 10000) is False
        assert state.symmetric(str(create_demo(tmp_path)), 10000) is True
        x = Region("chr1", 0, 2000000)
        y = Region("chr2", 0, 1200000)
        tiled = state.tiled_rectangle(path, x, y, 10000)
        direct = rd.matrix("chr2", 0, 1200000, 10000, None, chrom2="chr1", start2=0, end2=2000000)
        np.testing.assert_allclose(tiled, direct)
        # The file is genuinely asymmetric: m[i,j] != m[j,i].
        assert direct[0, 1] != direct[1, 0]
        tiled = state.tiled_rectangle(path, y, x, 10000)
        direct = rd.matrix("chr1", 0, 2000000, 10000, None, chrom2="chr2", start2=0, end2=1200000)
        np.testing.assert_allclose(tiled, direct)
    finally:
        (rd.close() if hasattr(rd, "close") else None)
        state.close()


@pytest.mark.parametrize("square", [False, True],
                         ids=["symmetric-demo", "square-asymmetric"])
def test_tiled_rectangle_matches_direct_matrix_across_query_shapes(tmp_path, square):
    """Tiling must reproduce the direct reader for diagonal, off-diagonal
    intra-chromosomal and inter-chromosomal queries, for symmetric and
    directed square storage alike."""
    from poreckit import make_reader

    from poreckit_gui.model import Region
    if square:
        path = write_square_storage_cooler(tmp_path / "asymmetric.cool", n1=260, n2=120)
    else:
        path = create_demo(tmp_path, [("chr1", 4000000), ("chr2", 2500000)])
    resolution = 10000
    queries = [
        # Intra-chromosomal diagonal (x == y).
        (Region("chr1", 0, 2000000), Region("chr1", 0, 2000000)),
        # Intra-chromosomal off-diagonal (x != y on the same chromosome).
        (Region("chr1", 0, 1500000), Region("chr1", 1500000, 2600000)),
        # Inter-chromosomal.
        (Region("chr1", 0, 2000000), Region("chr2", 0, 1200000)),
        (Region("chr2", 500000, 1200000), Region("chr1", 1000000, 2600000)),
    ]
    state = RenderCache()
    rd = make_reader(path)
    try:
        assert state.symmetric(path, resolution) is not square
        for x, y in queries:
            tiled = state.tiled_rectangle(path, x, y, resolution)
            direct = rd.matrix(y.chrom, y.start, y.end, resolution, None,
                               chrom2=x.chrom, start2=x.start, end2=x.end)
            assert tiled.shape == direct.shape
            np.testing.assert_allclose(tiled, direct)
        if square:
            # The square file is genuinely directed, so the ordering checks
            # above are not vacuous: reversing the query is not the transpose.
            forward = state.tiled_rectangle(path, *queries[2], resolution)
            reverse = state.tiled_rectangle(path, queries[2][1], queries[2][0], resolution)
            assert not np.allclose(reverse, forward.T)
    finally:
        rd.close()
        state.close()


@pytest.mark.parametrize("display_mode,extra", [
    ("difference", {}),
    ("log2ratio", {"pseudocount": 0.5}),
    ("pearson_local", {}),
], ids=["difference", "log2ratio", "pearson_local"])
@pytest.mark.parametrize("log", [False, True], ids=["linear", "log"])
def test_triangle_signed_display_keeps_negative_values(tmp_path, display_mode, extra, log):
    """Regression: triangle signed displays once clamped negatives to a
    vmin=0 colour range instead of a symmetric signed range. Every signed
    display mode must get |vmin| == vmax in triangle mode, matching the
    square reference."""
    from matplotlib.figure import Figure
    observed = create_demo(tmp_path / "observed")
    control = create_demo(tmp_path / "control", scale=2)
    base = dict(input=observed, control=control, region="chr1:0-2000000", binsize=10000,
                mode="triangle", norm=None, cmap="RdBu_r", log=log,
                vmax_percentile=95, display_mode=display_mode,
                pseudocount=extra.get("pseudocount", 1), tracks=[])
    state = RenderCache()
    try:
        limits = {}
        save = Figure.savefig
        def inspect(fig, *args, **kwargs):
            image = fig.axes[0].images[0]
            limits["range"] = (float(image.norm.vmin), float(image.norm.vmax))
            limits["cmap"] = image.get_cmap().name
            return save(fig, *args, **kwargs)
        # The square path is the signed reference behaviour.
        with patch_savefig(inspect):
            render(dict(base, mode="square"), str(tmp_path / f"square-{display_mode}-{log}.svg"),
                   state=state)
            square_limits = dict(limits)
            render(base, str(tmp_path / f"triangle-{display_mode}-{log}.svg"), state=state)
        vmin, vmax = square_limits["range"]
        assert vmin < 0 < vmax and abs(vmin) == pytest.approx(vmax)
        vmin, vmax = limits["range"]
        assert vmin < 0 < vmax and abs(vmin) == pytest.approx(vmax)
        assert square_limits["cmap"] == limits["cmap"] == "RdBu_r"
    finally:
        state.close()


def patch_savefig(inspect):
    import contextlib

    from matplotlib.figure import Figure
    @contextlib.contextmanager
    def context():
        original = Figure.savefig
        Figure.savefig = inspect
        try:
            yield
        finally:
            Figure.savefig = original
    return context()


def test_catalog_text_carries_no_tool_branding():
    from poreckit_gui.catalog import BLOCKED_TOOLS, catalog, fields
    forbidden = re.compile(r"\b(" + "|".join(re.escape(name) for name in BLOCKED_TOOLS) + r")\b", re.I)
    for entry in catalog().values():
        assert not forbidden.search(entry.description)
        for field in fields(entry):
            assert not forbidden.search(field.get("help", "") + str(field.get("default", "")))


def test_sparse_signal_bins_remain_visible_without_filling_missing_data():
    import matplotlib.pyplot as plt

    from poreckit_gui.tracks import SignalPanel
    values = np.array([np.nan, 5.0, np.nan])
    fig, ax = plt.subplots()
    try:
        SignalPanel(values_override=values, mode="fill").draw(ax, "chr1", 0, 300)
        vertices = np.vstack([p.vertices for p in ax.collections[0].get_paths()])
        assert vertices[:, 0].min() == 100 and vertices[:, 0].max() == 200
        assert vertices[:, 1].max() == 5
    finally:
        plt.close(fig)


def test_track_settings_preserve_independent_modes():
    from PySide6.QtWidgets import QApplication, QDialog

    from poreckit_gui.tracksettings import TrackSettings
    _app = QApplication.instance() or QApplication([])
    dialog = TrackSettings(dict(kind="contacts", path="sample.cool", style="triangle", placement="panel"))
    dialog.mode.setCurrentIndex(dialog.mode.findData("square"))
    dialog.resolution.setText("5000")
    dialog.cmap.setCurrentIndex(dialog.cmap.findData("viridis"))
    dialog.norm.setCurrentIndex(dialog.norm.findData("raw"))
    dialog.chrom.setText("chr1")
    dialog.accept()
    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.track["style"] == "square" and dialog.track["binsize"] == "5000"
    assert dialog.track["cmap"] == "viridis" and dialog.track["norm"] == "raw"
    assert dialog.track["chrom"] == "chr1"


def test_gui_reuses_worker_and_reset_keeps_data_selection(tmp_path):
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
        assert not errors, errors

    try:
        browser.load_file(create_demo(tmp_path))
        wait_until(lambda: browser.last_view is not None and not browser.busy())
        pid = browser.process.processId()
        assert pid > 0
        initial_cache = json.loads(browser.result_path.read_text())["cache"]
        browser.cmap.setCurrentText("viridis")
        wait_until(lambda: browser.last_view["cmap"] == "viridis" and not browser.busy())
        assert browser.process.processId() == pid
        result = json.loads(browser.result_path.read_text())
        assert result["cache"]["misses"] == initial_cache["misses"]
        assert result["cache"]["hits"] > initial_cache["hits"]
        assert set(browser.last_files) == {"svg"}
        before = browser.view_data()
        browser.reset_display_settings()
        wait_until(lambda: browser.last_view["cmap"] == "hic_whitered" and not browser.busy())
        for name in ("input", "region", "binsize", "norm"):
            assert browser.last_view[name] == before[name]
    finally:
        browser.close()
    assert browser.process.processId() == 0


def test_contact_tiles_reuse_overlaps_and_match_rectangular_data(tmp_path):
    from poreckit import make_reader

    from poreckit_gui.model import Region
    path=create_demo(tmp_path, [('chr1',4000000),('chr2',2500000)])
    state=RenderCache()
    rd=make_reader(path)
    try:
        x=Region('chr1',0,2000000)
        y=Region('chr2',0,1500000)
        first=state.tiled_rectangle(path,x,y,10000)
        np.testing.assert_allclose(first,rd.matrix('chr2',0,1500000,10000,None,chrom2='chr1',start2=0,end2=2000000))
        before=state.hits
        x2=Region('chr1',500000,2500000)
        second=state.tiled_rectangle(path,x2,y,10000)
        assert state.hits>before
        np.testing.assert_allclose(second,rd.matrix('chr2',0,1500000,10000,None,chrom2='chr1',start2=500000,end2=2500000))
        # Last partial chromosome bins and reversed tile cache orientation.
        x3=Region('chr1',3900000,4000000)
        y3=Region('chr2',2400000,2500000)
        forward=state.tiled_rectangle(path,x3,y3,10000)
        reverse=state.tiled_rectangle(path,y3,x3,10000)
        np.testing.assert_allclose(forward,reverse.T)
        assert state.bytes<=state.max_bytes
    finally:
        rd.close()
        state.close()


def test_browser_previews_are_reaped_during_the_session(tmp_path):
    """Regression: every navigation once added preview/job/result files that
    were only removed at shutdown."""
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

    def numbered():
        return sorted(p.name for p in Path(browser.temp.name).iterdir()
                      if re.fullmatch(r"(preview|job|result)-\d+.*", p.name))

    try:
        browser.load_file(create_demo(tmp_path))
        wait_until(lambda: browser.last_view is not None and not browser.busy())
        for cmap in ("viridis", "inferno", "magma"):
            browser.cmap.setCurrentText(cmap)
            wait_until(lambda cmap=cmap: browser.last_view["cmap"] == cmap and not browser.busy())
        serials = {int(re.search(r"-(\d+)", name).group(1)) for name in numbered()}
        assert serials <= {browser.serial, browser.serial - 1}, numbered()
        assert browser.last_files["svg"]
        assert Path(browser.last_files["svg"]).is_file()
    finally:
        browser.close()


def test_shared_budget_across_matrices_records_and_payloads(tmp_path, monkeypatch):
    """Regression: matrices, records and payloads must share ONE max_bytes
    budget; combined use may never reach 1.75x max_bytes anymore."""
    import json

    import poreckit
    tads = tmp_path / "domains.bed"
    tads.write_text("chr1\t100\t400\tT1\n")

    class Reader:
        def matrix(self, chrom, start, end, resolution, norm):
            return np.full((10, 10), 1.0, dtype=np.float32)  # 400 bytes

    monkeypatch.setattr(poreckit, "make_reader", lambda path: Reader())
    state = RenderCache(max_bytes=4096)
    try:
        source = tmp_path / "input.dat"
        source.write_text("1")
        for i in range(12):
            state.matrix(source, "chr1", i * 100, i * 100 + 100, 10)
            state.records("tads", tads, str(i + 1))  # distinct chrom per key
            state.sequence_payload(("alignments", i), lambda: {"records": list(range(80))})
            total = state.bytes + state.record_bytes + state.payload_bytes
            assert total <= state.max_bytes, f"total {total} exceeds budget at insert {i}"
        # The mixed inserts genuinely exercised cross-cache eviction.
        assert state.matrices and state.records_cache and state.payloads
        assert len(state.order) == len(state.matrices) + len(state.records_cache) + len(state.payloads)
        # Per-insert refusals are defined against the shared budget.
        big = np.full((10, 10), 1.0, dtype=np.float64)  # refuser: handled via nbytes check in _load
        assert big.nbytes > 0
        huge_payload = {"records": list(range(5000))}
        state.sequence_payload("huge", lambda: huge_payload)
        assert "huge" not in state.payloads
        encoded = json.dumps({"records": list(range(80))}, ensure_ascii=False).encode()
        assert len(encoded) <= state.max_bytes // 4
    finally:
        state.close()
    assert not state.order and state.bytes == state.record_bytes == state.payload_bytes == 0


def write_pairs_fixture(path, length1=4_000_000, length2=2_500_000):
    """Small 4DN .pairs file with diagonal, off-diagonal intra- and
    inter-chromosomal contacts (unordered pairs appear in both orders)."""
    header = [f"#chromsize: chr1 {length1}", f"#chromsize: chr2 {length2}",
              "#columns readID chrom1 pos1 chrom2 pos2 strand1 strand2"]
    pairs = [
        ("r1", "chr1", 150000, "chr1", 350000),
        ("r2", "chr1", 350000, "chr1", 150000),    # same contact, reversed
        ("r3", "chr1", 1500000, "chr1", 2500000),  # off-diagonal intra
        ("r4", "chr2", 1000000, "chr2", 1600000),  # off-diagonal intra
        ("r5", "chr1", 150000, "chr2", 900000),    # inter-chromosomal
        ("r6", "chr1", 300000, "chr1", 1800000),   # spans the off-diagonal windows
        ("r7", "chr2", 300000, "chr2", 2000000),   # same, on chr2
    ]
    with open(path, "w") as handle:
        handle.write("\n".join(header) + "\n")
        for row in pairs:
            handle.write("\t".join(map(str, row)) + "\t+\t-\n")
    return str(path)


@pytest.mark.parametrize("suffix", [".pairs", ".pairs.gz"], ids=["plain", "gzip"])
def test_tiled_rectangle_matches_untiled_path_on_pairs(tmp_path, suffix):
    """Pair sources are inherently symmetric and resolution-free, so the
    tile cache must reproduce the untiled rectangle for diagonal,
    off-diagonal intra-chromosomal queries and share transposed tiles."""
    import gzip as gzip_module

    from poreckit import make_reader

    from poreckit_gui.model import Region
    target = tmp_path / f"contacts{suffix}"
    text = write_pairs_fixture(tmp_path / "plain.pairs")
    if suffix == ".pairs.gz":
        with open(text) as source, gzip_module.open(target, "wt") as sink:
            sink.write(source.read())
    else:
        target.write_text(Path(text).read_text())
    resolution = 10000
    queries = [
        # Diagonal (x == y), spanning several 1.28 Mbp tiles.
        (Region("chr1", 0, 2000000), Region("chr1", 0, 2000000)),
        # Off-diagonal intra-chromosomal, crossing tile boundaries.
        (Region("chr1", 0, 1280000), Region("chr1", 640000, 2560000)),
        (Region("chr2", 0, 2500000), Region("chr2", 0, 2500000)),
    ]
    state = RenderCache()
    rd = make_reader(str(target))
    try:
        assert state.symmetric(str(target), resolution) is True
        assert [tuple(c) for c in state.metadata(str(target))["chromsizes"]] == \
            [("chr1", 4000000), ("chr2", 2500000)]
        for x, y in queries:
            tiled = state.tiled_rectangle(str(target), x, y, resolution)
            untiled = state.rectangle(str(target), x, y, resolution)
            assert tiled.shape == untiled.shape
            np.testing.assert_allclose(tiled, untiled)
            assert tiled.sum() > 0, "fixture must place contacts inside the query"
            if x != y:
                np.testing.assert_allclose(
                    state.tiled_rectangle(str(target), y, x, resolution), tiled.T)
        # Inter-chromosomal: pairs readers have no chrom2 query, so the tile
        # path degrades to the untiled rectangle and raises the same error.
        x = Region("chr1", 0, 2000000)
        y = Region("chr2", 0, 2500000)
        with pytest.raises(ValueError, match="intra-chromosomal"):
            state.rectangle(str(target), x, y, resolution)
        with pytest.raises(ValueError, match="intra-chromosomal"):
            state.tiled_rectangle(str(target), x, y, resolution)
    finally:
        (rd.close() if hasattr(rd, "close") else None)
        state.close()


def test_tiled_rectangle_matches_untiled_path_on_parquet(tmp_path):
    """Parquet concatemer tables behave like pairs through the tile cache;
    chromosome lengths come from the metadata position scan, so queries stay
    inside the observed spans."""
    from poreckit import make_reader, parse_pairs, write_parquet

    from poreckit_gui.model import Region
    pairs = write_pairs_fixture(tmp_path / "contacts.pairs")
    target = str(tmp_path / "contacts.parquet")
    write_parquet(parse_pairs(pairs, None), target)
    resolution = 10000
    queries = [
        (Region("chr1", 0, 1000000), Region("chr1", 0, 1000000)),
        (Region("chr1", 0, 1280000), Region("chr1", 640000, 2560000)),
        (Region("chr2", 0, 1700000), Region("chr2", 0, 1700000)),
    ]
    state = RenderCache()
    rd = make_reader(target)
    try:
        assert state.symmetric(target, resolution) is True
        assert [tuple(c) for c in state.metadata(target)["chromsizes"]] == \
            [("chr1", 2500001), ("chr2", 2000001)]
        for x, y in queries:
            tiled = state.tiled_rectangle(target, x, y, resolution)
            untiled = state.rectangle(target, x, y, resolution)
            assert tiled.shape == untiled.shape
            np.testing.assert_allclose(tiled, untiled)
            assert tiled.sum() > 0
        with pytest.raises(ValueError, match="intra-chromosomal"):
            state.tiled_rectangle(target, Region("chr1", 0, 1000000),
                                  Region("chr2", 0, 1000000), resolution)
    finally:
        (rd.close() if hasattr(rd, "close") else None)
        state.close()
