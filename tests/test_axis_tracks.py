
import pytest

from poreckit_gui.axistracks import placement
from poreckit_gui.cache import RenderCache
from poreckit_gui.demo import create_demo
from poreckit_gui.worker import render


def test_layout_follows_map_mode_and_preserves_explicit_placements():
    for kind in ("genes", "signal"):
        assert placement(dict(kind=kind, placement="auto"), True) == "both"
        assert placement(dict(kind=kind, placement="auto"), False) == "panel"
        assert placement(dict(kind=kind, placement="left"), False) == "panel"
        assert placement(dict(kind=kind, placement="panel"), True) == "panel"
    for kind in ("tads", "loops"):
        assert placement(dict(kind=kind, placement="auto"), True) == "overlay"
        assert placement(dict(kind=kind, placement="auto"), False) == "panel"


def test_square_axis_tracks_coordinate_projection_and_2d_boxes(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    path = create_demo(tmp_path)
    genes = tmp_path / "genes.gtf"
    genes.write_text('chr1\tdemo\texon\t200001\t400000\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GENE1";\n')
    signal = tmp_path / "signal.bg"
    signal.write_text("chr1\t0\t1000000\t2\nchr1\t1000000\t2000000\t8\n")
    tads = tmp_path / "tads.bed"
    tads.write_text("chr1\t100000\t500000\tT1\n")
    loops = tmp_path / "loops.bedpe"
    loops.write_text("chr1\t200000\t210000\tchr1\t700000\t710000\tL1\t5\n")
    tracks = [dict(kind=kind, path=str(source), name=name, placement="auto", color=color)
              for kind, source, name, color in [("genes", genes, "Genes", "#338855"),
                  ("signal", signal, "Signal", "#4477aa"), ("tads", tads, "Domains", "#d4a500"),
                  ("loops", loops, "Loops", "#299fca")]]
    view = dict(input=path, region="chr1:0-2000000", binsize=10000, mode="square", norm=None,
                cmap="hic_whitered", log=False, tracks=tracks, vmax_percentile=95)
    save = Figure.savefig
    captured = []

    def inspect(fig, *args, **kwargs):
        result = save(fig, *args, **kwargs)
        captured.append(fig)
        return result

    monkeypatch.setattr(Figure, "savefig", inspect)
    state = RenderCache()
    try:
        render(view, str(tmp_path / "square.svg"), state=state)
        fig = captured[-1]
        main = fig.axes[0]
        assert len(fig.axes) == 5
        top = [ax for ax in fig.axes if ax.get_label().startswith("top-track-")]
        left = [ax for ax in fig.axes if ax.get_label().startswith("left-track-")]
        assert len(top) == len(left) == 2
        main_pos = main.get_position()
        for ax in top:
            assert ax.get_xlim() == (0, 2000000)
            assert (ax.get_position().x0, ax.get_position().width) == pytest.approx((main_pos.x0, main_pos.width))
        for ax in left:
            assert ax.get_ylim() == (2000000, 0)
            assert (ax.get_position().y0, ax.get_position().height) == pytest.approx((main_pos.y0, main_pos.height))
        vertical_gene = next(ax for ax in left if ax.get_ylabel() == "Genes")
        exon = vertical_gene.patches[0]
        vertices = vertical_gene.transData.inverted().transform(exon.get_transform().transform(exon.get_path().vertices))
        assert vertices[:, 1].min() == pytest.approx(200000)
        assert vertices[:, 1].max() == pytest.approx(400000)
        vertical_signal = next(ax for ax in left if ax.get_ylabel() == "Signal")
        values = vertical_signal.lines[0].get_xdata()
        assert min(values) == 2 and max(values) == 8
        assert main.get_window_extent().width / main.get_window_extent().height == pytest.approx(1)
        assert main.get_xticklabels() and all(text.get_visible() for text in main.get_xticklabels())
        assert len(main.patches) == 2
        bounds = sorted(p.get_bbox().bounds for p in main.patches)
        assert bounds == [(200000, 700000, 10000, 10000), (700000, 200000, 10000, 10000)]
        domain = next(line for line in main.lines if len(line.get_xdata()) == 5)
        assert list(domain.get_xdata()) == [100000, 500000, 500000, 100000, 100000]
        assert list(domain.get_ydata()) == [100000, 100000, 500000, 500000, 100000]
        assert state.misses == 2  # one contact matrix and one shared X/Y signal query
        view["mode"] = "triangle"
        render(view, str(tmp_path / "triangle.svg"), state=state)
        fig = captured[-1]
        assert not any(ax.get_label().startswith(("top-track-", "left-track-")) for ax in fig.axes)
        assert len(fig.axes) == 6
        assert state.misses == 2
    finally:
        state.close()


def test_track_settings_exposes_axis_placements():
    from PySide6.QtWidgets import QApplication, QDialog

    from poreckit_gui.tracksettings import TrackSettings
    _app = QApplication.instance() or QApplication([])
    dialog = TrackSettings(dict(kind="genes", path="genes.gtf", placement="auto"))
    assert dialog.placement.currentData() == "auto"
    dialog.placement.setCurrentIndex(dialog.placement.findData("both"))
    dialog.accept()
    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.track["placement"] == "both"
