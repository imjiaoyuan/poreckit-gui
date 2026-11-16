import struct

import numpy as np
import pytest

from poreckit_gui.contactview import transform_matrix
from poreckit_gui.expected import expected_index, expected_matrix
from poreckit_gui.model import Region


def fake_hic(path, version):
    def text(value):
        return value.encode()+b"\0"
    def num(fmt, value):
        return struct.pack("<"+fmt, value)
    modern = version > 8
    header = b"HIC\0"+num("i", version)+num("q", 0)+text("test")
    header += b"\0"*16 if modern else b""
    header += num("i", 0)+num("i", 2)
    header += text("ALL")+num("q" if modern else "i", 40)
    header += text("chr1")+num("q" if modern else "i", 40000)
    header = header[:8]+num("q", len(header))+header[16:]
    def curve(norm, values, factor):
        data = text(norm) if norm else b""
        data += text("BP")+num("i", 10000)+num("q" if modern else "i", len(values))
        data += b"".join(num("f" if modern else "d", v) for v in values)
        data += num("i", 1)+num("i", 1)+num("f" if modern else "d", factor)
        return data
    footer = num("i", 0)+num("i", 1)+curve(None, [20, 10, 4, 2], 2)
    normalized = num("i", 1)+curve("KR", [12, 6, 3, 1.5], 3)
    path.write_bytes(header+num("q" if modern else "i", len(footer))+footer+normalized)


@pytest.mark.parametrize("version", [6, 7, 8, 9])
def test_hic_expected_curves_use_chromosome_factor_and_genomic_distance(tmp_path, version):
    path = tmp_path / "test.hic"
    fake_hic(path, version)
    x = Region("chr1", 10000, 40000)
    y = Region("1", 0, 20000)
    expected = expected_matrix(path, x, y, 10000)
    np.testing.assert_allclose(expected, [[5, 2, 1], [10, 5, 2]])
    np.testing.assert_allclose(expected_matrix(path, x, y, 10000, "KR"), [[2, 1, .5], [4, 2, 1]])
    # Values from the same bins are invariant when the viewport changes.
    sub = expected_matrix(path, Region("1", 20000, 40000), y, 10000)
    np.testing.assert_array_equal(sub, expected[:, 1:])
    with pytest.raises(ValueError, match="No file expected"):
        expected_matrix(path, x, y, 5000)
    with pytest.raises(ValueError, match="same chromosome"):
        expected_matrix(path, x, Region("2", 0, 20000), 10000)
    path.write_bytes(b"HIC\0")
    with pytest.raises(ValueError, match="Truncated"):
        expected_index(path)


def spec_hic(path, version):
    """Build a minimal but structurally faithful .hic file following the public .hic container spec."""
    def text(value):
        return value.encode()+b"\0"
    def num(fmt, value):
        return struct.pack("<"+fmt, value)
    modern = version > 8
    # Header: magic, version, footer position, genome id, [NVI for v9], attributes, chromosomes, resolutions.
    body = b""  # No matrix records; only the footer matters for expected values.
    header = b"HIC\0"+num("i", version)+num("q", 0)
    header += text("hg19")
    header += b"\0"*16 if modern else b""  # Normalization-vector index position and size (v9 only).
    header += num("i", 1)+text("software")+text("poreckit-gui-test")
    header += num("i", 2)
    header += text("chr1")+num("q" if modern else "i", 40000)
    header += text("chr2")+num("q" if modern else "i", 40000)
    header += num("i", 1)+num("i", 10000)+num("i", 0)
    def curve(norm, values, factor):
        data = text(norm) if norm else b""
        data += text("BP")+num("i", 10000)+num("q" if modern else "i", len(values))
        data += b"".join(num("f" if modern else "d", v) for v in values)
        data += num("i", 1)+num("i", 0)+num("f" if modern else "d", factor)
        return data
    frag_curve = text("FRAG")+num("i", 1)+num("q" if modern else "i", 2)
    frag_curve += b"".join(num("f" if modern else "d", v) for v in [7, 3])
    frag_curve += num("i", 1)+num("i", 1)+num("f" if modern else "d", 5)
    unnormalized = num("i", 2)+frag_curve+curve(None, [20, 10, 4, 2], 2)
    normalized = num("i", 1)+curve("KR", [12, 6, 3, 1.5], 3)
    master = num("i", 1)+text("0_0")+num("q", len(header))+num("i", 16)
    footer = num("q" if modern else "i", len(master)+len(unnormalized))
    footer += master+unnormalized
    trailer = num("i", 0)  # Empty normalization-vector index after the normalized expected vectors.
    header = header[:8]+num("q", len(header)+len(body))+header[16:]
    path.write_bytes(header+body+footer+normalized+trailer)


@pytest.mark.parametrize("version", [6, 7, 8, 9])
def test_hic_expected_footer_matches_format_spec(tmp_path, version):
    path = tmp_path / "spec.hic"
    spec_hic(path, version)
    chromosomes, curves = expected_index(path)
    assert chromosomes == ["chr1", "chr2"]
    # FRAG-unit vectors are ignored; both the raw and normalized BP curves are indexed.
    assert set(curves) == {("NONE", 10000), ("KR", 10000)}
    offset, count, dtype, factors = curves[("NONE", 10000)]
    assert count == 4
    assert dtype == ("<f4" if version > 8 else "<f8")
    assert factors == {0: 2.0}
    x = Region("chr1", 0, 20000)
    np.testing.assert_allclose(expected_matrix(path, x, x, 10000), [[10, 5], [5, 10]])
    np.testing.assert_allclose(expected_matrix(path, x, x, 10000, "KR"), [[4, 2], [2, 4]])


def test_hic_expected_rejects_unsupported_versions(tmp_path):
    path = tmp_path / "v5.hic"
    fake_hic(path, 5)
    with pytest.raises(ValueError, match="unsupported for version 5"):
        expected_index(path)
    path = tmp_path / "v10.hic"
    fake_hic(path, 10)
    with pytest.raises(ValueError, match="unsupported for version 10"):
        expected_index(path)


def test_expected_comparison_formulas():
    a = np.array([[4, 8], [8, 16]], dtype=np.float32)
    b = a*2
    e = np.full_like(a, 2)
    c = np.full_like(a, 8)
    np.testing.assert_allclose(transform_matrix(a, b, "oe", expected=e), a/2)
    np.testing.assert_allclose(transform_matrix(a, b, "control_oe", control_expected=c), b/8)
    np.testing.assert_allclose(transform_matrix(a, b, "difference_oe", expected=e, control_expected=c), a/4)
    np.testing.assert_allclose(transform_matrix(a, b, "ratio_oe", expected=e, control_expected=c), 2)
    np.testing.assert_allclose(transform_matrix(a, b, "log2oe", expected=e), np.log2(a/2))
    np.testing.assert_allclose(transform_matrix(a, b, "oe_p", 1, e), (a+1)/3)
    with pytest.raises(ValueError, match="requires expected"):
        transform_matrix(a, None, "oe")


def test_secondary_square_panel_uses_independent_y(tmp_path, monkeypatch):
    from poreckit_gui.contactview import ContactPanel
    from poreckit_gui.demo import create_demo
    from poreckit_gui.worker import render
    path = create_demo(tmp_path, [("chr1", 1200000), ("chr2", 800000)])
    panels = []
    original = ContactPanel.draw
    def draw(self, *args):
        panels.append(self)
        return original(self, *args)
    monkeypatch.setattr(ContactPanel, "draw", draw)
    render(dict(input=path, region="chr1:100000-900000", region_y="chr2:0-600000", mode="square",
                binsize=10000, norm=None, log=False, cmap="hic_whitered", vmax_percentile=95,
                tracks=[dict(kind="contacts", path=path, style="square", binsize=10000)]), str(tmp_path / "panel.svg"))
    assert len(panels) == 2
    assert panels[1].y.chrom == "chr2"
    assert panels[1].matrix.shape == (60, 80)
    np.testing.assert_array_equal(panels[0].matrix, panels[1].matrix)
