"""Pearson and compartment (A/B) whole-region display modes."""
import numpy as np
import pytest

from poreckit_gui.cache import RenderCache
from poreckit_gui.contactview import transform_matrix
from poreckit_gui.demo import create_demo
from poreckit_gui.matrixmodes import MIN_BINS, compartment_matrix, pearson_matrix
from poreckit_gui.worker import render


def demo_matrix(tmp_path):
    from poreckit import make_reader
    path = create_demo(tmp_path)
    reader = make_reader(path)
    try:
        return np.asarray(reader.matrix("chr1", 0, 2000000, 10000, None), dtype=np.float32)
    finally:
        reader.close()


def plaid_matrix(blocks=2, per_block=50, seed=3):
    """Dense square matrix with strong within-block, weak cross-block contacts."""
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(blocks), per_block)
    index = np.arange(labels.size)
    decay = 100 / (1 + np.abs(index[:, None] - index[None, :]) / 5)
    same = (labels[:, None] == labels[None, :]).astype(float)
    return rng.poisson(decay + 40 * same).astype(np.float32)


def test_pearson_matrix_is_symmetric_bounded_and_nan_safe(tmp_path):
    matrix = demo_matrix(tmp_path)
    result = pearson_matrix(matrix)
    assert result.dtype == np.float32
    assert result.shape == matrix.shape
    assert np.allclose(result, result.T, atol=1e-6)
    assert np.all(np.diag(result) == 1.0)
    assert np.isfinite(result).all()
    assert result.min() >= -1.0 and result.max() <= 1.0


def test_pearson_replaces_uncovered_rows_with_neutral_values():
    matrix = plaid_matrix()
    matrix[10] = 0.0  # bin with no contacts: NaN correlation profile
    matrix[:, 10] = 0.0
    result = pearson_matrix(matrix)
    assert np.isfinite(result).all()
    neutral = result.copy()
    np.fill_diagonal(neutral, 0.0)
    assert np.all(neutral[10] == 0.0) and np.all(neutral[:, 10] == 0.0)
    assert result[10, 10] == 1.0


def test_matrix_modes_do_not_modify_input(tmp_path):
    matrix = demo_matrix(tmp_path)
    matrix.setflags(write=False)
    for mode in (pearson_matrix, compartment_matrix):
        result = mode(matrix)
        assert result.dtype == np.float32
    assert np.isfinite(compartment_matrix(plaid_matrix())).all()


@pytest.mark.parametrize("mode", ["pearson", "compartment"])
def test_region_modes_reject_non_square_and_tiny_views(tmp_path, mode):
    matrix = demo_matrix(tmp_path)
    with pytest.raises(ValueError, match="identical X and Y regions"):
        transform_matrix(matrix[:, :100], None, mode)
    with pytest.raises(ValueError, match=f"at least {MIN_BINS} bins"):
        transform_matrix(matrix[:MIN_BINS - 1, :MIN_BINS - 1], None, mode)


def test_compartment_sign_partitions_plaid_blocks():
    matrix = plaid_matrix()
    shown = compartment_matrix(matrix)
    assert shown.min() >= -1.0 and shown.max() <= 1.0
    labels = np.repeat(np.arange(2), 50)
    same = (labels[:, None] == labels[None, :]).astype(np.float32) * 2 - 1
    assert (np.allclose(shown, same) or np.allclose(shown, -same))


def test_transform_matrix_routes_region_modes(tmp_path):
    matrix = demo_matrix(tmp_path)
    pearson = transform_matrix(matrix, None, "pearson")
    assert np.allclose(pearson, pearson.T, atol=1e-6)
    compartment = transform_matrix(matrix, None, "compartment")
    assert set(np.unique(compartment)) <= {-1.0, 0.0, 1.0}


@pytest.mark.parametrize("display_mode", ["pearson", "compartment"])
def test_worker_renders_region_modes(tmp_path, display_mode):
    path = create_demo(tmp_path)
    view = dict(input=path, region="chr1:0-2000000", binsize=10000, mode="square", norm=None,
                cmap="RdBu_r", log=False, vmax=None, vmax_percentile=None,
                display_mode=display_mode, tracks=[])
    state = RenderCache()
    try:
        result = render(view, str(tmp_path / f"{display_mode}.svg"), state=state)
        assert result["display"]["display_mode"] == display_mode
        with np.load(str(tmp_path / f"{display_mode}.data.npz"), allow_pickle=False) as data:
            shown = np.asarray(data["displayed"], dtype=np.float32)
        assert shown.shape[0] == shown.shape[1]
        assert np.isfinite(shown).all()
        if display_mode == "pearson":
            assert shown.min() >= -1.0 and shown.max() <= 1.0 and np.all(np.diag(shown) == 1.0)
        else:
            assert set(np.unique(shown)) <= {-1.0, 0.0, 1.0}
    finally:
        state.close()


def test_worker_rejects_region_mode_on_independent_axes(tmp_path):
    path = create_demo(tmp_path)
    view = dict(input=path, region="chr1:0-2000000", region_y="chr1:0-1000000",
                binsize=10000, mode="square", norm=None, cmap="RdBu_r", log=False,
                display_mode="pearson", tracks=[])
    with pytest.raises(ValueError, match="identical X and Y regions"):
        render(view, str(tmp_path / "guard.svg"), state=None)
