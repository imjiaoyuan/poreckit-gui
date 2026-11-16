"""Whole-region matrix display modes computed in the render worker.

The Pearson and compartment displays summarize one square, diagonal contact
block (identical X and Y region). Both need enough bins for a stable
per-distance normalization and both must tolerate uncovered bins, so the
guards and NaN policy live here next to the math.
"""
import numpy as np

# Below this size the per-distance expected profile is too noisy for a
# meaningful correlation or compartment partition.
MIN_BINS = 50


def _square(matrix):
    a = np.asarray(matrix, dtype=np.float32)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("The selected display needs identical X and Y regions")
    if a.shape[0] < MIN_BINS:
        raise ValueError(f"The selected display needs at least {MIN_BINS} bins; zoom out or lower the resolution")
    return a


def pearson_matrix(matrix):
    """Pearson correlation of per-distance-normalized rows over the visible square.

    Normalization and correlation are delegated to the poreckit routine
    (O/E rows, pairwise finite bins), so the numbers match the library's
    compartment analysis.

    NaN policy: bins left entirely uncovered yield NaN correlations; those
    entries are replaced with 0 (the neutral point of a diverging colormap)
    and the diagonal is set to 1, since every profile correlates perfectly
    with itself. Results are clamped to [-1, 1].
    """
    a = _square(matrix)
    from poreckit import contact_correlation
    result = np.asarray(contact_correlation(a), dtype=np.float64)
    result[~np.isfinite(result)] = 0.0
    np.fill_diagonal(result, 1.0)
    return np.clip(result, -1.0, 1.0).astype(np.float32)


def compartment_matrix(matrix):
    """A/B saddle view from the sign of the first compartment eigenvector.

    The eigenvector is computed by the poreckit routine (O/E -> Pearson ->
    power iteration, sign fixed by the bin marginals). Uncovered bins become
    0 (unassigned). The displayed matrix is the outer product of the bin
    signs: +1 marks pairs in the same compartment (A-A or B-B) and -1 marks
    cross-compartment (A-B) pairs, which reproduces the plaid structure of a
    compartmentalized map. Hover values report the pair relationship; the
    per-bin assignment is the sign of the matching row or column.
    """
    a = _square(matrix)
    from poreckit import compartment_eigenvector
    vector = np.asarray(compartment_eigenvector(a), dtype=np.float64)
    vector[~np.isfinite(vector)] = 0.0
    sign = np.sign(vector).astype(np.float32)
    return np.outer(sign, sign).astype(np.float32)
