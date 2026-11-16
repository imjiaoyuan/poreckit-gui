"""Data-backed two-dimensional contact view and basic control comparisons."""
from dataclasses import dataclass

import numpy as np
from poreckit import rasterize
from poreckit.tracks.base import BaseTrack


DISPLAY_MODES = {
    "observed": "Observed", "control": "Control", "split": "Observed / control split",
    "difference": "Observed − control", "ratio": "(Observed + p) / (Control + p)",
    "log2ratio": "log2((Observed + p) / (Control + p))",
    "oe_local": "Observed / expected (visible region)", "pearson_local": "Pearson (visible region)",
    "expected": "Expected (file)", "oe": "Observed / expected (file)",
    "log2oe": "log2(Observed / expected) (file)",
    "oe_p": "(Observed + p) / (Expected + p) (file)",
    "residual": "Observed − expected (file)",
    "control_expected": "Control expected (file)", "control_oe": "Control / expected (file)",
    "split_oe": "Observed / control O/E split (file)",
    "difference_oe": "Observed O/E − control O/E (file)",
    "ratio_oe": "Observed O/E / control O/E (file)",
}
EXPECTED_MODES = {"expected", "oe", "log2oe", "oe_p", "residual", "split_oe", "difference_oe", "ratio_oe"}
CONTROL_EXPECTED_MODES = {"control_expected", "control_oe", "split_oe", "difference_oe", "ratio_oe"}
CONTROL_MODES = {"control", "split", "difference", "ratio", "log2ratio"} | CONTROL_EXPECTED_MODES
SIGNED_MODES = {"difference", "log2ratio", "pearson_local", "log2oe", "residual", "difference_oe"}


def transform_matrix(observed, control, mode, pseudocount=1, expected=None, control_expected=None):
    a = np.asarray(observed, dtype=np.float32)
    def divide(numerator, denominator):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.divide(numerator, denominator, out=np.full_like(a, np.nan), where=denominator > 0)
    if mode in EXPECTED_MODES | CONTROL_EXPECTED_MODES:
        if pseudocount < 0:
            raise ValueError("Pseudocount must be nonnegative")
        if mode in EXPECTED_MODES and expected is None:
            raise ValueError("This display requires expected values stored in the file")
        if mode in CONTROL_EXPECTED_MODES and (control_expected is None or control is None):
            raise ValueError("This display requires control expected values stored in the file")
        if mode == "expected":
            return expected
        if mode == "control_expected":
            return control_expected
        if mode == "residual":
            return a-expected
        if mode == "control_oe":
            return divide(control, control_expected)
        if mode == "oe_p":
            return divide(a+pseudocount, expected+pseudocount)
        oe = divide(a, expected)
        if mode == "oe":
            return oe
        if mode == "log2oe":
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.where(oe > 0, np.log2(oe), np.nan)
        control_oe = divide(control, control_expected)
        if mode == "difference_oe":
            return oe-control_oe
        if mode == "ratio_oe":
            return divide(oe, control_oe)
        if mode == "split_oe":
            if a.shape[0] != a.shape[1]:
                raise ValueError("Split view needs a square, diagonal region")
            return np.where(np.triu(np.ones(a.shape, dtype=bool)), oe, control_oe)
    if mode == "observed":
        return a
    if mode in {"oe_local", "pearson_local"}:
        from poreckit import obs_exp
        oe = np.asarray(obs_exp(a), dtype=np.float32)
        if mode == "oe_local":
            return oe
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.corrcoef(oe).astype(np.float32)
    if control is None:
        raise ValueError("Load a control map for the selected display")
    b = np.asarray(control, dtype=np.float32)
    if a.shape != b.shape:
        raise ValueError("Observed and control matrices must have matching bins")
    if mode == "control":
        return b
    if mode == "difference":
        return a - b
    if mode == "split":
        if a.shape[0] != a.shape[1]:
            raise ValueError("Split view needs a square, diagonal region")
        return np.where(np.triu(np.ones(a.shape, dtype=bool)), a, b)
    if mode in {"ratio", "log2ratio"}:
        if pseudocount < 0:
            raise ValueError("Pseudocount must be nonnegative")
        with np.errstate(divide="ignore", invalid="ignore"):
            result = np.divide(a + pseudocount, b + pseudocount,
                               out=np.full_like(a, np.nan), where=(b + pseudocount) != 0)
            if mode == "log2ratio":
                result = np.where(result > 0, np.log2(result), np.nan)
        return result
    raise ValueError(f"Unknown matrix display: {mode}")


@dataclass
class ContactPanel(BaseTrack):
    x: object = None
    y: object = None
    matrix: object = None
    cmap: str = "hic_whitered"
    log: bool = False
    vmax: float | None = None
    percentile: float | None = 95
    signed: bool = False
    display: str = "observed"
    grid: bool = False
    resolution: int = 0
    matrix_x_start: int | None = None
    matrix_y_start: int | None = None

    def region_tuple(self):
        return self.x.chrom, self.x.start, self.x.end

    def draw(self, ax, *args):
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
        data = np.asarray(self.matrix, dtype=np.float32)
        positive = data[np.isfinite(data) & (data > 0)]
        if self.x == self.y and data.shape[0] == data.shape[1] and not self.signed:
            positive = data[np.triu_indices_from(data, k=1)]
            positive = positive[np.isfinite(positive) & (positive > 0)]
        if self.signed:
            finite = np.abs(data[np.isfinite(data)])
            upper = self.vmax or (float(np.percentile(finite, self.percentile or 100)) if finite.size else 1)
            upper = max(upper, 1e-12)
            lower = -upper
        else:
            upper = self.vmax
            if upper is None and self.percentile is not None:
                upper = float(np.percentile(positive, self.percentile)) if positive.size else 1
            lower = 0
        image = np.asarray(rasterize(data, 800))
        if self.log:
            image = np.sign(image) * np.log1p(np.abs(image))
            lower = np.sign(lower) * np.log1p(abs(lower))
            upper = None if upper is None else np.log1p(upper)
        x0 = self.x.start if self.matrix_x_start is None else self.matrix_x_start
        y0 = self.y.start if self.matrix_y_start is None else self.matrix_y_start
        im = ax.imshow(np.ma.masked_invalid(image), origin="upper", aspect="equal", interpolation="nearest",
                       extent=[x0, x0+data.shape[1]*self.resolution if self.resolution else self.x.end,
                               y0+data.shape[0]*self.resolution if self.resolution else self.y.end, y0],
                       cmap=plt.get_cmap(self.cmap).with_extremes(bad="#e5e5e5"), vmin=lower, vmax=upper)
        ax.set_xlim(self.x.start, self.x.end)
        ax.set_ylim(self.y.end, self.y.start)
        fmt = FuncFormatter(lambda v, _: f"{v / 1e6:g}")
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)
        ax.set_xlabel(f"{self.x.chrom} (Mb)")
        ax.set_ylabel(f"{self.y.chrom} (Mb)")
        ax.set_title(self.name)
        if self.grid:
            ax.grid(alpha=.25, linewidth=.4)
        cax = ax.inset_axes([1.025, .45, .018, .45])
        bar = ax.figure.colorbar(im, cax=cax)
        unit = "contact" if self.display in {"observed", "control", "split"} else self.display.replace("_local", " (region)")
        bar.set_label(("signed log1p " if self.log and self.signed else "log1p " if self.log else "") + unit, fontsize=8)
        bar.ax.tick_params(labelsize=7)
        bar.outline.set_visible(False)
