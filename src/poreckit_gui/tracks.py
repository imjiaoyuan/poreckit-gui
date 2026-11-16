"""Browser-specific panels built on poreckit's public track and reader APIs."""
from dataclasses import dataclass

import numpy as np
from matplotlib.ticker import MaxNLocator
from poreckit import TrackLoop
from poreckit.tracks.base import BaseTrack


def query_signal(path, chrom, start, end, bins, index=None):
    from .remote import source_path
    if not source_path(path).lower().endswith(('.bedgraph','.bg','.bedgraph.gz','.bg.gz')):
        from .bigwig import query_bigwig
        return query_bigwig(path,chrom,start,end,bins)
    from .model import Region
    from .sequenceio import tabix_lines
    values=np.zeros(bins,dtype=float)
    covered=np.zeros(bins,dtype=float)
    width=(end-start)/bins
    for line in tabix_lines(path,Region(chrom,start,end),index):
        if not line.strip() or line.startswith(('#','track ','browser ')):
            continue
        c=line.split()
        if c[0].removeprefix('chr')!=chrom.removeprefix('chr'):
            continue
        a,b,value=max(start,int(c[1])),min(end,int(c[2])),float(c[3])
        if b<=a or not np.isfinite(value):
            continue
        lo,hi=max(0,int((a-start)/width)),min(bins,int(np.ceil((b-start)/width)))
        pos=np.arange(lo,hi)
        overlap=np.minimum(b,start+(pos+1)*width)-np.maximum(a,start+pos*width)
        values[lo:hi]+=overlap*value
        covered[lo:hi]+=overlap
    return np.divide(values,covered,out=np.full(bins,np.nan),where=covered>0)


@dataclass
class SquareLoopOverlay(TrackLoop):
    y_chrom: str = ""

    def __post_init__(self):
        if self.loop_path:
            from poreckit import parse_loops
            self.loops = list(parse_loops(self.loop_path))
        self.y_chrom = self.y_chrom or self.chrom

    def draw(self, ax, chrom=None, start=None, end=None):
        ymin, ymax = sorted(ax.get_ylim())
        def same(a, b):
            return a.removeprefix("chr") == b.removeprefix("chr")
        seen = set()
        for c1, s1, e1, c2, s2, e2, _score in self.loops:
            for a, lo1, hi1, b, lo2, hi2 in ((c1, s1, e1, c2, s2, e2), (c2, s2, e2, c1, s1, e1)):
                key = (lo1, hi1, lo2, hi2)
                if key in seen or not same(a, self.chrom) or not same(b, self.y_chrom):
                    continue
                if hi1 <= start or lo1 >= end or hi2 <= ymin or lo2 >= ymax:
                    continue
                seen.add(key)
                self._mark(ax, lo1, hi1, lo2, hi2)

    def _mark(self, ax, s1, e1, s2, e2):
        x, y = (s1 + e1) / 2, (s2 + e2) / 2
        if self.style == "point":
            ax.scatter([x], [y], s=36, facecolors="none", edgecolors=self.colour, linewidths=self.lw, zorder=5)
        else:
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((s1, s2), max(1, e1-s1), max(1, e2-s2), fill=False,
                                   edgecolor=self.colour, linewidth=self.lw, zorder=5))


@dataclass
class SignalPanel(BaseTrack):
    path: str = ""
    chrom: str = ""
    name: str = "Signal"
    height: float = 0.25
    mode: str = "fill"
    color: str = "#397ca8"
    values_override: object = None

    def values(self, start, end):
        return np.asarray(self.values_override if self.values_override is not None else query_signal(self.path, self.chrom, start, end, 800), dtype=float)

    def draw(self, ax, chrom=None, start=None, end=None):
        values = self.values(start, end)
        if values.size:
            x = start + (np.arange(values.size) + 0.5) * (end - start) / values.size
            if self.mode == "fill":
                edges = np.linspace(start, end, values.size + 1)
                stepped_x = np.column_stack((edges[:-1], edges[1:])).ravel()
                stepped_y = np.repeat(values, 2)
                ax.fill_between(stepped_x, stepped_y, 0, color=self.color, alpha=0.7, linewidth=0)
                ax.plot(stepped_x, stepped_y, color=self.color, linewidth=0.6)
            elif self.mode == "line":
                ax.plot(x, values, color=self.color, linewidth=0.8)
                finite = np.isfinite(values)
                isolated = finite & ~np.r_[False, finite[:-1]] & ~np.r_[finite[1:], False]
                ax.scatter(x[isolated], values[isolated], s=6, color=self.color)
            elif self.mode == "bar":
                ax.bar(x, values, width=(end - start) / values.size, color=self.color, linewidth=0)
            elif self.mode == "points":
                ax.scatter(x, values, s=3, color=self.color)
            else:
                raise ValueError(f"Unknown signal mode: {self.mode}")
            ax.axhline(0, color="#888888", linewidth=0.4)
            if not np.isfinite(values).any():
                ax.text(.5, .5, "No signal in this region", transform=ax.transAxes, ha="center", fontsize=8)
        ax.yaxis.set_major_locator(MaxNLocator(2))
        ax.tick_params(axis="y", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def draw_vertical(self, ax, chrom=None, start=None, end=None):
        values = self.values(start, end)
        if values.size:
            y = start + (np.arange(values.size) + .5) * (end-start) / values.size
            if self.mode == "fill":
                edges = np.linspace(start, end, values.size + 1)
                ys = np.column_stack((edges[:-1], edges[1:])).ravel()
                vs = np.repeat(values, 2)
                ax.fill_betweenx(ys, vs, 0, color=self.color, alpha=.7, linewidth=0)
                ax.plot(vs, ys, color=self.color, linewidth=.6)
            elif self.mode == "line":
                ax.plot(values, y, color=self.color, linewidth=.8)
                finite = np.isfinite(values)
                isolated = finite & ~np.r_[False, finite[:-1]] & ~np.r_[finite[1:], False]
                ax.scatter(values[isolated], y[isolated], s=6, color=self.color)
            elif self.mode == "bar":
                ax.barh(y, values, height=(end-start)/values.size, color=self.color, linewidth=0)
            elif self.mode == "points":
                ax.scatter(values, y, s=3, color=self.color)
            else:
                raise ValueError(f"Unknown signal mode: {self.mode}")
            ax.axvline(0, color="#888888", linewidth=.4)
        ax.set_ylim(end, start)
        ax.invert_xaxis()
        ax.xaxis.set_major_locator(MaxNLocator(2))
