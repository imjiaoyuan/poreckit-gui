"""Genomic annotations arranged along the two axes of a square contact map."""
from __future__ import annotations

from matplotlib.transforms import Affine2D, Bbox


def placement(track, square):
    kind = track["kind"]
    selected = track.get("placement", "overlay" if kind in {"tads", "loops"} else "panel")
    if selected == "auto":
        selected = ("both" if kind in {"genes", "signal"} else "overlay") if square else "panel"
        if kind == "contacts":
            selected = "panel"
    if not square and selected in {"top", "left", "both"}:
        selected = "panel"
    if selected not in {"panel", "overlay", "top", "left", "both"}:
        raise ValueError(f"Unknown track placement: {selected}")
    if selected in {"top", "left", "both"} and kind not in {"genes", "signal"}:
        raise ValueError("Only genes and signals can be placed along matrix axes")
    return selected


def thickness(track):
    return max(.35, float(track.get("height", .3)) * 1.8)


class _VerticalAxes:
    """Let the original gene layout measure its future rotated dimensions."""
    def __init__(self, ax):
        self.ax = ax

    def __getattr__(self, name):
        return getattr(self.ax, name)

    def get_window_extent(self):
        box = self.ax.get_window_extent()
        return Bbox.from_bounds(box.x0, box.y0, box.height, box.width)


def draw_vertical_genes(track, ax, chrom, start, end):
    ax.set_xlim(start, end)
    track.draw(_VerticalAxes(ax), chrom, start, end)
    lane_limits = ax.get_ylim()
    swap = Affine2D().from_values(0, 1, 1, 0, 0, 0)
    for artist in [*ax.lines, *ax.patches]:
        artist.set_transform(swap + ax.transData)
    for text in ax.texts:
        text.set_transform(swap + ax.transData)
        text.set_rotation(90)
    ax.set_xlim(*lane_limits)
    ax.set_ylim(end, start)


class AxisTracks:
    def __init__(self, tracks, square):
        self.top = [t for t in tracks if placement(t, square) in {"top", "both"}]
        self.left = [t for t in tracks if placement(t, square) in {"left", "both"}]
        self.top_inches = sum(thickness(t) + .18 for t in self.top) + (.35 if self.top else 0)
        self.left_inches = sum(thickness(t) + .18 for t in self.left) + (.6 if self.left else 0)
        self.objects = {}

    def register(self, track, obj, y_obj=None):
        self.objects[id(track)] = (obj, y_obj or obj)

    def reserve(self, fig):
        if self.top or self.left:
            width, height = fig.get_size_inches()
            fig.subplots_adjust(left=(1 + self.left_inches) / width, right=1 - .85 / width,
                                top=1 - (.6 + self.top_inches) / height, bottom=.65 / height)

    def draw(self, fig, heatmap_ax, chrom, start, end, y_region=None):
        if not self.top and not self.left:
            return
        fig.canvas.draw()
        position = heatmap_ax.get_position()
        width, height = fig.get_size_inches()
        cursor = position.y1 + .35 / height
        for track in reversed(self.top):
            size = thickness(track) / height
            ax = fig.add_axes([position.x0, cursor, position.width, size], label=f"top-track-{id(track)}")
            obj = self.objects[id(track)][0]
            ax.set_xlim(start, end)
            obj.draw(ax, chrom, start, end)
            ax.set_xlim(start, end)
            self._decorate(ax, track, vertical=False)
            cursor += size + .18 / height
        cursor = position.x0 - .6 / width
        for track in reversed(self.left):
            size = thickness(track) / width
            ax = fig.add_axes([cursor - size, position.y0, size, position.height], label=f"left-track-{id(track)}")
            obj = self.objects[id(track)][1]
            y_chrom, y_start, y_end = (y_region.chrom, y_region.start, y_region.end) if y_region else (chrom, start, end)
            if track["kind"] == "genes":
                draw_vertical_genes(obj, ax, y_chrom, y_start, y_end)
            else:
                obj.draw_vertical(ax, y_chrom, y_start, y_end)
            ax.set_ylim(y_end, y_start)
            self._decorate(ax, track, vertical=True)
            cursor -= size + .18 / width

    @staticmethod
    def _decorate(ax, track, vertical):
        from pathlib import Path
        name = track.get("name") or Path(track["path"]).name
        for spine in ax.spines.values():
            spine.set_visible(False)
        if vertical:
            ax.tick_params(axis="y", left=False, labelleft=False)
            if track["kind"] == "genes":
                ax.set_xticks([])
            ax.tick_params(axis="x", labelsize=6)
            ax.set_ylabel(name, fontsize=7)
        else:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
            ax.tick_params(axis="y", labelsize=6)
            if track["kind"] == "genes":
                ax.set_yticks([])
            ax.set_title(name, loc="left", fontsize=7, pad=2)
