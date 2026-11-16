"""Independent genomic axes and navigation operations."""
from dataclasses import dataclass

from .model import Region


@dataclass(frozen=True)
class Viewport:
    x: Region
    y: Region

    @classmethod
    def from_view(cls, view):
        x = Region.parse(view["region"])
        return cls(x, Region.parse(view.get("region_y") or x.text()))

    @property
    def diagonal(self):
        return self.x == self.y

    def zoom(self, factor, fx=.5, fy=.5, lengths=None):
        lengths = lengths or {}
        def axis(region, fraction):
            span = region.end - region.start
            centre = region.start + span * fraction
            return region.window(centre - span * factor * fraction, span * factor, lengths.get(region.chrom))
        return Viewport(axis(self.x, fx), axis(self.y, fy))

    def pan(self, dx, dy, lengths=None):
        lengths = lengths or {}
        return Viewport(self.x.pan(dx, lengths.get(self.x.chrom)), self.y.pan(dy, lengths.get(self.y.chrom)))

    def select(self, x0, x1, y0, y1, lengths=None):
        lengths = lengths or {}
        def axis(region, a, b):
            lo, hi = sorted((max(0, min(1, a)), max(0, min(1, b))))
            span = region.end - region.start
            return region.window(region.start + lo * span, max(1, (hi-lo)*span), lengths.get(region.chrom))
        return Viewport(axis(self.x, x0, x1), axis(self.y, y0, y1))
