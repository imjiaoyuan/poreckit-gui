"""Portable session data and zero-based, half-open genomic navigation."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    def __post_init__(self):
        if not self.chrom or self.start < 0 or self.end <= self.start:
            raise ValueError("Region must satisfy 0 ≤ start < end, e.g. chr1:0-2,000,000")

    @classmethod
    def parse(cls, text):
        match = re.fullmatch(r"([^:\s]+):([\d,]+)-([\d,]+)", text.strip())
        if not match:
            raise ValueError("Region format: chr1:0-2,000,000 (zero-based, half-open)")
        chrom, start, end = match.groups()
        return cls(chrom, int(start.replace(",", "")), int(end.replace(",", "")))

    def text(self):
        return f"{self.chrom}:{self.start:,}-{self.end:,}"

    def window(self, start, span, length=None):
        span = max(1, round(span))
        if length is not None:
            span = min(span, length)
        start = max(0, round(start))
        if length is not None:
            start = min(start, length - span)
        return Region(self.chrom, start, start + span)

    def zoom(self, factor, length=None):
        span = (self.end - self.start) * factor
        return self.window((self.start + self.end - span) / 2, span, length)

    def pan(self, fraction, length=None):
        span = self.end - self.start
        return self.window(self.start + span * fraction, span, length)


def validate_view(view):
    region = Region.parse(view["region"])
    region_y = Region.parse(view.get("region_y") or view["region"])
    resolution = view.get("binsize")
    if resolution is not None and (not isinstance(resolution, int) or resolution <= 0):
        raise ValueError("Resolution must be a positive integer in bp or auto")
    if resolution and (region.end - region.start + resolution - 1) // resolution > 2500:
        raise ValueError("This region exceeds 2500 bins. Select a coarser resolution or a smaller region to limit memory use.")
    if view.get("mode") not in {"triangle", "square"}:
        raise ValueError("Unknown map mode")
    if view["mode"] == "triangle" and region_y != region:
        raise ValueError("Triangle view needs identical X and Y regions; use square view for independent axes")
    if resolution and (region_y.end - region_y.start + resolution - 1) // resolution > 2500:
        raise ValueError("The Y region exceeds 2500 bins; select a coarser resolution")
    if not Path(view["input"]).is_file():
        raise ValueError("Load a valid contact data file first")
    return region
