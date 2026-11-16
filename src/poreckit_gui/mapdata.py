"""Portable map payloads for exact hover values and genomic navigation."""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


def svg_plot_bounds(path):
    root = ET.parse(path).getroot()
    node = next((e for e in root.iter() if e.attrib.get("id") == "contact-map-bounds"), None)
    if node is None:
        return None
    shape = next((e for e in node.iter() if e.tag.endswith("path")), None)
    if shape is None:
        return None
    coords = np.asarray([float(n) for n in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", shape.attrib["d"])])
    points = coords.reshape(-1, 2)
    box = [float(np.min(points[:, 0])), float(np.min(points[:, 1])), float(np.ptp(points[:, 0])), float(np.ptp(points[:, 1]))]
    viewbox = [float(v) for v in root.attrib["viewBox"].split()]
    return {"bounds": box, "svg_size": viewbox[2:]}


def save_map_data(path, observed, displayed, control, x, y, resolution):
    path = str(Path(path).with_suffix(".data.npz"))
    values = {"observed": observed}
    if displayed is not observed:
        values["displayed"] = displayed
    if control is not None:
        values["control"] = control
    np.savez_compressed(path, **values)
    return {"path": path, "x": [x.chrom, x.start, x.end], "y": [y.chrom, y.start, y.end], "binsize": resolution}


class MapData:
    def __init__(self, spec):
        self.spec = spec
        with np.load(spec["path"], allow_pickle=False) as arrays:
            self.observed = arrays["observed"]
            self.displayed = arrays["displayed"] if "displayed" in arrays else self.observed
            self.control = arrays["control"] if "control" in arrays else None

    def at(self, fx, fy):
        x, y = self.spec["x"], self.spec["y"]
        bx = min(x[2]-1, max(x[1], int(x[1] + fx * (x[2]-x[1]))))
        by = min(y[2]-1, max(y[1], int(y[1] + fy * (y[2]-y[1]))))
        size = self.spec["binsize"]
        col, row = (bx-x[1]) // size, (by-y[1]) // size
        if row >= self.observed.shape[0] or col >= self.observed.shape[1]:
            return None
        result = {"x": (x[0], x[1]+col*size, min(x[2], x[1]+(col+1)*size)),
                  "y": (y[0], y[1]+row*size, min(y[2], y[1]+(row+1)*size)),
                  "observed": float(self.observed[row, col]), "value": float(self.displayed[row, col])}
        if self.control is not None:
            result["control"] = float(self.control[row, col])
        if self.spec.get("chromosomes"):
            def genomic(index):
                c = next(c for c in self.spec["chromosomes"] if c["offset"] <= index < c["offset"]+c["count"])
                lo = (index-c["offset"])*self.spec["genomic_binsize"]
                return c["name"], lo, min(c["length"], lo+self.spec["genomic_binsize"])
            result.update(x=genomic(col), y=genomic(row))
        return result
