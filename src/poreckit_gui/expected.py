"""Read distance-dependent expected curves from local HIC file footers."""
import functools
import struct
from pathlib import Path

import numpy as np


class BinaryReader:
    def __init__(self, stream):
        self.stream = stream
        self.size = stream.seek(0, 2)
        stream.seek(0)

    def read(self, n):
        value = self.stream.read(n)
        if len(value) != n:
            raise ValueError("Truncated HIC expected-value data")
        return value

    def number(self, fmt):
        return struct.unpack("<"+fmt, self.read(struct.calcsize("<"+fmt)))[0]

    def count(self, fmt="i", maximum=1000000):
        value = self.number(fmt)
        if not 0 <= value <= maximum:
            raise ValueError("Invalid HIC expected-value count")
        return value

    def text(self):
        value = bytearray()
        while len(value) < 65536:
            char = self.read(1)
            if char == b"\0":
                return value.decode("utf-8")
            value.extend(char)
        raise ValueError("Invalid HIC string")

    def seek(self, position):
        if not 0 <= position <= self.size:
            raise ValueError("Invalid HIC expected-value offset")
        self.stream.seek(position)


@functools.lru_cache(maxsize=4)
def _index(path, mtime, size, inode):
    with open(path, "rb") as stream:
        rd = BinaryReader(stream)
        if rd.read(4) != b"HIC\0":
            raise ValueError("Expected values require a HIC file")
        version = rd.number("i")
        if not 6 <= version <= 9:
            raise ValueError(f"HIC expected values are unsupported for version {version}")
        master = rd.number("q")
        rd.text()  # Assembly identifier.
        modern = version > 8
        if modern:
            rd.read(16)  # Normalization-vector index position and size.
        for _ in range(rd.count()):
            rd.text()
            rd.text()
        chromosomes = []
        for _ in range(rd.count()):
            chromosomes.append(rd.text())
            rd.number("q" if modern else "i")
        rd.seek(master)
        nbytes = rd.count("q" if modern else "i", maximum=rd.size)
        norm_section = stream.tell()+nbytes
        for _ in range(rd.count()):
            rd.text()
            rd.read(12)  # Matrix block position and byte count.
        curves = {}
        dtype = "<f4" if modern else "<f8"
        width = 4 if modern else 8

        def section(normalized=False):
            for _ in range(rd.count()):
                norm = rd.text() if normalized else "NONE"
                unit = rd.text()
                resolution = rd.number("i")
                count = rd.count("q" if modern else "i", maximum=rd.size//width)
                offset = stream.tell()
                rd.seek(offset+count*width)
                factors = {}
                for _ in range(rd.count()):
                    chromosome = rd.number("i")
                    factors[chromosome] = rd.number("f" if modern else "d")
                if unit == "BP" and count:
                    curves[(norm.upper(), resolution)] = (offset, count, dtype, factors)
        section()
        if norm_section < rd.size:
            rd.seek(norm_section)
            section(normalized=True)
        return chromosomes, curves


def expected_index(path):
    path = Path(path).resolve()
    stat = path.stat()
    return _index(str(path), stat.st_mtime_ns, stat.st_size, stat.st_ino)


def expected_matrix(path, x, y, resolution, norm=None):
    if x.chrom.removeprefix("chr") != y.chrom.removeprefix("chr"):
        raise ValueError("Distance-dependent expected values require the same chromosome on both axes")
    chromosomes, curves = expected_index(path)
    key = ((norm or "NONE").upper(), resolution)
    if key not in curves:
        raise ValueError(f"No file expected values for {key[0]} at {resolution:,} bp")
    chromosome = next((i for i, name in enumerate(chromosomes)
                       if name.removeprefix("chr") == x.chrom.removeprefix("chr")), None)
    if chromosome is None:
        raise ValueError("Chromosome is missing from the HIC header")
    offset, count, dtype, factors = curves[key]
    factor = factors.get(chromosome, 1)
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError("Invalid chromosome expected-value normalization factor")
    xs = np.arange(x.start//resolution, (x.end+resolution-1)//resolution)
    ys = np.arange(y.start//resolution, (y.end+resolution-1)//resolution)
    distances = np.minimum(np.abs(ys[:, None]-xs[None, :]), count-1)
    lo, hi = int(distances.min()), int(distances.max())+1
    with open(path, "rb") as stream:
        stream.seek(offset+lo*np.dtype(dtype).itemsize)
        curve = np.frombuffer(stream.read((hi-lo)*np.dtype(dtype).itemsize), dtype=dtype)
    if curve.size != hi-lo:
        raise ValueError("Truncated HIC expected-value curve")
    result = np.asarray(curve[distances-lo]/factor, dtype=np.float32)
    result[~np.isfinite(result) | (result < 0)] = np.nan
    return result
