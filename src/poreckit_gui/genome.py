"""Bounded whole-genome contact aggregation and rendering."""
import math
import mmap
import tempfile
from pathlib import Path

import numpy as np

MAX_BINS = 1500
CHUNK = 262144


def genome_layout(chromsizes, resolution):
    result, offset = [], 0
    for name, length in chromsizes:
        count = math.ceil(length/resolution)
        result.append(dict(name=name, length=length, offset=offset, count=count))
        offset += count
    if not result or offset > MAX_BINS:
        raise ValueError(f"Whole-genome view needs 1–{MAX_BINS} bins; select a coarser overview resolution")
    return result, offset


def cooler_overview(path, resolution, layout, source_resolution):
    """Sum stored pixels into genomic bins without materializing a fine matrix."""
    import h5py
    total = layout[-1]["offset"]+layout[-1]["count"]
    result = np.zeros((total, total), dtype=np.float64)
    with h5py.File(path, "r") as source, tempfile.TemporaryFile() as backing:
        g = source[f"resolutions/{source_resolution}"] if "resolutions" in source else source
        names = [v.decode() if isinstance(v, bytes) else str(v) for v in g["chroms/name"][:]]
        mapping = {names.index(c["name"]): c for c in layout}
        nbin = len(g["bins/start"])
        if not nbin:
            raise ValueError("The contact map contains no bins")
        # Disk-backed bin lookup also bounds memory for very fine input files.
        backing.truncate(nbin*4)
        memory = mmap.mmap(backing.fileno(), nbin*4)
        target = np.ndarray(nbin, dtype=np.int32, buffer=memory)
        target[:] = -1
        try:
            for start in range(0, nbin, CHUNK):
                stop = min(start+CHUNK, nbin)
                chrom = g["bins/chrom"][start:stop]
                if chrom.dtype.kind in "SUO":
                    ids = {name: i for i, name in enumerate(names)}
                    chrom = np.array([ids[v.decode() if isinstance(v, bytes) else str(v)] for v in chrom])
                starts = g["bins/start"][start:stop]
                for cid in np.unique(chrom):
                    if cid not in mapping:
                        continue
                    c = mapping[cid]
                    selected = chrom == cid
                    positions = starts[selected]//resolution
                    if np.any(positions < 0) or np.any(positions >= c["count"]):
                        raise ValueError("A source bin lies outside its chromosome")
                    target[start+np.flatnonzero(selected)] = c["offset"]+positions
            storage = g.attrs.get("storage-mode", "symmetric-upper")
            if isinstance(storage, bytes):
                storage = storage.decode()
            if storage not in {"symmetric-upper", "square"}:
                raise ValueError(f"Unsupported contact storage mode: {storage}")
            for start in range(0, len(g["pixels/count"]), CHUNK):
                stop = start+CHUNK
                b1, b2 = g["pixels/bin1_id"][start:stop], g["pixels/bin2_id"][start:stop]
                if np.any(b1 < 0) or np.any(b2 < 0) or np.any(b1 >= nbin) or np.any(b2 >= nbin):
                    raise ValueError("Invalid contact bin index")
                a, b = target[b1], target[b2]
                values = np.asarray(g["pixels/count"][start:stop], dtype=np.float64)
                valid = (a >= 0) & (b >= 0) & np.isfinite(values)
                a, b, values = a[valid], b[valid], values[valid]
                np.add.at(result, (a, b), values)
                if storage == "symmetric-upper":
                    off = a != b
                    np.add.at(result, (b[off], a[off]), values[off])
        finally:
            del target
            memory.close()
    return result


def genome_data(path, chromsizes, resolutions, resolution=None, state=None):
    lengths = [length for name, length in chromsizes]
    if not lengths:
        raise ValueError("This file does not provide chromosome lengths")
    target = max(1, math.ceil(sum(lengths)/600))
    cool = str(path).lower().endswith((".cool", ".mcool"))
    available = sorted(resolutions)
    if resolution is None:
        if cool:
            base = available[-1]
            resolution = next((r for r in available if r >= target), math.ceil(target/base)*base)
        else:
            resolution = next((n for n in available if n >= target), available[-1] if available else target)
    if not isinstance(resolution, int) or resolution <= 0:
        raise ValueError("Overview resolution must be a positive integer")
    if not cool and available and resolution not in available:
        raise ValueError("Select a resolution stored in this contact map")
    layout, n = genome_layout(chromsizes, resolution)
    def compute():
        if cool:
            choices = [r for r in available if resolution % r == 0]
            if not choices:
                raise ValueError("Overview resolution must be a multiple of a stored resolution")
            return cooler_overview(path, resolution, layout, max(choices))
        from poreckit import make_reader
        from .model import Region
        rd = state.reader(path) if state else make_reader(path)
        try:
            result = np.zeros((n, n), dtype=np.float32)
            for i, a in enumerate(layout):
                x = Region(a["name"], 0, a["length"])
                for b in layout[i:]:
                    y = Region(b["name"], 0, b["length"])
                    if state:
                        block = state.rectangle(path, x, y, resolution)
                    elif x == y:
                        block = rd.matrix(x.chrom, 0, x.end, resolution, None)
                    else:
                        block = rd.matrix(y.chrom, 0, y.end, resolution, None, chrom2=x.chrom, start2=0, end2=x.end)
                    if block.shape != (b["count"], a["count"]):
                        raise ValueError("Whole-genome matrix dimensions do not match chromosome bins")
                    xs, ys = slice(a["offset"], a["offset"]+a["count"]), slice(b["offset"], b["offset"]+b["count"])
                    result[ys, xs] = block
                    if x != y:
                        result[xs, ys] = block.T
            return result
        finally:
            if state is None and hasattr(rd, "close"):
                rd.close()
    if state:
        matrix = state._load((state.version(path), "genome", resolution, tuple(map(tuple, chromsizes))), compute)
    else:
        matrix = compute()
    return matrix, layout, resolution


def draw_genome(matrix, layout, resolution, output, cmap="hic_whitered", log=False, name="Whole genome"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from poreckit.layout.style import register_hic_colormaps
    from .mapdata import save_map_data, svg_plot_bounds
    from .model import Region
    register_hic_colormaps()
    n = len(matrix)
    fig, ax = plt.subplots(figsize=(9, 8), dpi=140)
    try:
        finite = matrix[np.isfinite(matrix) & (matrix > 0)]
        vmax = float(np.percentile(finite, 99)) if finite.size else 1
        shown = np.log1p(matrix) if log else matrix
        im = ax.imshow(np.ma.masked_invalid(shown), origin="upper", extent=[0, n, n, 0], interpolation="nearest",
                       cmap=plt.get_cmap(cmap).with_extremes(bad="#e5e5e5"), vmin=0, vmax=np.log1p(vmax) if log else vmax)
        ticks = [c["offset"]+c["count"]/2 for c in layout]
        labels = [c["name"] for c in layout]
        ax.set_xticks(ticks, labels, rotation=90, fontsize=7)
        ax.set_yticks(ticks, labels, fontsize=7)
        for c in layout[1:]:
            ax.axvline(c["offset"], color="#777777", lw=.35)
            ax.axhline(c["offset"], color="#777777", lw=.35)
        ax.set_title(f"{name}\nRaw contacts · {resolution:,} bp / bin")
        bar = fig.colorbar(im, ax=ax, fraction=.025, pad=.02)
        bar.set_label("log1p contact" if log else "contact")
        ax.patch.set_gid("contact-map-bounds")
        fig.savefig(output, bbox_inches="tight", dpi=180)
        result = {"output": output, "binsize": resolution, "chromosomes": layout, "log": log, "cmap": cmap, "name": name}
        if Path(output).suffix.lower() == ".svg":
            spec = save_map_data(output, matrix, matrix, None, Region("genome", 0, n), Region("genome", 0, n), 1)
            spec.update(svg_plot_bounds(output) or {})
            spec.update(chromosomes=layout, genomic_binsize=resolution)
            result["map_data"] = spec
        return result
    finally:
        plt.close(fig)
