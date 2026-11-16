"""Explicit, non-executable data references for plotting forms and presets."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def load_value(value, context=None, resources=None, seen=None):
    context = context or {}
    seen = set() if seen is None else seen
    recurse = lambda v: load_value(v, context, resources, seen)
    if isinstance(value, list):
        return [recurse(v) for v in value]
    if not isinstance(value, dict):
        return value
    if "$array" in value:
        return np.asarray(recurse(value["$array"]))
    if "$file" in value:
        path = Path(value["$file"]).expanduser().resolve()
        if path in seen:
            raise ValueError(f"Circular data reference: {path}")
        seen.add(path)
        try:
            if path.suffix.lower() == ".json":
                result = recurse(json.loads(path.read_text(encoding="utf-8")))
            elif path.suffix.lower() == ".npy":
                result = np.load(path, allow_pickle=False)
            elif path.suffix.lower() == ".npz":
                with np.load(path, allow_pickle=False) as archive:
                    result = {k: archive[k] for k in archive.files}
                if "key" not in value and len(result) == 1:
                    result = next(iter(result.values()))
            elif path.suffix.lower() in {".csv", ".tsv", ".txt"}:
                result = np.loadtxt(path, delimiter="," if path.suffix.lower() == ".csv" else None)
            else:
                raise ValueError("Structured data files must be JSON, NPY, NPZ or numeric CSV/TSV/TXT")
            if "key" in value:
                result = result[value["key"]]
            return result
        finally:
            seen.remove(path)
    if "$reader" in value:
        from poreckit import make_reader
        reader = make_reader(value["$reader"])
        if resources is not None and hasattr(reader, "close"):
            resources.callback(reader.close)
        return reader
    if "$contacts" in value:
        from poreckit import read_contacts
        return read_contacts(value["$contacts"])
    if "$matrix" in value:
        from poreckit import make_reader
        from .model import Region
        spec = {**context, **(value["$matrix"] if isinstance(value["$matrix"], dict) else {})}
        if not spec.get("input") or not spec.get("region"):
            raise ValueError("Load contact data in Genome Browser, or provide input and region in the matrix reference")
        region = Region.parse(spec["region"])
        reader = make_reader(spec["input"])
        try:
            from .cache import source_chromosome
            chrom = source_chromosome(reader, region.chrom)
            resolution = spec.get("binsize")
            if resolution is None:
                target = max(1, (region.end - region.start + 249) // 250)
                available = sorted(reader.resolutions())
                resolution = next((n for n in available if n >= target), available[-1] if available else target)
            if resolution <= 0 or (region.end - region.start + resolution - 1) // resolution > 2500:
                raise ValueError("Matrix exceeds 2500 bins; use a smaller region or coarser resolution")
            return reader.matrix(chrom, region.start, region.end, resolution, spec.get("norm") or None)
        finally:
            if hasattr(reader, "close"):
                reader.close()
    return {k: recurse(v) for k, v in value.items()}
