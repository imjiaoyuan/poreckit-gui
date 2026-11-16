"""Bounded, file-versioned reader/matrix caches owned by a plotting worker."""
from collections import OrderedDict
from pathlib import Path


def source_chromosome(reader, chrom):
    if hasattr(reader, "chromsizes"):
        names = {name for name, length in reader.chromsizes()}
        alternate = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if chrom not in names and alternate in names:
            return alternate
    return chrom


class RenderCache:
    """Bounded, file-versioned reader/matrix caches owned by a plotting worker.

    Matrices, annotation records and sequence payloads share ONE total memory
    budget of ``max_bytes``: ``bytes + record_bytes + payload_bytes`` never
    exceeds ``max_bytes``. When an insert would push the combined total past
    the budget, the oldest entries are evicted across all three caches until
    it fits again. A single item larger than ``max_bytes`` (matrices),
    ``max_bytes // 2`` (records) or ``max_bytes // 4`` (payloads) is refused
    rather than stored. Readers are bounded separately by ``max_readers``.
    """

    def __init__(self, max_bytes=64 * 1024**2, max_readers=2):
        self.max_bytes = max_bytes
        self.max_readers = max_readers
        self.readers = OrderedDict()
        self.matrices = OrderedDict()
        self.metadata_cache = OrderedDict()
        self.records_cache = OrderedDict()
        self.record_bytes = 0
        self.bytes = 0
        self.hits = 0
        self.misses = 0
        self.payloads = OrderedDict()
        self.payload_bytes = 0
        # Global insertion order across the three bounded caches, so eviction
        # is oldest-first over the shared budget instead of per-cache.
        self.order = OrderedDict()
        self.storage_cache = OrderedDict()

    def _touch(self, name, key):
        self.order.move_to_end((name, key))

    def _evict_oldest(self):
        (name, key), _ = self.order.popitem(last=False)
        cache = getattr(self, name)
        if name == "matrices":
            self.bytes -= cache.pop(key).nbytes
        elif name == "records_cache":
            self.record_bytes -= cache.pop(key)[1]
        else:
            self.payload_bytes -= len(cache.pop(key))

    def _store(self, name, key, entry, size):
        """Insert into one bounded cache, evicting the oldest entries across
        all three caches (matrices, records, payloads) so the shared total
        ``bytes + record_bytes + payload_bytes`` stays within ``max_bytes``."""
        while self.order and self.bytes + self.record_bytes + self.payload_bytes + size > self.max_bytes:
            self._evict_oldest()
        getattr(self, name)[key] = entry
        self.order[(name, key)] = None
        if name == "matrices":
            self.bytes += size
        elif name == "records_cache":
            self.record_bytes += size
        else:
            self.payload_bytes += size

    def sequence_payload(self, key, compute):
        """Keep bounded query results so display edits do not rescan alignments."""
        import json
        if key in self.payloads:
            self.payloads.move_to_end(key)
            self._touch("payloads", key)
            return json.loads(self.payloads[key])
        result = compute()
        encoded = json.dumps(result, ensure_ascii=False).encode("utf8")
        if len(encoded) <= self.max_bytes // 4:
            self._store("payloads", key, encoded, len(encoded))
        return result

    def version(self, path):
        path = Path(path).resolve()
        stat = path.stat()
        return str(path), stat.st_mtime_ns, stat.st_size, stat.st_ino

    def reader(self, path):
        from poreckit import make_reader
        version = self.version(path)
        for key in list(self.readers):
            if key[0] == version[0] and key != version:
                self._close(self.readers.pop(key))
        if version not in self.readers:
            reader = make_reader(str(path))
            if hasattr(reader, "cache_bytes"):
                reader.cache_bytes = 16 * 1024**2
            self.readers[version] = reader
        self.readers.move_to_end(version)
        while len(self.readers) > self.max_readers:
            self._close(self.readers.popitem(last=False)[1])
        return self.readers[version]

    def matrix(self, path, chrom, start, end, resolution, norm=None):
        key = (self.version(path), chrom, start, end, resolution, norm)
        def compute():
            reader = self.reader(path)
            return reader.matrix(source_chromosome(reader, chrom), start, end, resolution, norm)
        return self._load(key, compute)

    def signal(self, path, chrom, start, end, bins=800, index=None):
        from .tracks import query_signal
        key = (self.version(path), "signal", chrom, start, end, bins,
               self.version(index) if index else None)
        return self._load(key, lambda: query_signal(path, chrom, start, end, bins, index))

    def metadata(self, path):
        key = self.version(path)
        if key not in self.metadata_cache:
            from .worker import inspect_file
            self.metadata_cache[key] = inspect_file(path)
        self.metadata_cache.move_to_end(key)
        while len(self.metadata_cache) > 4:
            self.metadata_cache.popitem(last=False)
        return self.metadata_cache[key]

    def rectangle(self, path, x, y, resolution, norm=None):
        if x == y:
            return self.matrix(path, x.chrom, x.start, x.end, resolution, norm)
        key = (self.version(path), "rectangle", x, y, resolution, norm)
        def compute():
            reader = self.reader(path)
            return reader.matrix(source_chromosome(reader, y.chrom), y.start, y.end, resolution, norm,
                                 chrom2=source_chromosome(reader, x.chrom), start2=x.start, end2=x.end)
        return self._load(key, compute)

    def symmetric(self, path, resolution):
        """True when the source stores a symmetric matrix, so a tile may be
        shared with its transpose. Cooler `storage-mode='square'` files keep
        both directions independently and must be queried as directed (x, y)."""
        key = (self.version(path), resolution)
        if key not in self.storage_cache:
            name = str(path).lower()
            # .hic only stores the upper triangle, so symmetry is guaranteed;
            # cooler defaults to symmetric-upper unless a group says 'square'.
            symmetric = name.endswith(".hic")
            if name.endswith((".cool", ".mcool", ".h5")):
                import h5py
                try:
                    with h5py.File(path, "r") as source:
                        group = (source["resolutions"][str(resolution)]
                                 if "resolutions" in source else source)
                        symmetric = group.attrs.get("storage-mode", "symmetric-upper") != "square"
                except OSError:
                    symmetric = False
            self.storage_cache[key] = symmetric
            while len(self.storage_cache) > 16:
                self.storage_cache.popitem(last=False)
        return self.storage_cache[key]

    def tiled_rectangle(self, path, x, y, resolution, norm=None, tile_bins=128):
        """Reuse aligned source tiles across overlapping genomic viewports."""
        import math
        import numpy as np
        from .model import Region
        reader = self.reader(path)
        if not hasattr(reader, "chromsizes") or x.start % resolution or y.start % resolution:
            return self.rectangle(path, x, y, resolution, norm)
        lengths = dict(reader.chromsizes())
        def axis(region):
            name = source_chromosome(reader, region.chrom)
            length = lengths[name]
            step = tile_bins * resolution
            return [Region(name, a, min(length, a+step))
                    for a in range(region.start//step*step, region.end, step)]
        key = (self.version(path), "tiled-view", x, y, resolution, norm)
        symmetric = self.symmetric(path, resolution)
        def compute():
            result = np.empty((math.ceil((y.end-y.start)/resolution),
                               math.ceil((x.end-x.start)/resolution)), dtype=np.float32)
            for tx in axis(x):
                for ty in axis(y):
                    # Only symmetric source tiles may share one cache entry.
                    transpose = symmetric and (tx.chrom, tx.start) > (ty.chrom, ty.start)
                    block = self.rectangle(path, ty, tx, resolution, norm).T if transpose else self.rectangle(path, tx, ty, resolution, norm)
                    xa, xb = max(x.start, tx.start), min(x.end, tx.end)
                    ya, yb = max(y.start, ty.start), min(y.end, ty.end)
                    dest_x = slice((xa-x.start)//resolution, math.ceil((xb-x.start)/resolution))
                    dest_y = slice((ya-y.start)//resolution, math.ceil((yb-y.start)/resolution))
                    src_x = slice((xa-tx.start)//resolution, math.ceil((xb-tx.start)/resolution))
                    src_y = slice((ya-ty.start)//resolution, math.ceil((yb-ty.start)/resolution))
                    result[dest_y, dest_x] = block[src_y, src_x]
            return result
        return self._load(key, compute)

    def records(self, kind, path, chrom=None):
        """Reuse public annotation parsers while tracking their Python heap size."""
        from poreckit import parse_genes, parse_tads, parse_loops
        key = self.version(path), kind, chrom
        if key in self.records_cache:
            self.records_cache.move_to_end(key)
            self._touch("records_cache", key)
            return self.records_cache[key][0]
        path = str(path)
        same = lambda name: name.removeprefix("chr") == chrom.removeprefix("chr")
        if kind == "genes":
            values = parse_genes(path, chrom)
            if not values:
                alternate = chrom[3:] if chrom.startswith("chr") else "chr"+chrom
                values = parse_genes(path, alternate)
        elif kind == "tads":
            values = [(chrom, *v[1:]) for v in parse_tads(path) if same(v[0])]
        elif kind == "loops":
            values = list(parse_loops(path))
            if chrom is not None:
                values = [(chrom, *v[1:3], chrom, *v[4:]) for v in values if same(v[0]) and same(v[3])]
        else:
            raise ValueError(f"Unsupported annotation type: {kind}")
        values = tuple(values)
        import sys
        pending, seen, size = [values], set(), 0
        while pending:
            value = pending.pop()
            if id(value) in seen:
                continue
            seen.add(id(value))
            size += sys.getsizeof(value)
            if isinstance(value, (tuple, list)):
                pending.extend(value)
        if size <= self.max_bytes // 2:
            self._store("records_cache", key, (values, size), size)
        return values

    def _load(self, key, compute):
        import numpy as np
        if key in self.matrices:
            self.hits += 1
            self.matrices.move_to_end(key)
            self._touch("matrices", key)
            return self.matrices[key]
        self.misses += 1
        matrix = np.asarray(compute(), dtype=np.float32)
        if matrix.nbytes <= self.max_bytes:
            matrix.setflags(write=False)
            self._store("matrices", key, matrix, matrix.nbytes)
        return matrix

    @staticmethod
    def _close(reader):
        if hasattr(reader, "close"):
            reader.close()

    def close(self):
        for reader in self.readers.values():
            self._close(reader)
        self.readers.clear()
        self.matrices.clear()
        self.metadata_cache.clear()
        self.records_cache.clear()
        self.record_bytes = 0
        self.storage_cache.clear()
        self.bytes = 0
        self.payloads.clear()
        self.payload_bytes = 0
        self.order.clear()
