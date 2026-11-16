"""HTTP byte-range reader for remote `.hic` contact maps.

Mirrors the binary layout of poreckit's Rust `HicReader`
(`src/io/hic.rs`): a 4-byte magic, version and master-index position in the
head, then the master-index footer (chromosome table lives in the head,
per-pair MatrixZoomData entries in the footer), and on demand individual
zlib-compressed block byte ranges. Everything is read through
:mod:`poreckit_gui.remote`, so ranged reads share its bounded 64 MiB block
cache, ETag/Last-Modified validation and eviction discipline; parsed index
layouts are memoised per file identity with the same bounded discipline as
`expected._index`. Dense assembly reuses poreckit's public
`build_dense_matrix` / `build_dense_matrix_rect`, exactly like the Rust
reader, so results are identical to local `.hic` queries.

The reader exposes the same surface the GUI uses for local contact maps:
`chromsizes()`, `resolutions()`, `meta()` and
`matrix(chrom, start, end, bin_size, norm, chrom2=, start2=, end2=)` with
zero-based half-open coordinates.
"""
from __future__ import annotations

import functools
import io
import threading
import zlib

import numpy as np

from .expected import BinaryReader
from .remote import identity, is_remote, open_binary

# Bound on simultaneous HTTP range requests across all remote .hic readers.
RANGE_SLOTS = threading.BoundedSemaphore(8)


def _bounded(method):
    """Serialise ranged fetches through a global semaphore."""
    def wrapper(self, *args):
        with RANGE_SLOTS:
            return method(self, *args)
    return wrapper


class _Cursor:
    """BinaryReader over a seekable remote (or local) binary stream."""

    def __init__(self, path):
        self._raw = open_binary(path)
        self.reader = BinaryReader(self._raw)

    def close(self):
        self._raw.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    @_bounded
    def number(self, fmt):
        return self.reader.number(fmt)

    @_bounded
    def count(self, fmt="i", maximum=1000000):
        return self.reader.count(fmt, maximum)

    @_bounded
    def text(self):
        return self.reader.text()

    @_bounded
    def read(self, n):
        return self.reader.read(n)

    def seek(self, position):
        self.reader.seek(position)

    def tell(self):
        return self._raw.tell()


def _plausible(value, limit, message):
    if not 0 <= value <= limit:
        raise ValueError(f"Implausible {message} — corrupt .hic?")
    return value


def _norm_candidates(norm):
    """Alias resolution mirroring poreckit's `norm_candidates`."""
    if norm in {"NONE", "none", "raw"}:
        return []
    if norm in {"balanced", "BALANCED", "ICE", "ice", "weight"}:
        return ["GW_SCALE", "GW_ICE", "ICE", "KR", "VC_SQRT", "VC"]
    return [norm]


class RemoteHicReader:
    """Remote `.hic` contact-map reader over HTTP(S) byte ranges."""

    def __init__(self, path):
        if not is_remote(path) or not str(path).lower().endswith(".hic"):
            raise ValueError(f"RemoteHicReader needs a remote .hic URL, got {path!r}")
        self.path = str(path)
        self._layout(identity(self.path))  # fail fast on bad magic or missing range support

    # The parsed master index and chromosome table are memoised per file
    # identity with a small bounded LRU, exactly like expected._index.
    @functools.lru_cache(maxsize=4)  # noqa: B019 -- per-instance bounded LRU on an intentional long-lived reader; a manual dict would change no behavior but this keeps the cache bounds explicit
    def _layout(self, ident):
        return _parse_layout(ident[0] if isinstance(ident, tuple) else ident)

    def close(self):
        pass  # streams are opened per query; nothing persistent to release

    def chromsizes(self, bin_size=None):
        """`(name, length)` chromosome table in file order (resolution-free)."""
        layout = self._layout(identity(self.path))
        return [(name, length) for name, length, _ in layout.chroms]

    def resolutions(self):
        return list(self._layout(identity(self.path)).resolutions)

    def meta(self):
        """Container facts as a plain dict, mirroring poreckit's `HicMeta`."""
        layout = self._layout(identity(self.path))
        return {
            "format": ".hic",
            "version": layout.version,
            "assembly": layout.genome_id,
            "attributes": dict(layout.attributes),
            "chromsizes": [(name, length) for name, length, _ in layout.chroms],
            "resolutions": list(layout.resolutions),
            "resolutions_all_match": layout.resolutions_all_match,
            "block_sizes": list(layout.block_sizes),
            "n_chrom_pairs": layout.n_chrom_pairs,
            "n_matrix_zoom_data": layout.n_matrix_zoom_data,
            "nvi": {"position": layout.nvi_position, "present": layout.norm_idx_present},
            "norms": {name: list(res) for name, res in layout.norms.items()},
        }

    def matrix(self, chrom, start, end, bin_size, norm=None,
               chrom2=None, start2=None, end2=None):
        """Dense `(n, n)` float32 intra-chromosomal matrix, or `(n1, n2)`
        inter-chromosomal when `chrom2` is given."""
        start2 = start if start2 is None else start2
        end2 = end if end2 is None else end2
        result = self._query(chrom, start, end, bin_size, norm,
                             None if chrom2 is None else chrom2, start2, end2)
        if isinstance(result, _DenseResult):
            return result.matrix
        from poreckit import build_dense_matrix, build_dense_matrix_rect
        if chrom2 is None:
            bin1, bin2, counts, weights, sel = result
            return build_dense_matrix(bin1, bin2, counts, sel, weights)
        bin1, bin2, counts, wa, wb, sel_a, sel_b, swapped = result
        if wa is not None and wb is not None:
            weights = np.concatenate([wa, wb])
            offset = len(wa)
        else:
            weights, offset = None, 0
        sel_bs = sel_b + offset
        if swapped:
            return build_dense_matrix_rect(bin2 + offset, bin1, counts, sel_bs, sel_a, weights)
        return build_dense_matrix_rect(bin1, bin2 + offset, counts, sel_a, sel_bs, weights)

    # ---- query core, mirroring hic.rs HicReader::query ----

    def _query(self, chrom, start, end, bin_size, norm, chrom2=None, start2=None, end2=None):
        if not bin_size:
            raise ValueError("bin_size must be > 0")
        if end <= start or (chrom2 is not None and end2 is not None and end2 <= start2):
            raise ValueError(f"empty .hic query range {chrom}:{start}-{end}")
        layout = self._layout(identity(self.path))
        names = {name: idx for name, _, idx in layout.chroms}
        if chrom not in names:
            raise ValueError(f"chrom {chrom!r} not in .hic")
        cidx = names[chrom]
        if chrom2 is None:
            ca, cb, swapped, s_a, e_a, s_b, e_b = cidx, cidx, False, start, end, start, end
        else:
            if chrom2 not in names:
                raise ValueError(f"chrom {chrom2!r} not in .hic")
            cidx2 = names[chrom2]
            if cidx2 == cidx:
                raise ValueError("same-chromosome query: omit chrom2 (intra path)")
            ca, cb = min(cidx, cidx2), max(cidx, cidx2)
            swapped = cidx > cidx2
            s_a, e_a = (start2, end2) if swapped else (start, end)
            s_b, e_b = (start, end) if swapped else (start2, end2)

        pair_pos = layout.pair_pos.get(f"{ca}_{cb}")
        if pair_pos is None:
            # The master index only lists chrom pairs that carry blocks; a
            # valid pair without one simply has no contacts — answer zeros.
            na = -(-e_a // bin_size) - s_a // bin_size
            nb = -(-e_b // bin_size) - s_b // bin_size
            if chrom2 is None:
                shape = (na, na)
            elif swapped:
                shape = (nb, na)
            else:
                shape = (na, nb)
            return _DenseResult(np.zeros(shape, dtype=np.float32))

        wanted = _norm_candidates(norm) if norm else []
        norm_positions = self._norm_vectors(layout, wanted, chrom2 is not None,
                                            ca, cb, bin_size)

        blocks, block_bin_count, block_col_count = self._zoom_blocks(layout, pair_pos, bin_size)
        lo_a, hi_a = s_a // bin_size, -(-e_a // bin_size)
        lo_b, hi_b = s_b // bin_size, -(-e_b // bin_size)

        bin1, bin2, counts = [], [], []
        for number in _wanted_blocks(layout.version, chrom2 is not None,
                                     block_bin_count, block_col_count,
                                     lo_a, hi_a, lo_b, hi_b):
            entry = blocks.get(number)
            if entry is not None:
                more = _decode_block(self, entry[0], entry[1], layout.version)
                bin1.extend(more[0])
                bin2.extend(more[1])
                counts.extend(more[2])
        bin1_arr = np.asarray(bin1, dtype=np.int64)
        bin2_arr = np.asarray(bin2, dtype=np.int64)
        counts_arr = np.asarray(counts, dtype=np.float32)
        sel_a = np.arange(lo_a, hi_a, dtype=np.int64)
        sel_b = np.arange(lo_b, hi_b, dtype=np.int64)
        weights = [_read_norm_vector(self, pos, layout.version) for pos in norm_positions]
        if chrom2 is None:
            return bin1_arr, bin2_arr, counts_arr, weights[0], sel_a
        return bin1_arr, bin2_arr, counts_arr, weights[0], weights[1], sel_a, sel_b, swapped

    def _norm_vectors(self, layout, wanted, has_rect, ca, cb, bin_size):
        """Per-side (filePosition, sizeInBytes) for the requested norm, or
        [None] / [None, None] for raw counts. A candidate only wins when every
        queried chromosome carries it."""
        if not wanted:
            return [None] * (2 if has_rect else 1)
        sides = [min(ca, cb)] if not has_rect else [ca, cb]
        found: list[tuple[int, int] | None] = [None] * len(sides)
        for candidate in wanted:
            trial: list[tuple[int, int] | None] = [None] * len(sides)
            for entry in layout.norm_entries:
                nt, chridx, unit, resolution, fpos, size = entry
                if nt != candidate or unit != "BP" or resolution != bin_size:
                    continue
                for i, side in enumerate(sides):
                    if chridx == side:
                        trial[i] = (fpos, size)
            if all(trial):
                found = trial
                break
        if not all(found):
            available = sorted({nt for nt, _, unit, res, _, _ in layout.norm_entries
                                if unit == "BP" and res == bin_size})
            hint = ("no normalization vectors stored at this resolution" if not available
                    else f"available: {', '.join(available)} (use NONE/raw for observed counts)")
            raise ValueError(f"norm {wanted[0]} not found at {bin_size}bp — {hint}")
        return found

    def _zoom_blocks(self, layout, pair_pos, bin_size):
        """Block index of the BP `bin_size` zoom level of one chrom pair."""
        with _Cursor(self.path) as cursor:
            cursor.seek(pair_pos)
            cursor.count()  # c1
            cursor.count()  # c2
            nres = _plausible(cursor.count(), 64, "zoom level count")
            blocks, block_bin_count, block_col_count = {}, 0, 0
            for _ in range(nres):
                unit = cursor.text()
                cursor.number("i")  # zoomIndex
                cursor.number("4f")  # sumCounts, occupiedCell, stdDev, percent95
                binsize = cursor.number("i")
                bbc = cursor.number("i")
                bcc = cursor.number("i")
                nblocks = _plausible(cursor.count(), 1 << 24, "block count")
                is_target = unit == "BP" and binsize == bin_size
                for _ in range(nblocks):
                    bn = cursor.number("i")
                    fpos = cursor.number("q")
                    size = cursor.number("i")
                    if is_target:
                        blocks[bn] = (fpos, size)
                if is_target:
                    block_bin_count, block_col_count = bbc, bcc
        if not block_bin_count:
            raise ValueError(f"BP resolution {bin_size} not in .hic")
        return blocks, block_bin_count, block_col_count


class _DenseResult:
    """Pre-built dense result for pairs without stored blocks."""

    def __init__(self, matrix):
        self.matrix = matrix


def _parse_layout(path):
    """Header + master index + per-pair zoom headers + norm index, in one
    bounded pass over ranged reads."""
    with _Cursor(path) as cursor:
        if cursor.read(4) != b"HIC\0":
            raise ValueError("Not a .hic file (bad magic or truncated)")
        version = cursor.number("i")
        if not 6 <= version <= 9:
            raise ValueError(f"Unsupported .hic version {version}")
        modern = version > 8
        master_pos = cursor.number("q")
        genome_id = cursor.text()
        nvi_position = 0
        if modern:
            nvi_position = cursor.number("q")
            cursor.number("q")  # nviLength
        attributes = []
        for _ in range(_plausible(cursor.count(), 1 << 16, "attribute count")):
            attributes.append((cursor.text(), cursor.text()))
        chroms = []
        for index in range(_plausible(cursor.count(), 65536, "chromosome count")):
            name = cursor.text()
            length = cursor.number("q") if modern else cursor.number("i")
            chroms.append((name, int(length), index))

        # ---- master index footer ----
        cursor.seek(master_pos)
        cursor.number("q" if modern else "i")  # nBytes
        pair_pos = {}
        for _ in range(_plausible(cursor.count(), 1 << 20, "master-index size")):
            key = cursor.text()
            fpos = cursor.number("q")
            cursor.number("i")  # sizeInBytes
            parts = key.split("_")
            if len(parts) == 2 and all(part.lstrip("-").isdigit() for part in parts):
                pair_pos[key] = fpos
        after_entries = cursor.tell()

        # ---- per-pair MatrixZoomData headers: resolutions & block sizes ----
        n_matrix_zoom_data = 0
        resolutions: list[int] = []
        resolutions_all_match = True
        block_sizes: list[tuple[int, int, int]] = []
        for pos in pair_pos.values():
            cursor.seek(pos)
            cursor.count()
            cursor.count()  # c1, c2
            nres = _plausible(cursor.count(), 64, "zoom level count")
            n_matrix_zoom_data += nres
            res = set()
            for _ in range(nres):
                unit = cursor.text()
                cursor.number("i")
                cursor.number("4f")  # zoomIndex + stats
                binsize = cursor.number("i")
                bbc = cursor.number("i")
                bcc = cursor.number("i")
                for _ in range(_plausible(cursor.count(), 1 << 24, "block count")):
                    cursor.number("i")
                    cursor.number("q")
                    cursor.number("i")
                if unit == "BP":
                    res.add(binsize)
                    if not any(row[0] == binsize for row in block_sizes):
                        block_sizes.append((binsize, bbc, bcc))
            if not resolutions:
                resolutions = sorted(res)
            elif sorted(res) != resolutions:
                resolutions_all_match = False
            resolutions = sorted(set(resolutions) | res)

        # ---- norm index: v9 NVI position, v8 footer tail, v9 inline fallback ----
        norm_entries = []  # (type, chrIdx, unit, resolution, filePosition, size)
        norm_sources = []  # (type, unit, resolution) pairs feeding meta().norms
        norm_idx_present = False
        if nvi_position:
            norm_idx_present = True
            cursor.seek(nvi_position)
            for _ in range(_plausible(cursor.count(), 1 << 22, "norm-index size")):
                entry = _norm_entry(cursor, modern)
                norm_entries.append(entry)
                norm_sources.append((entry[0], entry[2], entry[3]))
        elif not modern:
            norm_idx_present = True
            cursor.seek(after_entries)
            for _ in range(_plausible(cursor.count(), 1 << 22, "expected-value size")):
                cursor.text()  # unit
                cursor.number("i")  # resolution
                for _ in range(_plausible(cursor.count(), 1 << 24, "expected count")):
                    cursor.number("d")
                for _ in range(_plausible(cursor.count(), 1 << 24, "expected pairs")):
                    cursor.number("id")
            for _ in range(_plausible(cursor.count(), 1 << 22, "norm-combo count")):
                nt, unit = cursor.text(), cursor.text()
                resolution = cursor.number("i")
                norm_sources.append((nt, unit, resolution))
                for _ in range(_plausible(cursor.count(), 1 << 24, "norm vector")):
                    cursor.number("d")
                for _ in range(_plausible(cursor.count(), 1 << 24, "norm map")):
                    cursor.number("id")
            for _ in range(_plausible(cursor.count(), 1 << 22, "norm-vector count")):
                norm_entries.append(_norm_entry(cursor, modern))
        else:
            # v9 without NVI: the norm index follows the master entries inline;
            # poreckit's meta reports no norms or NVI in that case, queries
            # still resolve through the inline index.
            cursor.seek(after_entries)
            for _ in range(_plausible(cursor.count(), 1 << 22, "norm-index size")):
                norm_entries.append(_norm_entry(cursor, modern))

    norms: dict = {}
    for nt, unit, resolution in norm_sources:
        if unit == "BP":
            values = norms.setdefault(nt, [])
            if resolution not in values:
                values.append(resolution)
    for values in norms.values():
        values.sort()
    layout = _Layout(version=version, genome_id=genome_id, attributes=attributes,
                     chroms=chroms, pair_pos=dict(sorted(pair_pos.items())),
                     resolutions=sorted(resolutions),
                     resolutions_all_match=resolutions_all_match,
                     block_sizes=block_sizes, n_chrom_pairs=len(pair_pos),
                     n_matrix_zoom_data=n_matrix_zoom_data,
                     nvi_position=nvi_position or None,
                     norm_idx_present=norm_idx_present, norms=norms,
                     norm_entries=norm_entries)
    return layout


def _norm_entry(cursor, modern):
    nt = cursor.text()
    chridx = cursor.number("i")
    unit = cursor.text()
    resolution = cursor.number("i")
    fpos = cursor.number("q")
    size = cursor.number("q" if modern else "i")
    return nt, chridx, unit, resolution, fpos, size


class _Layout:
    def __init__(self, **fields):
        self.__dict__.update(fields)


def _read_norm_vector(owner, position, version):
    """Norm weights as a divisive vector of reciprocals (NaN stays NaN)."""
    if position is None:
        return None
    fpos, _size = position
    with _Cursor(owner.path) as cursor:
        cursor.seek(fpos)
        nvalues = _plausible(cursor.number("q" if version > 8 else "i"),
                             1 << 24, "norm vector length")
        width = 4 if version > 8 else 8
        raw = cursor.read(nvalues * width)
    weights = np.frombuffer(raw, dtype="<f4" if version > 8 else "<f8").astype(np.float32)
    nonzero = weights != 0
    weights[nonzero] = 1.0 / weights[nonzero]
    return weights


def _wanted_blocks(version, is_rect, bbc, bcc, lo_a, hi_a, lo_b, hi_b):
    """Block numbers covering the query window: rectangular grid for v8 and
    all inter-chromosomal queries, diagonal-depth blocks for v9 intra."""
    if not (version > 8 and not is_rect):
        return [br * bcc + bc
                for br in range(max(lo_a, 0) // bbc, max(hi_a - 1, 0) // bbc + 1)
                for bc in range(max(lo_b, 0) // bbc, max(hi_b - 1, 0) // bbc + 1)]
    step = max(bbc // 2, 1)
    from math import log, sqrt
    def distance(x, y):
        n = abs(x - y) / sqrt(2.0) / bbc
        return int(log(1.0 + n) / log(2.0))
    wanted: set = set()
    b1 = lo_a
    while b1 < hi_a:
        b2 = min(b1 + step, hi_a)
        b3 = max(b1, lo_b)
        while b3 < hi_b:
            b4 = min(b3 + step, hi_b)
            _push_depth_range(wanted, b1, b2, b3, b4, bbc, bcc, distance)
            b3 = b4
        if b3 >= hi_b and hi_b > lo_b:
            _push_depth_range(wanted, b1, b2, hi_b - 1, hi_b, bbc, bcc, distance)
        b1 = b2
    return sorted(wanted)


def _push_depth_range(wanted, b1, b2, b3, b4, bbc, bcc, distance):
    pad_lo = (b1 + b3) // 2 // bbc
    pad_hi = (b2 + b4) // 2 // bbc + 1
    d_near = distance(b1, max(b4 - 1, 0))
    d_far = distance(max(b2 - 1, 0), b3)
    includes_diag = b1 + 1 > max(b4 - 1, 0) or b3 + 1 > max(b2 - 1, 0)
    near = 0 if includes_diag else min(d_near, d_far)
    far = max(d_near, d_far) + 1
    for pad in range(pad_lo, pad_hi + 1):
        for depth in range(near, far + 1):
            wanted.add(depth * bcc + pad)


def _decode_block(owner, fpos, size, version):
    """One zlib-compressed block body → (bin1, bin2, count) columns."""
    with _Cursor(owner.path) as cursor:
        cursor.seek(fpos)
        raw = cursor.read(size)
    try:
        body = zlib.decompress(raw)
    except zlib.error as exc:
        raise ValueError(f"zlib inflate failed for block at {fpos}: {exc}") from None
    stream = BinaryReader(io.BytesIO(body))
    nrecords = stream.number("i")
    if version < 7:
        bin1, bin2, counts = [], [], []
        for _ in range(nrecords):
            bin1.append(stream.number("i"))
            bin2.append(stream.number("i"))
            counts.append(stream.number("f"))
        return bin1, bin2, counts
    bin_x_offset = stream.number("i")
    bin_y_offset = stream.number("i")
    use_short = stream.number("B")
    if version > 8:
        use_short_x = stream.number("B")
        use_short_y = stream.number("B")
    else:
        use_short_x = use_short_y = use_short
    typ = stream.number("B")
    bin1, bin2, counts = [], [], []
    if typ == 1:
        # Offsets are base + delta (never accumulated). Widths follow the
        # per-axis flags on v9; v7/v8 have no flags and are always u16.
        rows_wide = version > 8 and use_short_y != 0
        cols_wide = version > 8 and use_short_x != 0
        row_count = stream.number("I" if rows_wide else "H")
        for _ in range(row_count):
            dy = stream.number("I" if rows_wide else "H")
            bin_y = bin_y_offset + dy
            col_count = stream.number("I" if cols_wide else "H")
            for _ in range(col_count):
                dx = stream.number("I" if cols_wide else "H")
                bin1.append(bin_x_offset + dx)
                bin2.append(bin_y)
                counts.append(_read_value(stream, use_short))
    elif typ == 2:
        # Dense row-major values (row = binY, col = binX); NaN (-32768 in
        # short mode) marks empty cells and must not become pixels.
        n_pts = stream.number("i")
        width = max(stream.number("h"), 1)
        for i in range(n_pts):
            value = _read_value(stream, use_short)
            if value != value:  # NaN
                continue
            row, col = divmod(i, width)
            bin1.append(bin_x_offset + col)
            bin2.append(bin_y_offset + row)
            counts.append(value)
    else:
        raise ValueError(f"Unknown .hic block record type {typ} at block offset {fpos}")
    return bin1, bin2, counts


def _read_value(stream, use_short):
    if use_short == 0:
        value = stream.number("h")
        return float("nan") if value == -32768 else float(value)
    return stream.number("f")
