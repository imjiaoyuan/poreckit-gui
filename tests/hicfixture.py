"""Minimal .hic builder mirroring poreckit's Rust test fixtures (hic.rs)."""
import struct
import zlib


def wstr(buf, s):
    buf.extend(s.encode() + b"\0")


def dense_block(vals, w):
    body = struct.pack('<i', 1) + struct.pack('<i', 0) + struct.pack('<i', 0)
    body += bytes([0, 0, 0, 2])  # useShort, useShortBinX, useShortBinY, type=dense
    body += struct.pack('<i', len(vals)) + struct.pack('<h', w)
    body += struct.pack(f'<{len(vals)}h', *vals)
    return zlib.compress(body)


def dense_block_v8(vals, w):
    body = struct.pack('<i', 1) + struct.pack('<i', 0) + struct.pack('<i', 0)
    body += bytes([0, 2])  # useShort, type=dense
    body += struct.pack('<i', len(vals)) + struct.pack('<h', w)
    body += struct.pack(f'<{len(vals)}h', *vals)
    return zlib.compress(body)


def mzd(c1, c2, block_pos, block_len, binsize=1, bbc=2, bcc=2):
    b = bytearray(struct.pack('<ii', c1, c2) + struct.pack('<i', 1))
    wstr(b, 'BP')
    b += struct.pack('<i', 0) + struct.pack('<4f', 0, 0, 0, 0)
    b += struct.pack('<iii', binsize, bbc, bcc) + struct.pack('<i', 1)
    b += struct.pack('<i', 0) + struct.pack('<q', block_pos) + struct.pack('<i', block_len)
    return bytes(b)


def build(w1=(1.0, 1.0), w2=(1.0, 1.0), blk_intra=None, blk_inter=None,
          version=9, norm_names=('VC', 'VC')):
    """v9 layout: header, footer, NVI norm index, MZDs, blocks, f32 vectors.
    v8 layout: header, footer (master entries, two empty expected layers, norm
    combos, blob index), MZDs, blocks, f64 vectors."""
    modern = version > 8
    dense = dense_block if version > 8 else dense_block_v8
    blk_intra = blk_intra or dense([5, 0, 0, 5], 2)
    blk_inter = blk_inter or dense([3, 0, 4, 7], 2)
    mzd_len = 67
    header_len = (4 + 4 + 8 + 5 + (16 if modern else 0) + 4 + 4
                  + 2 * (5 + (8 if modern else 4)))
    entry_len = 4 + 2 * (4 + 8 + 4)
    nbytes_len = 8 if modern else 4

    def norm_entry(name, chridx, pos, size):
        b = bytearray()
        wstr(b, name)
        b += struct.pack('<i', chridx)
        wstr(b, 'BP')
        b += struct.pack('<i', 1)
        b += struct.pack('<q', pos)
        b += struct.pack('<q' if modern else '<i', size)
        return b

    if modern:
        tail = bytearray(struct.pack('<i', 2))
        tail += norm_entry(norm_names[0], 0, 0, 14)
        tail += norm_entry(norm_names[1], 1, 0, 14)
    else:
        tail = bytearray()
        tail += struct.pack('<i', 0)  # expected-value layer: zero entries
        tail += struct.pack('<i', 2)  # norm combos (no positions here)
        for name in norm_names:
            wstr(tail, name)
            wstr(tail, 'BP')
            tail += struct.pack('<iii', 1, 0, 0)
        tail += struct.pack('<i', 2)  # vector-blob index
        tail += norm_entry(norm_names[0], 0, 0, 14)
        tail += norm_entry(norm_names[1], 1, 0, 14)

    intra_pos = header_len + nbytes_len + entry_len + len(tail)
    inter_pos = intra_pos + mzd_len + len(blk_intra)
    vec_size = 16 if modern else 20
    vec1_pos = inter_pos + mzd_len + len(blk_inter)
    vec2_pos = vec1_pos + vec_size

    # patch the norm-index positions now that the layout is known
    if modern:
        at = 4
        for i, pos in enumerate((vec1_pos, vec2_pos)):
            at += len(norm_names[i]) + 1 + 4 + 3 + 4
            struct.pack_into('<q', tail, at, pos)
            at += 8 + 8
    else:
        at = 4  # empty expected layer
        at += 4  # combo count
        for name in norm_names:
            at += len(name) + 1 + 3 + 4 + 4 + 4
        at += 4  # blob count
        for i, pos in enumerate((vec1_pos, vec2_pos)):
            at += len(norm_names[i]) + 1 + 4 + 3 + 4
            struct.pack_into('<q', tail, at, pos)
            at += 8 + 4

    chrom_len_fmt = '<q' if modern else '<i'
    buf = bytearray()
    buf += b'HIC\0' + struct.pack('<i', version) + struct.pack('<q', header_len)
    wstr(buf, 'hg19')
    if modern:
        buf += struct.pack('<qq', 0, 0)  # nviPosition, nviLength
    buf += struct.pack('<i', 0) + struct.pack('<i', 2)
    for name in ('chr1', 'chr2'):
        wstr(buf, name)
        buf += struct.pack(chrom_len_fmt, 2)
    assert len(buf) == header_len
    footer_len = nbytes_len + entry_len + len(tail)
    buf += struct.pack('<q' if modern else '<i', footer_len)
    buf += struct.pack('<i', 2)
    for key, pos in (('0_0', intra_pos), ('0_1', inter_pos)):
        wstr(buf, key)
        buf += struct.pack('<q', pos) + struct.pack('<i', 0)
    assert len(buf) == header_len + nbytes_len + entry_len
    buf += tail
    assert len(buf) == intra_pos, (len(buf), intra_pos)
    buf += mzd(0, 0, intra_pos + mzd_len, len(blk_intra))
    buf += blk_intra
    buf += mzd(0, 1, inter_pos + mzd_len, len(blk_inter))
    buf += blk_inter
    for w in (w1, w2):
        if modern:
            buf += struct.pack('<q', 2) + struct.pack('<2f', *w)
        else:
            buf += struct.pack('<i', 2) + struct.pack('<2d', *w)
    return bytes(buf)
