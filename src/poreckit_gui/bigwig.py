"""Exact regional bigWig queries from local files or HTTP byte ranges."""
import math
import struct
import zlib

import numpy as np

from .remote import open_binary


def query_bigwig(path, chrom, start, end, bins):
    from .sequenceio import bigwig_chromsizes, chromosome, exact, version
    chrom = chromosome(dict(bigwig_chromsizes(version(path))), chrom)
    with open_binary(path) as source:
        head = exact(source, 64)
        endian = '<' if struct.unpack_from('<I', head)[0] == 0x888FFC26 else '>'
        tree, _, index = struct.unpack_from(endian+'QQQ', head, 8)
        compressed = struct.unpack_from(endian+'I', head, 52)[0]
        source.seek(tree)
        _, _, key_size, _, _, _ = struct.unpack(endian+'IIIIQQ', exact(source, 32))
        pending, visited, cid = [tree+32], set(), None
        while pending:
            node = pending.pop()
            if node in visited:
                raise ValueError('Cyclic bigWig chromosome tree')
            visited.add(node)
            source.seek(node)
            leaf, _, n = struct.unpack(endian+'BBH', exact(source, 4))
            for _ in range(n):
                key = exact(source, key_size).rstrip(b'\0').decode()
                value = exact(source, 8)
                if leaf and key == chrom:
                    cid = struct.unpack(endian+'II', value)[0]
                elif not leaf:
                    pending.append(struct.unpack(endian+'Q', value)[0])
        if cid is None:
            raise ValueError('Chromosome absent from bigWig')
        source.seek(index)
        header = exact(source, 48)
        if struct.unpack_from(endian+'I', header)[0] != 0x2468ACE0:
            raise ValueError('Invalid bigWig spatial index')
        pending, visited, blocks = [index+48], set(), set()
        while pending:
            node = pending.pop()
            if node in visited:
                raise ValueError('Cyclic bigWig spatial index')
            visited.add(node)
            source.seek(node)
            leaf, _, n = struct.unpack(endian+'BBH', exact(source, 4))
            for _ in range(n):
                sc, sb, ec, eb = struct.unpack(endian+'IIII', exact(source, 16))
                offset = struct.unpack(endian+'Q', exact(source, 8))[0]
                size = struct.unpack(endian+'Q', exact(source, 8))[0] if leaf else 0
                if (sc < cid or sc == cid and sb < end) and (ec > cid or ec == cid and eb > start):
                    if leaf:
                        blocks.add((offset, size))
                    else:
                        pending.append(offset)
        total, covered = np.zeros(bins), np.zeros(bins)
        width = (end-start)/bins
        for offset, size in sorted(blocks):
            source.seek(offset)
            data = exact(source, size)
            if compressed:
                if compressed>64*1024*1024:
                    raise ValueError('Oversized bigWig decompression buffer')
                decoder=zlib.decompressobj()
                data=decoder.decompress(data,64*1024*1024+1)
                if len(data)>64*1024*1024 or not decoder.eof:
                    raise ValueError('Invalid compressed bigWig data block')
            at = 0
            while at < len(data):
                section_chrom, origin, _, step, span, typ, _, count = struct.unpack_from(endian+'IIIIIBBH', data, at)
                at += 24
                for i in range(count):
                    if typ == 1:
                        a, b, value = struct.unpack_from(endian+'IIf', data, at)
                        at += 12
                    elif typ == 2:
                        a, value = struct.unpack_from(endian+'If', data, at)
                        b = a+span
                        at += 8
                    elif typ == 3:
                        value = struct.unpack_from(endian+'f', data, at)[0]
                        a, b = origin+i*step, origin+i*step+span
                        at += 4
                    else:
                        raise ValueError('Unsupported bigWig data section')
                    a, b = max(a, start), min(b, end)
                    if section_chrom != cid or b <= a or not math.isfinite(value):
                        continue
                    lo, hi = max(0, int((a-start)/width)), min(bins, math.ceil((b-start)/width))
                    positions = np.arange(lo, hi)
                    overlap = np.minimum(b, start+(positions+1)*width)-np.maximum(a, start+positions*width)
                    total[lo:hi] += overlap*value
                    covered[lo:hi] += overlap
        return np.divide(total, covered, out=np.full(bins,np.nan), where=covered>0).astype(np.float32)
