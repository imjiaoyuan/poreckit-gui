"""Local, bounded genomic queries. Coordinates are zero-based and half-open.

BAM/BAI, BGZF/TBI, FASTA/FAI and VCF are read without optional Python
packages. CRAM and BCF use samtools/bcftools, with argument arrays only.
"""
from __future__ import annotations

import gzip
import json
import math
import shutil
import struct
import subprocess
import zlib
from functools import lru_cache
from pathlib import Path

import numpy as np

from .remote import exists, identity, is_remote, open_binary, sidecar, source_path

MAX_RECORDS = 100_000
BASES = '=ACMGRSVTWYHKDBN'
CIGAR = 'MIDNSHP=XB'
# BGZF blocks decompress to at most 64 KiB by spec; anything larger means a
# corrupt or hostile file. Queries also get a cumulative decompressed-read cap.
MAX_BGZF_BLOCK = 64 * 1024
MAX_BGZF_QUERY_BYTES = 256 * 1024**2
# Index allocations are capped so a corrupt header cannot imply enormous reads.
MAX_INDEX_REFS = 1_000_000
MAX_INDEX_ENTRIES = 10_000_000  # bins + chunks combined


def exact(stream, size):
    if size < 0 or size > 256 * 1024 * 1024:
        raise ValueError('Invalid genomic record size')
    data = stream.read(size)
    if len(data) != size:
        raise ValueError('Truncated genomic file')
    return data


def version(path):
    return identity(path)


def chromosome(names, requested):
    if requested in names:
        return requested
    other = requested[3:] if requested.startswith('chr') else 'chr' + requested
    if other in names:
        return other
    raise ValueError(f'Chromosome {requested!r} is absent from this file')


class BGZF:
    def __init__(self, path):
        self.name = str(path)
        self.file = open_binary(path)
        self.block = b''
        self.consumed = 0
        self.offset = self.next = self.pos = 0
        self.seek(0)

    def close(self):
        self.file.close()

    def seek(self, virtual):
        offset, position = virtual >> 16, virtual & 65535
        self.file.seek(offset)
        head = self.file.read(12)
        if not head:
            self.block = b''
            self.offset = self.next = offset
            self.pos = 0
            return
        if head[:4] != b'\x1f\x8b\x08\x04':
            raise ValueError('An indexed BAM/VCF must use BGZF compression')
        extra = exact(self.file, struct.unpack_from('<H', head, 10)[0])
        at, size = 0, None
        while at + 4 <= len(extra):
            key, length = extra[at:at+2], struct.unpack_from('<H', extra, at+2)[0]
            if key == b'BC' and length == 2:
                size = struct.unpack_from('<H', extra, at+4)[0] + 1
            at += 4 + length
        if size is None or size < 12 + len(extra) + 8:
            raise ValueError('Invalid BGZF block')
        try:
            payload = head + extra + exact(self.file, size-12-len(extra))
            self.block = gzip.decompress(payload)
        except (ValueError, OSError, EOFError, zlib.error) as exc:
            raise ValueError(f'{self.name}: corrupt or truncated BGZF block ({exc})') from exc
        if len(self.block) > MAX_BGZF_BLOCK:
            raise ValueError(
                f'{self.name}: BGZF block decompresses to {len(self.block)} bytes, '
                f'above the {MAX_BGZF_BLOCK}-byte BGZF limit')
        if position > len(self.block):
            raise ValueError('Index points beyond BGZF block')
        self.offset, self.next, self.pos = offset, offset+size, position

    def tell(self):
        return (self.offset << 16) | self.pos if self.pos < len(self.block) else self.next << 16

    def _take(self, data):
        """Track cumulative decompressed reads so one query cannot run away."""
        self.consumed += len(data)
        if self.consumed > MAX_BGZF_QUERY_BYTES:
            raise ValueError(
                f'{self.name}: BGZF query exceeded the {MAX_BGZF_QUERY_BYTES}-byte read cap')
        return data

    def read(self, size):
        result = bytearray()
        while len(result) < size:
            if self.pos == len(self.block):
                self.seek(self.next << 16)
                if not self.block:
                    break
            take = min(size-len(result), len(self.block)-self.pos)
            result.extend(self.block[self.pos:self.pos+take])
            self.pos += take
        return self._take(bytes(result))

    def readline(self):
        result = bytearray()
        while True:
            if self.pos == len(self.block):
                self.seek(self.next << 16)
                if not self.block:
                    return self._take(bytes(result))
            end = self.block.find(b'\n', self.pos)
            if end >= 0:
                result.extend(self.block[self.pos:end+1])
                self.pos = end+1
                return self._take(bytes(result))
            result.extend(self.block[self.pos:])
            self.pos = len(self.block)


@lru_cache(maxsize=8)
def _bam_header(identity):
    f = BGZF(identity[0])
    try:
        if exact(f, 4) != b'BAM\x01':
            raise ValueError('Invalid BAM header')
        text = exact(f, struct.unpack('<i', exact(f, 4))[0]).decode('utf8', 'replace')
        refs = []
        for _ in range(struct.unpack('<i', exact(f, 4))[0]):
            name = exact(f, struct.unpack('<i', exact(f, 4))[0]).rstrip(b'\0').decode()
            refs.append((name, struct.unpack('<i', exact(f, 4))[0]))
        return refs, text
    finally:
        f.close()


@lru_cache(maxsize=8)
def _region_index(identity):
    path = identity[0]
    with open_binary(path) as probe:
        compressed=probe.read(2)==b'\x1f\x8b'
    from contextlib import ExitStack
    with ExitStack() as stack:
        raw=stack.enter_context(open_binary(path))
        f=stack.enter_context(gzip.GzipFile(fileobj=raw)) if compressed else raw
        magic = exact(f, 4)
        if magic not in {b'BAI\x01', b'TBI\x01', b'CSI\x01'}:
            raise ValueError('Unsupported genomic index')
        names = []
        min_shift, depth = 14, 5
        if magic == b'CSI\x01':
            min_shift, depth, auxlen = struct.unpack('<iii', exact(f, 12))
            aux = exact(f, auxlen)
            if len(aux) >= 28:
                names = aux[28:].rstrip(b'\0').decode().split('\0')
            nref = struct.unpack('<i', exact(f, 4))[0]
        else:
            nref = struct.unpack('<i', exact(f, 4))[0]
            if magic == b'TBI\x01':
                fields = struct.unpack('<7i', exact(f, 28))
                names = exact(f, fields[-1]).rstrip(b'\0').decode().split('\0')
        if not 0 <= depth <= 10 or not 0 <= min_shift <= 30:
            raise ValueError('Invalid genomic index geometry')
        if not 0 <= nref <= MAX_INDEX_REFS:
            raise ValueError(
                f'{path}: index declares {nref} references, above the {MAX_INDEX_REFS} limit')
        refs = []
        total_bins = total_chunks = 0
        for _ in range(nref):
            bins = {}
            for _ in range(struct.unpack('<i', exact(f, 4))[0]):
                key = struct.unpack('<I', exact(f, 4))[0]
                if magic == b'CSI\x01':
                    exact(f, 8)  # loffset; chunks are sufficient for correctness
                n = struct.unpack('<i', exact(f, 4))[0]
                total_bins += 1
                total_chunks += max(n, 0)
                if total_bins > MAX_INDEX_ENTRIES or total_chunks > MAX_INDEX_ENTRIES:
                    raise ValueError(
                        f'{path}: index declares more than {MAX_INDEX_ENTRIES} bins or chunks')
                bins[key] = [struct.unpack('<QQ', exact(f, 16)) for _ in range(n)]
            if magic != b'CSI\x01':
                n = struct.unpack('<i', exact(f, 4))[0]
                total_chunks += max(n, 0)
                if total_chunks > MAX_INDEX_ENTRIES:
                    raise ValueError(
                        f'{path}: index declares more than {MAX_INDEX_ENTRIES} bins or chunks')
                exact(f, n*8)
            refs.append(bins)
        return refs, names, min_shift, depth


def chunks(index, ref, start, end):
    refs, _, shift, depth = _region_index(version(index))
    if not 0 <= ref < len(refs):
        return []
    selected = []
    for level in range(depth+1):
        bits = shift+3*(depth-level)
        offset = ((1 << (3*level))-1)//7
        lo, hi = offset+(start >> bits), offset+((end-1) >> bits)
        # Iterate stored bins, never a potentially enormous empty range.
        for key, values in refs[ref].items():
            if lo <= key <= hi:
                selected.extend(values)
    result: list[tuple] = []
    for a, b in sorted(selected):
        if result and a <= result[-1][1]:
            result[-1] = (result[-1][0], max(b, result[-1][1]))
        else:
            result.append((a, b))
    return result


def index_path(path):
    from urllib.error import HTTPError
    suffixes=('.bai','.csi') if source_path(path).lower().endswith('.bam') else ('.tbi','.csi')
    candidates=[sidecar(path,suffix) for suffix in suffixes]
    if source_path(path).lower().endswith('.bam'):
        candidates.insert(1,sidecar(path,'.bai',True))
    for candidate in candidates:
        try:
            if exists(candidate):
                return candidate
        except HTTPError as exc:
            if exc.code != 403:
                raise
    return None


def cigar_parts(text):
    if text=='*':
        return []
    import re
    parts = [(op, int(n)) for n, op in re.findall(r'(\d+)([MIDNSHP=XB])', text)]
    if ''.join(f'{n}{op}' for op, n in parts) != text:
        raise ValueError(f'Invalid CIGAR: {text}')
    return parts


def read_geometry(start, cigar):
    ref, query, blocks, edits = start, 0, [], []
    for op, size in cigar:
        if op in 'M=X':
            blocks.append([ref, ref+size, query])
            ref += size
            query += size
        elif op in 'DN':
            edits.append([op, ref, size, query])
            ref += size
        elif op in 'IS':
            edits.append([op, ref, size, query])
            query += size
        elif op not in 'HP':
            raise ValueError(f'Unsupported CIGAR operation {op}')
    return ref, blocks, edits


def bam_tags(data):
    tags: dict[str, object] = {}
    at=0
    formats={'c':('b',1),'C':('B',1),'s':('h',2),'S':('H',2),'i':('i',4),'I':('I',4),'f':('f',4),'d':('d',8)}
    while at<len(data):
        if at+3>len(data):
            raise ValueError('Truncated BAM auxiliary tag')
        name=data[at:at+2].decode('ascii')
        kind=chr(data[at+2])
        at+=3
        if kind=='A':
            value: object = chr(data[at])
            at+=1
        elif kind in formats:
            fmt,size=formats[kind]
            value=struct.unpack_from('<'+fmt,data,at)[0]
            at+=size
        elif kind in 'ZH':
            end=data.find(b'\0',at)
            if end<0:
                raise ValueError('Unterminated BAM auxiliary string')
            value=data[at:end].decode('utf8','replace')
            at=end+1
        elif kind=='B':
            subtype=chr(data[at])
            fmt,size=formats[subtype]
            count=struct.unpack_from('<I',data,at+1)[0]
            at+=5
            if count>len(data)//size:
                raise ValueError('Invalid BAM auxiliary array length')
            value=list(struct.unpack_from('<'+str(count)+fmt,data,at))
            at+=count*size
        else:
            raise ValueError('Unsupported BAM auxiliary type '+kind)
        tags[name]=value
    return tags


def bam_record(data, refs):
    rid, start, binmqnl, flagnc, nseq, mate, mate_start, tlen = struct.unpack_from('<iiIIiiii', data)
    nname, ncigar, flag = binmqnl & 255, flagnc & 65535, flagnc >> 16
    at = 32
    name = data[at:at+nname].rstrip(b'\0').decode('utf8', 'replace')
    at += nname
    cigar = []
    for value in struct.unpack_from('<'+'I'*ncigar, data, at):
        if value & 15 >= len(CIGAR):
            raise ValueError('Invalid BAM CIGAR operation')
        cigar.append((CIGAR[value & 15], value >> 4))
    at += 4*ncigar
    packed = data[at:at+(nseq+1)//2]
    seq = ''.join(BASES[v >> 4]+BASES[v & 15] for v in packed)[:nseq]
    at += (nseq+1)//2
    qualities = list(data[at:at+nseq])
    tags=bam_tags(data[at+nseq:])
    if 'CG' in tags:
        cigar=[(CIGAR[v & 15],v >> 4) for v in tags.pop('CG')]
    end, blocks, edits = read_geometry(start, cigar)
    return dict(name=name, chrom=refs[rid][0] if rid >= 0 else '*', start=start, end=end,
                flag=flag, mapq=(binmqnl >> 8) & 255, cigar=''.join(f'{n}{op}' for op,n in cigar),
                sequence=seq, qualities=qualities, blocks=blocks, edits=edits,
                mate_chrom=refs[mate][0] if mate >= 0 else '*', mate_start=mate_start,
                template_length=tlen, tags=tags)


def bam_reads(path, region, options):
    refs, _ = _bam_header(version(path))
    chrom = chromosome(dict(refs), region.chrom)
    rid = [n for n, _ in refs].index(chrom)
    index = options.get('index') or index_path(path)
    if not index or not exists(index):
        raise ValueError('BAM needs a BAI/CSI index. Use Build index or provide a sidecar index.')
    f = BGZF(path)
    try:
        for a,b in chunks(index, rid, region.start, region.end):
            f.seek(a)
            while f.tell() < b:
                raw = f.read(4)
                if not raw:
                    break
                data = exact(f, struct.unpack('<i', raw)[0])
                record_ref, start = struct.unpack_from('<ii', data)
                if record_ref != rid or start >= region.end:
                    continue
                try:
                    record = bam_record(data, refs)
                except (IndexError, KeyError, struct.error, UnicodeDecodeError) as exc:
                    raise ValueError(
                        f'{path}: malformed BAM record at virtual offset {f.tell()} ({exc})'
                    ) from exc
                if record['end'] > region.start:
                    yield record
    finally:
        f.close()


def tool_output(tool, args):
    program = shutil.which(tool)
    if not program:
        raise ValueError(f'{tool} is required for this format. Install it and add it to PATH.')
    # Redirect stderr to a temporary file so a full error pipe cannot deadlock.
    import tempfile
    with tempfile.TemporaryFile() as errors:
        p = subprocess.Popen([program, *args], stdout=subprocess.PIPE, stderr=errors, text=True)
        lines = p.stdout
        assert lines is not None  # stdout=PIPE always yields a stream
        try:
            yield from lines
            if p.wait() != 0:
                errors.seek(0)
                raise ValueError(errors.read(16000).decode('utf8', 'replace'))
        finally:
            lines.close()
            if p.poll() is None:
                p.terminate()
                try:
                    p.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    p.kill()
                    p.wait()


def sam_record(line):
    cols = line.rstrip().split('\t')
    if len(cols) < 11:
        raise ValueError('Truncated SAM alignment')
    start = int(cols[3])-1
    end, blocks, edits = read_geometry(start, cigar_parts(cols[5]))
    return dict(name=cols[0], chrom=cols[2], start=start, end=end, flag=int(cols[1]),
                mapq=int(cols[4]), cigar=cols[5], sequence=cols[9],
                qualities=[ord(c)-33 for c in cols[10]] if cols[10] != '*' else [], blocks=blocks, edits=edits,
                mate_chrom=cols[2] if cols[6]=='=' else cols[6], mate_start=int(cols[7])-1,
                template_length=int(cols[8]), tags=cols[11:])


def cram_reads(path, region, options):
    refs = inspect_sequence(path)['chromsizes']
    chrom = chromosome(dict(refs), region.chrom)
    args = ['view']
    if options.get('reference'):
        args += ['-T', options['reference']]
    args += [str(path), f'{chrom}:{region.start+1}-{region.end}']
    for line in tool_output('samtools', args):
        yield sam_record(line)


def fasta_identity(path):
    """Cache identity for FASTA content: the FASTA version plus the .fai's own
    version, because fixing the index alone must invalidate cached bases."""
    index = sidecar(path,'.fai') if is_remote(path) else str(path)+'.fai'
    try:
        if exists(index):
            return version(path), version(index)
    except (OSError, ValueError):
        pass
    return version(path), None


@lru_cache(maxsize=8)
def fasta_index(identity, index_identity=None):
    if is_remote(identity[0]):
        records=[]
        for line in text_lines(sidecar(identity[0],'.fai')):
            cols=line.rstrip().split('\t')
            records.append((cols[0],*map(int,cols[1:5])))
        if not records:
            raise ValueError('Remote FASTA requires a nonempty .fai index')
        return records
    p = Path(identity[0])
    if Path(str(p)+'.fai').is_file():
        records = []
        for line in Path(str(p)+'.fai').read_text().splitlines():
            c = line.split('\t')
            records.append((c[0], *map(int,c[1:5])))
        return records
    if p.suffix.lower() == '.gz':
        raise ValueError('Compressed FASTA needs BGZF .fai/.gzi indexes and samtools')
    result = []
    with p.open('rb') as f:
        name, length, offset, bases, width, previous = None, 0, 0, 0, 0, 0
        while line := f.readline():
            if line.startswith(b'>'):
                if name is not None:
                    result.append((name, length, offset, bases, width))
                name = line[1:].split()[0].decode()
                length, offset, bases, width, previous = 0, f.tell(), 0, 0, 0
            elif name is not None:
                n = len(line.rstrip(b'\r\n'))
                if bases and previous != bases:
                    raise ValueError('Uneven FASTA lines require a valid .fai index')
                if not bases:
                    bases, width = n, len(line)
                length += n
                previous = n
        if name is not None:
            result.append((name, length, offset, bases, width))
    if not result or any(n <= 0 or b <= 0 for _,n,_,b,_ in result):
        raise ValueError('FASTA contains no sequence')
    return result


def reference_sequence(path, region):
    records = fasta_index(*fasta_identity(path))
    chrom = chromosome({r[0] for r in records}, region.chrom)
    _, length, offset, bases, width = next(r for r in records if r[0]==chrom)
    if region.end > length:
        raise ValueError('Region extends beyond reference sequence')
    if source_path(path).lower().endswith('.gz'):
        if not exists(sidecar(path,'.gzi')):
            raise ValueError('Compressed FASTA needs an existing .gzi index')
        lines = tool_output('samtools', ['faidx', str(path), f'{chrom}:{region.start+1}-{region.end}'])
        return ''.join(line.strip() for line in lines if not line.startswith('>')).upper()
    result, pos = bytearray(), region.start
    with open_binary(path) as f:
        while pos < region.end:
            count = min(region.end-pos, bases-pos % bases)
            f.seek(offset+(pos//bases)*width+pos % bases)
            result.extend(exact(f,count))
            pos += count
    return result.decode('ascii').upper()


@lru_cache(maxsize=8)
def bigwig_chromsizes(identity):
    with open_binary(identity[0]) as f:
        head = exact(f,64)
        endian = '<' if struct.unpack_from('<I', head)[0] == 0x888FFC26 else '>'
        if struct.unpack_from(endian+'I', head)[0] != 0x888FFC26:
            raise ValueError('Invalid bigWig magic')
        root = struct.unpack_from(endian+'Q',head,8)[0]
        f.seek(root)
        magic, _, key_size, value_size, count, _ = struct.unpack(endian+'IIIIQQ',exact(f,32))
        if magic != 0x78CA8C91 or value_size != 8 or key_size > 65536:
            raise ValueError('Invalid bigWig chromosome tree')
        pending, visited, result = [root+32], set(), []
        while pending:
            node = pending.pop()
            if node in visited:
                raise ValueError('Cyclic bigWig chromosome tree')
            visited.add(node)
            f.seek(node)
            leaf, _, n = struct.unpack(endian+'BBH',exact(f,4))
            children = []
            for _ in range(n):
                key = exact(f,key_size).rstrip(b'\0').decode()
                if leaf:
                    _, size = struct.unpack(endian+'II',exact(f,8))
                    result.append((key,size))
                else:
                    children.append(struct.unpack(endian+'Q',exact(f,8))[0])
            pending.extend(reversed(children))
        if len(result) != count:
            raise ValueError('Incomplete bigWig chromosome tree')
        return result


def file_kind(path):
    name = source_path(path).lower()
    for suffix, kind in [('.bw','signal'),('.bigwig','signal'),('.bedgraph','signal'),('.bg','signal'),
                         ('.bam','alignment'),('.cram','alignment'),('.sam','alignment'),
                         ('.fa','reference'),('.fasta','reference'),('.fna','reference'),
                         ('.vcf','variants'),('.bcf','variants'),('.gtf','genes'),('.gff','genes'),
                         ('.gff3','genes'),('.bed','intervals')]:
        if name.removesuffix('.gz').endswith(suffix):
            return kind
    raise ValueError('Unsupported sequence/track format')


def text_lines(path):
    import io
    from contextlib import ExitStack
    with ExitStack() as stack:
        raw=stack.enter_context(open_binary(path))
        binary=stack.enter_context(gzip.GzipFile(fileobj=raw)) if source_path(path).lower().endswith('.gz') else raw
        text=stack.enter_context(io.TextIOWrapper(binary,encoding='utf8'))
        yield from text


def vcf_span(cols):
    """Zero-based half-open span of a VCF data row. POS is 1-based; symbolic
    records must honor INFO/END so navigation matches interval_records()."""
    start = int(cols[1])-1
    end = start+max(1,len(cols[3]))
    for field in cols[7].split(';'):
        if field.startswith('END='):
            try:
                stop = int(field[4:])
            except ValueError:
                continue
            if stop > start:
                end = max(end, stop)
    return start, end


def inspect_sequence(path, chromsizes=None, index=None):
    path = str(path) if is_remote(path) else str(Path(path).resolve())
    kind = file_kind(path)
    refs = []
    if kind == 'reference':
        refs = [(r[0],r[1]) for r in fasta_index(*fasta_identity(path))]
    elif source_path(path).lower().endswith(('.bw','.bigwig')):
        refs = bigwig_chromsizes(version(path))
    elif source_path(path).lower().endswith('.bam'):
        refs,_ = _bam_header(version(path))
    elif source_path(path).lower().endswith(('.cram','.sam')):
        lines=tool_output('samtools',['view','-H',path]) if source_path(path).lower().endswith('.cram') else text_lines(path)
        for line in lines:
            if not line.startswith('@'):
                lines.close()
                break
            if line.startswith('@SQ\t'):
                fields = dict(part.split(':',1) for part in line.rstrip().split('\t')[1:] if ':' in part)
                refs.append((fields['SN'],int(fields['LN'])))
    elif chromsizes:
        refs=list(chromsizes)
    else:
        import re
        sizes = {}
        declared = set()
        lines = tool_output('bcftools',['view','-h',path]) if source_path(path).lower().endswith('.bcf') else text_lines(path)
        try:
            # Header metadata pass: declared contig lengths end the scan at
            # #CHROM so a huge indexed body is never read line by line.
            first_data = None
            for line in lines:
                if line.startswith('##FASTA'):
                    break
                if line.startswith(('track ','browser ')):
                    continue
                if line.startswith('##contig='):
                    name = re.search(r'ID=([^,>]+)',line)
                    length = re.search(r'length=(\d+)',line,re.I)
                    if name and length:
                        sizes[name[1]] = int(length[1])
                        declared.add(name[1])
                elif line.startswith('#CHROM'):
                    if sizes:
                        break
                elif line.startswith('#') or not line.strip():
                    continue
                else:
                    first_data = line  # headerless or missing lengths: infer below
                    break
            if first_data is not None or not sizes:
                from itertools import chain
                rows = chain(() if first_data is None else (first_data,), lines)
                for i,line in enumerate(rows):
                    c = line.rstrip('\r\n').split('\t')
                    end = vcf_span(c)[1] if kind=='variants' else int(c[4] if kind=='genes' else c[2])
                    # Declared header lengths stay authoritative over inference.
                    if c[0] not in declared:
                        sizes[c[0]] = max(sizes.get(c[0],0),end)
                    if i >= 1_000_000:
                        raise ValueError('Text track metadata exceeds one million rows; load a reference genome first and add this as a track')
            refs = list(sizes.items())
        finally:
            lines.close()
    if not refs:
        raise ValueError('No chromosome lengths or intervals found in this file')
    return dict(path=path, kind=kind, name=Path(source_path(path)).name, chromsizes=refs,
                lengths_exact=kind in {'reference','alignment'} or source_path(path).lower().endswith(('.bw','.bigwig')),
                index=index or (index_path(path) if kind in {'alignment','variants'} or source_path(path).lower().endswith('.gz') and kind!='reference' else None), visible=True, color='#397ca8', scale='auto', mode='fill')


def tabix_lines(path, region, index=None):
    if not source_path(path).lower().endswith('.gz'):
        yield from text_lines(path)
        return
    index = index or index_path(path)
    if not index or not source_path(path).lower().endswith('.gz'):
        yield from text_lines(path)
        return
    _, names, _, _ = _region_index(version(index))
    chrom = chromosome(names,region.chrom)
    f = BGZF(path)
    try:
        for a,b in chunks(index,names.index(chrom),region.start,region.end):
            f.seek(a)
            while f.tell() < b:
                line = f.readline()
                if not line:
                    break
                yield line.decode('utf8')
    finally:
        f.close()


@lru_cache(maxsize=8)
def variant_samples(identity):
    path=identity[0]
    lines=tool_output('bcftools',['view','-h',path]) if source_path(path).lower().endswith('.bcf') else text_lines(path)
    try:
        for line in lines:
            if line.startswith('##FASTA'):
                break
            if line.startswith(('track ','browser ')):
                continue
            if line.startswith('#CHROM'):
                return line.rstrip().split('\t')[9:]
            if not line.startswith('#'):
                break
        return []
    finally:
        lines.close()


def interval_records(track, region, limit=MAX_RECORDS):
    path, kind = track['path'], track['kind']
    if source_path(path).lower().endswith('.bcf'):
        chrom = chromosome(dict(inspect_sequence(path)['chromsizes']), region.chrom)
        lines = tool_output('bcftools',['view','-H','-r',f'{chrom}:{region.start+1}-{region.end}',path])
    else:
        lines = tabix_lines(path,region,track.get('index'))
    samples, result, scanned = variant_samples(version(path)) if kind=='variants' else [], [], 0
    try:
        for line in lines:
            if line.startswith('##FASTA'):
                break
            if line.startswith(('track ','browser ')):
                continue
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
            if line.startswith('#') or not line.strip():
                continue
            scanned += 1
            if scanned > 2_000_000:
                raise ValueError('Unindexed text query exceeds two million records; bgzip/tabix the file for regional access')
            c = line.rstrip('\r\n').split('\t')
            if c[0].removeprefix('chr') != region.chrom.removeprefix('chr'):
                continue
            if kind == 'variants':
                info = dict(v.split('=',1) if '=' in v else (v,True) for v in c[7].split(';') if v!='.')
                start, end = vcf_span(c)
                if end <= region.start or start >= region.end:
                    continue
                if track.get('pass_only') and c[6] not in {'PASS','.'}:
                    continue
                fmt = c[8].split(':') if len(c)>8 else []
                record = dict(start=start,end=end,name=c[2],ref=c[3],alt=c[4],qual=c[5],filter=c[6],info=info,
                              samples={samples[i] if i<len(samples) else str(i+1):dict(zip(fmt,v.split(':'),strict=False)) for i,v in enumerate(c[9:])})
            else:
                start, end = (int(c[3])-1,int(c[4])) if kind=='genes' else (int(c[1]),int(c[2]))
                if end <= region.start or start >= region.end:
                    continue
                record = dict(start=start,end=end,name=c[8] if kind=='genes' else (c[3] if len(c)>3 else ''),
                              strand=c[6] if kind=='genes' else (c[5] if len(c)>5 else '.'),raw=line.rstrip())
                if kind=='intervals' and len(c)>=12:
                    sizes=list(map(int,c[10].rstrip(',').split(',')))
                    starts=list(map(int,c[11].rstrip(',').split(',')))
                    if len(sizes)!=int(c[9]) or len(starts)!=int(c[9]):
                        raise ValueError('BED12 block count does not match block coordinates')
                    record['blocks']=[[start+a,start+a+n] for a,n in zip(starts,sizes,strict=False)]
            result.append(record)
            if len(result) > limit:
                return dict(records=result[:limit],truncated=True,scanned=scanned)
            if scanned >= 2_000_000:
                raise ValueError('Unindexed text query exceeds two million records; bgzip/tabix the file for regional access')
        return dict(records=result,truncated=False,scanned=scanned)
    finally:
        lines.close()


def alignment_query(track, region, bins, reference=''):
    options = {**track, 'reference':reference}
    if source_path(track['path']).lower().endswith('.bam'):
        reads = bam_reads(track['path'],region,options)
    elif source_path(track['path']).lower().endswith('.cram'):
        reads = cram_reads(track['path'],region,options)
    else:
        reads = (sam_record(raw) for raw in text_lines(track['path']) if not raw.startswith('@'))
    coverage = np.zeros(bins,dtype=float)
    delta = np.zeros(bins+1,dtype=float)
    import random
    random_source=random.Random(0)
    record_limit=max(0,int(track.get('record_limit',0)))
    junctions: dict = {}
    selected: list = []
    seen, accepted, truncated = 0,0,False
    excluded = 4 | (0 if track.get('duplicates') else 1024) | (0 if track.get('secondary') else 256) | (0 if track.get('supplementary',True) else 2048) | 512
    max_reads = min(10000,max(1,int(track.get('max_reads',1000))))
    span = region.end-region.start
    try:
        for read in reads:
            seen += 1
            if seen%100000==0:
                print('@poreckit-progress:'+json.dumps(dict(message=f'Reading alignments · {seen:,} records · {accepted:,} accepted')),flush=True)
            if record_limit and seen > record_limit:
                truncated=True
                break
            if read['chrom'].removeprefix('chr') != region.chrom.removeprefix('chr') or read['end'] <= region.start or read['start'] >= region.end:
                continue
            if read['flag'] & excluded or read['mapq'] < int(track.get('mapq',0)):
                continue
            accepted += 1
            for a,b,_ in read['blocks']:
                a,b = max(a,region.start),min(b,region.end)
                if b<=a:
                    continue
                lo,hi = (a-region.start)*bins/span, (b-region.start)*bins/span
                first,last = max(0,int(lo)), min(bins-1,int(math.ceil(hi)-1))
                if first==last:
                    coverage[first] += hi-lo
                else:
                    coverage[first] += first+1-lo
                    coverage[last] += hi-last
                    delta[first+1] += 1
                    delta[last] -= 1
            for op,pos,size,_ in read['edits']:
                if op=='N':
                    key=(pos,pos+size)
                    junctions[key]=junctions.get(key,0)+1
            if len(selected)<max_reads:
                selected.append(read)
            else:
                slot=random_source.randrange(accepted)
                if slot<max_reads:
                    selected[slot]=read
        sort = track.get('sort','position')
        selected.sort(key={'name':lambda r:r['name'],'strand':lambda r:(bool(r['flag'] & 16),r['start']),
                           'mapq':lambda r:(-r['mapq'],r['start'])}.get(sort,lambda r:r['start']))
        coverage+=np.cumsum(delta[:-1])
        return dict(records=selected,coverage=coverage.tolist(),accepted=accepted,scanned=seen,
                    junctions=[dict(start=a,end=b,count=count) for (a,b),count in sorted(junctions.items(),key=lambda value:-value[1])[:2000]],
                    truncated=truncated,reads_hidden=max(0,accepted-len(selected)))
    finally:
        reads.close()


def query_sequence(job, state=None):
    from .model import Region
    from .tracks import query_signal
    region = Region.parse(job['region'])
    bins = min(4000,max(1,min(region.end-region.start,int(job.get('bins',1000)))))
    reference = job.get('reference','')
    sequence,reference_error='', ''
    if reference and region.end-region.start<=20000:
        try:
            sequence=reference_sequence(reference,region)
        except (ValueError,OSError) as exc:
            reference_error=str(exc)
    result = dict(region=region.text(),sequence=sequence,tracks=[],reference=reference,reference_error=reference_error)
    shared: dict = {}
    for track in job['tracks']:
        if not track.get('visible',True) or track['kind']=='reference':
            continue
        try:
            if reference and track.get('lengths_exact') and track.get('chromsizes'):
                reference_lengths={v[0]:v[1] for v in fasta_index(*fasta_identity(reference))}
                source_lengths=dict(track['chromsizes'])
                ref_name=chromosome(reference_lengths,region.chrom)
                source_name=chromosome(source_lengths,region.chrom)
                if reference_lengths[ref_name]!=source_lengths[source_name]:
                    raise ValueError('Reference and track chromosome lengths differ; select matching genome assemblies')
            if track['kind']=='signal':
                values = np.asarray(query_signal(track['path'],region.chrom,region.start,region.end,bins,track.get('index')),dtype=float)
                finite = values[np.isfinite(values)]
                lo,hi = (min(0.,float(finite.min())),max(0.,float(finite.max()))) if finite.size else (0.,1.)
                if track.get('scale')=='fixed':
                    lo,hi = float(track.get('ymin',0)),float(track.get('ymax',1))
                    if not lo<hi:
                        raise ValueError('Signal minimum must be below maximum')
                group = track.get('group','') if track.get('scale')=='shared' else ''
                if group:
                    a,b=shared.get(group,(0.,0.))
                    shared[group]=(min(a,lo),max(b,hi))
                data=dict(values=[float(v) if np.isfinite(v) else None for v in values],ymin=lo,ymax=hi)
            elif track['kind']=='alignment':
                if state is None:
                    data=alignment_query(track,region,bins,reference)
                else:
                    defaults=dict(mapq=0,duplicates=False,secondary=False,supplementary=True,
                                  max_reads=1000,record_limit=0,sort='position')
                    key=('alignments',version(track['path']),region,bins,
                         version(reference) if reference else None,
                         version(track['index']) if track.get('index') else None,
                         tuple((k,track.get(k,default)) for k,default in defaults.items()))
                    data=state.sequence_payload(key,lambda track=track:alignment_query(track,region,bins,reference))
            else:
                data=interval_records(track,region,limit=10000)
                if track['kind']=='genes':
                    data['records']=gene_models(data['records'])
            result['tracks'].append(dict(track=track,data=data))
        except (ValueError,OSError,struct.error) as exc:
            result['tracks'].append(dict(track=track,error=str(exc)))
    for entry in result['tracks']:
        t=entry['track']
        if 'data' in entry and t.get('scale')=='shared' and t.get('group') in shared:
            entry['data']['ymin'],entry['data']['ymax']=shared[t['group']]
    return result


def gene_models(records):
    """Assemble exon/CDS rows into transcript models for GTF and GFF3."""
    from .genesearch import attributes
    models,genes={},{}
    for record in records:
        cols=record['raw'].split('\t')
        attrs=attributes(cols[8])
        feature=cols[2].lower()
        if feature=='gene':
            key=attrs.get('gene_id') or attrs.get('ID') or attrs.get('Name')
            genes[key]={**record,'name':attrs.get('gene_name') or attrs.get('Name') or key}
        if feature not in {'transcript','mrna','exon','cds','ncrna'}:
            continue
        ids=([attrs['transcript_id']] if attrs.get('transcript_id') else
             attrs.get('Parent','').split(',') if feature in {'exon','cds'} else
             [attrs.get('ID') or attrs.get('Name') or record['name']])
        for key in ids:
            if not key:
                continue
            if key not in models:
                models[key]={**record,'name':attrs.get('gene_name') or attrs.get('Name') or key,
                             'id':key,'gene':attrs.get('gene_id') or (attrs.get('Parent') if feature in {'transcript','mrna'} else ''),
                             'blocks':[],'cds':[]}
            model=models[key]
            model['start'],model['end']=min(model['start'],record['start']),max(model['end'],record['end'])
            if feature in {'exon','cds'}:
                model['blocks' if feature=='exon' else 'cds'].append([record['start'],record['end']])
    if not models:
        return list(genes.values()) or records
    represented=set()
    for model in models.values():
        represented.add(model['gene'])
        if model['gene'] in genes:
            model['name']=genes[model['gene']]['name']+' · '+model['id']
        model['blocks']=sorted(model['blocks'] or model['cds'] or [[model['start'],model['end']]])
    return [*models.values(),*[record for key,record in genes.items() if key not in represented and not any(record['start']<=m['start'] and record['end']>=m['end'] for m in models.values())]]
