"""File-versioned gene-name indexes for GTF, GFF and BED annotations."""
import functools
import gzip
import re
from pathlib import Path
from typing import Any
from urllib.parse import unquote


def attributes(text):
    if re.match(r"\s*[^\s;=]+=", text):
        return {key.strip(): unquote(value.strip()) for item in text.split(";")
                if "=" in item for key, value in [item.split("=", 1)]}
    return {match[0]: match[1] or match[2].strip() for match in
            re.findall(r'(?:^|;)\s*([^\s;]+)\s+(?:"([^"]*)"|([^;]*))', text)}


@functools.lru_cache(maxsize=2)
def _index(path, mtime, size, inode):
    genes: dict[tuple[str, str], dict[str, Any]] = {}
    fallback: dict[tuple[str, str], dict[str, Any]] = {}
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as stream:
        for line_no, line in enumerate(stream, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 3:
                cols = line.split()
            if len(cols) >= 9 and not cols[2].isdigit():
                attrs = attributes(cols[8])
                kind = cols[2].lower()
                is_gene = kind == "gene" or kind.endswith("_gene") or kind == "pseudogene"
                gtf_id = attrs.get("gene_id") or attrs.get("gene_name")
                key = gtf_id or attrs.get("ID") or attrs.get("Name")
                if not key:
                    continue
                # Explicit GFF genes and GTF gene IDs define gene bounds;
                # transcript-only GFF files remain searchable as a fallback.
                destination = genes if is_gene or gtf_id else fallback
                if destination is fallback and kind not in {"mrna", "transcript", "ncrna"}:
                    continue
                try:
                    start, end = int(cols[3])-1, int(cols[4])
                except ValueError:
                    raise ValueError(f"Invalid annotation coordinates at {Path(path).name}:{line_no}") from None
                name = attrs.get("gene_name") or attrs.get("Name") or key
                aliases = {key, name}
                for field in ("gene_id", "transcript_id", "ID", "Alias", "gene", "Dbxref"):
                    aliases.update(v for v in attrs.get(field, "").split(",") if v)
                aliases |= {v.split(":", 1)[1] for v in aliases if ":" in v}
            else:
                try:
                    start, end = int(cols[1]), int(cols[2])
                except (ValueError, IndexError):
                    raise ValueError(f"Expected GTF, GFF or BED coordinates at {Path(path).name}:{line_no}") from None
                name = cols[3] if len(cols) > 3 and cols[3] not in {"", "."} else f"{cols[0]}:{start}-{end}"
                key = f"{name}:{start}:{end}"
                aliases, destination = {name}, genes
            if start < 0 or end <= start:
                raise ValueError(f"Invalid annotation interval at {Path(path).name}:{line_no}")
            identity = cols[0], key
            if identity not in destination:
                destination[identity] = dict(chrom=cols[0], start=start, end=end, name=name, aliases=set())
            record = destination[identity]
            record["start"], record["end"] = min(record["start"], start), max(record["end"], end)
            if name != key:
                record["name"] = name
            record["aliases"].update(aliases)
    values = list(genes.values() if genes else fallback.values())
    return tuple({**v, "aliases": tuple(sorted(v["aliases"])),
                  "search": tuple(a.casefold() for a in v["aliases"])} for v in values)


def gene_index(path):
    path = Path(path).resolve()
    stat = path.stat()
    return _index(str(path), stat.st_mtime_ns, stat.st_size, stat.st_ino)


def search_genes(paths, query, chromsizes, limit=200):
    query = query.strip().casefold()
    if not query:
        raise ValueError("Enter a gene name, ID or alias")
    chromosomes = {name.removeprefix("chr"): (name, length) for name, length in chromsizes}
    results, seen = [], set()
    for path in dict.fromkeys(paths):
        for gene in gene_index(path):
            if not any(query in alias for alias in gene["search"]):
                continue
            canonical = chromosomes.get(gene["chrom"].removeprefix("chr"))
            if not canonical:
                continue
            chrom, length = canonical
            start, end = gene["start"], min(gene["end"], length)
            if end <= start:
                continue
            key = chrom, start, end, gene["name"]
            if key in seen:
                continue
            seen.add(key)
            rank = 0 if query in gene["search"] else 1 if any(a.startswith(query) for a in gene["search"]) else 2
            results.append(dict(chrom=chrom, start=start, end=end, name=gene["name"], aliases=list(gene["aliases"]),
                                source=str(path), rank=rank))
    results.sort(key=lambda value: (value["rank"], value["name"].casefold(), value["chrom"], value["start"]))
    return dict(results=results[:limit], total=len(results), query=query)
