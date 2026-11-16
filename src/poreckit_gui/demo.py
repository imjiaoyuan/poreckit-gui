"""Small, explicitly synthetic data for trying the browser without downloads."""
from pathlib import Path


def create_demo(folder, chromosomes=None, scale=1):
    import h5py
    import numpy as np

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "synthetic-demo.cool"
    step = 10000
    chromosomes = chromosomes or [("chr1", 2000000)]
    lengths = [length for name, length in chromosomes]
    sizes = [(length+step-1)//step for length in lengths]
    n = sum(sizes)
    i, j = np.triu_indices(n)
    rng = np.random.default_rng(7)
    count = rng.poisson(70 / (1 + abs(i - j) / 3)
                        + 18 * ((i // 40) == (j // 40))).astype("int32")
    count *= scale
    with h5py.File(path, "w") as f:
        f.attrs.update({"format": "HDF5::Cooler", "format-version": 3,
                        "bin-size": step, "bin-type": "fixed", "storage-mode": "symmetric-upper",
                        "assembly": "synthetic-demo"})
        chroms = f.create_group("chroms")
        chroms["name"] = np.array([name.encode() for name, length in chromosomes])
        chroms["length"] = np.array(lengths, dtype="int64")
        bins = f.create_group("bins")
        bins["chrom"] = np.repeat(np.arange(len(sizes), dtype="int32"), sizes)
        bins["start"] = np.concatenate([np.arange(size, dtype="int64") * step for size in sizes])
        bins["end"] = np.concatenate([np.minimum((np.arange(size, dtype="int64")+1)*step, length) for size, length in zip(sizes, lengths)])
        pixels = f.create_group("pixels")
        pixels["bin1_id"], pixels["bin2_id"], pixels["count"] = i, j, count
        indexes = f.create_group("indexes")
        indexes["chrom_offset"] = np.r_[0, np.cumsum(sizes)].astype("int64")
        indexes["bin1_offset"] = np.searchsorted(i, np.arange(n + 1)).astype("int64")
    return str(path)


def create_track_demo(directory):
    """Synthetic reference, signal, reads and variants; no biological meaning."""
    root=Path(directory)
    root.mkdir(parents=True,exist_ok=True)
    sequence='ACGT'*500
    reference=root/'synthetic.fa'
    reference.write_text('>chr1\n'+'\n'.join(sequence[i:i+80] for i in range(0,len(sequence),80))+'\n')
    signal=root/'synthetic.bedgraph'
    signal.write_text('chr1\t0\t100\t2\nchr1\t100\t160\t8\nchr1\t160\t250\t-3\nchr1\t250\t2000\t1\n')
    reads=root/'synthetic.sam'
    a=sequence[100:140]
    b=sequence[110:150]
    b=b[:15]+('A' if b[15]!='A' else 'T')+b[16:]
    inserted=sequence[120:130]+'AA'+sequence[130:150]
    reads.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n'+
                     f'pair1\t99\tchr1\t101\t60\t40M\t=\t171\t110\t{a}\t'+('I'*40)+'\tNM:i:0\n'+
                     f'read2\t16\tchr1\t111\t45\t40M\t*\t0\t0\t{b}\t'+('I'*40)+'\tNM:i:1\n'+
                     f'insertion\t0\tchr1\t121\t60\t10M2I20M\t*\t0\t0\t{inserted}\t'+('I'*32)+'\tNM:i:2\n'+
                     f'pair1\t147\tchr1\t171\t60\t40M\t=\t101\t-110\t{sequence[170:210]}\t'+('I'*40)+'\tNM:i:0\n')
    variants=root/'synthetic.vcf'
    variants.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=2000>\n'+
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleA\tSampleB\n'+
                        'chr1\t126\tsynthetic_snv\tC\tA\t60\tPASS\tDP=20\tGT:DP\t0/1:20\t0/0:18\n'+
                        'chr1\t130\tsynthetic_insertion\tC\tCAA\t50\tPASS\tDP=12\tGT:DP\t0/1:12\t1/1:10\n')
    genes=root/'synthetic.gtf'
    genes.write_text('chr1\tsynthetic\tgene\t81\t230\t.\t+\t.\tgene_id "DEMO1"; gene_name "SyntheticGene";\n'+
                     'chr1\tsynthetic\texon\t81\t140\t.\t+\t.\tgene_id "DEMO1"; transcript_id "DEMO1.1";\n'+
                     'chr1\tsynthetic\texon\t171\t230\t.\t+\t.\tgene_id "DEMO1"; transcript_id "DEMO1.1";\n')
    return [str(p) for p in (reference,signal,reads,variants,genes)]
