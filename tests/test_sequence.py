import gzip
import json
import shutil
import struct
import time
import zlib
from pathlib import Path

import numpy as np
import pytest

samtools_required = pytest.mark.skipif(shutil.which('samtools') is None,
                                       reason='samtools is not installed on this runner')

from poreckit_gui.model import Region
from poreckit_gui.sequenceio import (
    BGZF,alignment_query,bam_reads,bigwig_chromsizes,inspect_sequence,
    interval_records,query_sequence,reference_sequence,version,
)


def bgzf(data):
    compressor=zlib.compressobj(wbits=-15)
    compressed=compressor.compress(data)+compressor.flush()
    size=18+len(compressed)+8
    header=b'\x1f\x8b\x08\x04'+b'\0'*4+b'\0\xff'+struct.pack('<H',6)+b'BC'+struct.pack('<HH',2,size-1)
    return header+compressed+struct.pack('<II',zlib.crc32(data),len(data))


def write_bam(tmp_path):
    """Small independent BAM fixture including clipping, insertion and intron."""
    path=tmp_path/'reads.bam'
    text=b'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n'
    header=b'BAM\x01'+struct.pack('<i',len(text))+text+struct.pack('<ii',1,5)+b'chr1\0'+struct.pack('<i',2000)
    records=[]
    for name,start,flag,mapq,cigar,sequence in [
        ('read1',100,0,60,[(0,4),(1,2),(0,3),(2,1),(0,3),(3,5),(0,2)],'ACGTAACGTACGTA'),
        ('read2',102,16,20,[(4,2),(0,8)],'TTAAAAAAAA'),
        ('duplicate',105,1024,60,[(0,5)],'AAAAA'),
        ('secondary',108,256,60,[(0,4)],'CCCC'),
    ]:
        name=name.encode()+b'\0'
        code='=ACMGRSVTWYHKDBN'
        encoded=bytes((code.index(sequence[i])<<4)|(code.index(sequence[i+1]) if i+1<len(sequence) else 0) for i in range(0,len(sequence),2))
        payload=struct.pack('<iiIIiiii',0,start,(4681<<16)|(mapq<<8)|len(name),(flag<<16)|len(cigar),len(sequence),-1,-1,0)
        payload+=name+b''.join(struct.pack('<I',(n<<4)|op) for op,n in cigar)+encoded+bytes([35]*len(sequence))+b'NMi'+struct.pack('<i',1)
        records.append(struct.pack('<i',len(payload))+payload)
    data=header+b''.join(records)
    path.write_bytes(bgzf(data)+bgzf(b''))
    # Root bin covers all records. Chunk starts at a virtual offset inside the block.
    index=b'BAI\x01'+struct.pack('<iiIiQQi',1,1,0,1,len(header),len(data),0)
    Path(str(path)+'.bai').write_bytes(index)
    return path


def wait(app,predicate,timeout=15):
    deadline=time.monotonic()+timeout
    while not predicate():
        app.processEvents()
        time.sleep(.01)
        assert time.monotonic()<deadline,'GUI task timed out'


def test_bgzf_seeks_and_reads_across_blocks(tmp_path):
    p=tmp_path/'blocks'
    first=bgzf(b'abc\ndef')
    p.write_bytes(first+bgzf(b'ghi\njkl')+bgzf(b''))
    f=BGZF(p)
    try:
        assert f.readline()==b'abc\n'
        assert f.readline()==b'defghi\n'
        f.seek((len(first)<<16)|4)
        assert f.read(3)==b'jkl'
    finally:
        f.close()


def test_indexed_bam_filter_cigar_and_coverage(tmp_path):
    p=write_bam(tmp_path)
    meta=inspect_sequence(p)
    assert meta['chromsizes']==[('chr1',2000)]
    r=Region('1',100,120)
    data=alignment_query(meta,r,100)
    assert data['accepted']==2
    assert data['records'][0]['cigar']=='4M2I3M1D3M5N2M'
    assert data['records'][0]['sequence']=='ACGTAACGTACGTA'
    assert data['records'][0]['end']==118
    assert data['records'][0]['edits']==[['I',104,2,4],['D',107,1,9],['N',111,5,12]]
    # Mean depth: reference aligned bases 12+8 over a 20 bp window.
    assert np.mean(data['coverage'])==pytest.approx(1.)
    assert len(list(bam_reads(p,Region('chr1',117,118),{})))==1
    data=alignment_query({**meta,'mapq':30,'duplicates':True},r,100)
    assert [v['name'] for v in data['records']]==['read1','duplicate']
    Path(str(p)+'.bai').unlink()
    with pytest.raises(ValueError,match='BAI/CSI'):
        alignment_query(meta,r,100)


def test_fasta_regions_vcf_genotypes_and_signed_shared_signals(tmp_path):
    fasta=tmp_path/'ref.fa'
    fasta.write_text('>chr1\nACGTAC\nGTACGT\nAC\n>chr2\nTTTT\n')
    assert reference_sequence(fasta,Region('1',4,13))=='ACGTACGTA'
    assert inspect_sequence(fasta)['chromsizes']==[('chr1',14),('chr2',4)]
    vcf=tmp_path/'variants.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=14>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleA\n'
                   'chr1\t3\trs1\tG\tA\t60\tPASS\tDP=20\tGT:DP\t0/1:20\n'
                   'chr1\t8\tindel\tTAC\tT\t10\tLowQual\t.\tGT:DP\t1/1:4\n')
    variants=interval_records(dict(path=str(vcf),kind='variants'),Region('chr1',0,14))
    assert variants['records'][0]['start']==2
    assert variants['records'][0]['samples']['SampleA']=={'GT':'0/1','DP':'20'}
    assert len(interval_records(dict(path=str(vcf),kind='variants',pass_only=True),Region('chr1',0,14))['records'])==1
    tracks=[]
    for i,value in enumerate([-3,7]):
        p=tmp_path/f'{i}.bg'
        p.write_text(f'chr1\t0\t14\t{value}\n')
        tracks.append(dict(path=str(p),kind='signal',name=str(i),scale='shared',group='RNA'))
    result=query_sequence(dict(region='chr1:0-14',reference=str(fasta),tracks=tracks))
    assert result['sequence']=='ACGTACGTACGTAC'
    assert all(t['data']['ymin']==-3 and t['data']['ymax']==7 for t in result['tracks'])


def test_tabix_header_preserves_sample_names(tmp_path):
    p=tmp_path/'variants.vcf.gz'
    header=b'##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPatient\n'
    row=b'chr1\t11\tv\tA\tG\t50\tPASS\t.\tGT\t0/1\n'
    p.write_bytes(bgzf(header+row)+bgzf(b''))
    names=b'chr1\0'
    index=b'TBI\x01'+struct.pack('<8i',1,2,1,2,0,35,0,len(names))+names
    index+=struct.pack('<iIiQQi',1,0,1,len(header),len(header)+len(row),0)
    Path(str(p)+'.tbi').write_bytes(gzip.compress(index))
    track=dict(path=str(p),kind='variants')
    result=interval_records(track,Region('1',10,11))
    assert result['records'][0]['samples']=={'Patient':{'GT':'0/1'}}


def test_independent_browser_opens_bam_tracks_exports_and_restores(tmp_path,monkeypatch):
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    from poreckit_gui.dataview import export_scene
    app=QApplication.instance() or QApplication([])
    b=Browser()
    p=write_bam(tmp_path)
    signal=tmp_path/'sample.bg'
    signal.write_text('chr1\t0\t2000\t3\n')
    ref=tmp_path/'ref.fa'
    ref.write_text('>chr1\n'+'A'*2000+'\n')
    seq=b.sequence_browser
    try:
        b.load_file(str(p))
        wait(app,lambda:not seq.busy() and seq.result is not None)
        assert not b.input_path
        assert b.tabs.currentWidget() is seq
        assert seq.result['tracks'][0]['data']['accepted']==2
        seq.add_files([signal,ref])
        wait(app,lambda:not seq.busy() and len(seq.result['tracks'])==2)
        seq.navigate(Region('chr1',95,130))
        wait(app,lambda:not seq.busy())
        assert len(seq.result['sequence'])==35
        for suffix in ['svg','png','pdf']:
            dest=tmp_path/f'preview.{suffix}'
            export_scene(seq.canvas.scene(),str(dest))
            assert dest.stat().st_size>100
        data=seq.session_data(tmp_path)
        seq.restore(json.loads(json.dumps(data)),tmp_path)
        wait(app,lambda:not seq.busy())
        assert seq.region==Region('chr1',95,130)
        assert seq.tracks.count()==2
        assert seq.reference==str(ref)
        seq.canvas.selected.emit(seq.result['tracks'][0]['data']['records'][0])
        assert 'read1' in seq.detail.toPlainText()
    finally:
        b.close()


def test_track_demo_clone_lives_after_original_closes_and_links_loci(tmp_path):
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    app=QApplication.instance() or QApplication([])
    original=Browser()
    clone=None
    try:
        original.open_track_demo()
        seq=original.sequence_browser
        wait(app,lambda:not seq.busy() and seq.result is not None)
        assert not any('error' in entry for entry in seq.result['tracks'])
        models=next(t['data']['records'] for t in seq.result['tracks'] if t['track']['kind']=='genes')
        assert len(models)==1 and len(models[0]['blocks'])==2
        clone=original.windows.create(duplicate=True)
        other=clone.sequence_browser
        wait(app,lambda:not other.busy() and other.result is not None)
        assert Path(other.reference).is_relative_to(other.temp.name)
        seq.track_link.setChecked(True)
        other.track_link.setChecked(True)
        seq.navigate(Region('chr1',90,220))
        wait(app,lambda:not other.busy() and not seq.busy())
        assert other.region==seq.region==Region('chr1',90,220)
        original.close()
        original=None
        other.zoom(.5)
        wait(app,lambda:not other.busy())
        assert not any('error' in entry for entry in other.result['tracks'])
        assert 'ACGT' in other.result['sequence']
        with pytest.raises(ValueError,match='temporary'):
            other.session_data(tmp_path)
    finally:
        if original:
            original.close()
        if clone:
            clone.close()


def test_reference_mismatch_reports_track_error_and_queries_survive_bad_tracks(tmp_path):
    bam=write_bam(tmp_path)
    ref=tmp_path/'short.fa'
    ref.write_text('>chr1\n'+'A'*1000+'\n')
    signal=tmp_path/'signal.bg'
    signal.write_text('chr1\t0\t1000\t-2\n')
    result=query_sequence(dict(region='chr1:100-120',reference=str(ref),tracks=[inspect_sequence(bam),inspect_sequence(signal)]))
    assert 'lengths differ' in result['tracks'][0]['error']
    assert result['tracks'][1]['data']['ymin']==-2
    assert len(result['tracks'][1]['data']['values'])==20


def test_high_depth_coverage_is_complete_and_read_display_is_sampled(tmp_path,monkeypatch):
    import poreckit_gui.sequenceio as module
    count=100_100
    def reads(*args):
        for i in range(count):
            yield dict(name=str(i),chrom='chr1',start=0,end=10,flag=0,mapq=60,
                       blocks=[[0,10,0]],edits=[['N',10,5,10]])
    monkeypatch.setattr(module,'bam_reads',reads)
    track=dict(path='virtual.bam',max_reads=20)
    full=alignment_query(track,Region('chr1',0,20),20)
    assert full['accepted']==count and not full['truncated']
    assert full['coverage'][:10]==[count]*10
    assert full['coverage'][10:]==[0]*10
    assert len(full['records'])==20
    assert max(int(r['name']) for r in full['records'])>count*.8
    assert full['junctions']==[dict(start=10,end=15,count=count)]
    limited=alignment_query({**track,'record_limit':100},Region('chr1',0,20),20)
    assert limited['truncated'] and limited['accepted']==100


def test_display_changes_reuse_query_but_filters_and_file_edits_invalidate(tmp_path,monkeypatch):
    import poreckit_gui.sequenceio as module
    from poreckit_gui.cache import RenderCache
    p=write_bam(tmp_path)
    track=inspect_sequence(p)
    original=module.alignment_query
    calls=[]
    def query(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(module,'alignment_query',query)
    state=RenderCache()
    try:
        job=dict(region='chr1:100-120',tracks=[track],bins=20)
        first=query_sequence(job,state)
        track.update(color='#ffffff',pairs=True,junction_arcs=True,baseq=30)
        second=query_sequence(job,state)
        assert len(calls)==1
        assert first['tracks'][0]['data']==second['tracks'][0]['data']
        # Returned records can be modified without poisoning cached payloads.
        second['tracks'][0]['data']['records'].clear()
        assert query_sequence(job,state)['tracks'][0]['data']['records']
        track['mapq']=30
        assert query_sequence(job,state)['tracks'][0]['data']['accepted']==1
        assert len(calls)==2
        import os
        stamp=p.stat().st_mtime_ns
        os.utime(p,ns=(stamp,stamp+1000000))
        query_sequence(job,state)
        assert len(calls)==3
        assert state.payload_bytes<=state.max_bytes//4
    finally:
        state.close()
    assert not state.payloads and state.payload_bytes==0


def test_vcf_length_inference_uses_half_open_coordinates_and_end(tmp_path):
    """Regression: inferred VCF lengths used POS+len(REF) (off by one) and
    ignored INFO/END for symbolic SVs."""
    vcf=tmp_path/'inferred.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   'chr1\t10\tsnp\tA\tG\t50\tPASS\t.\n'
                   'chr2\t10\tdel\tN\t<DEL>\t50\tPASS\tEND=1000\n'
                   'chr3\t10\tins\tACG\tT\t50\tPASS\t.\n')
    sizes=dict(inspect_sequence(vcf)['chromsizes'])
    assert sizes['chr1']==10      # POS-1+len(REF), not POS+len(REF)
    assert sizes['chr2']==1000    # symbolic SV honors INFO/END
    assert sizes['chr3']==12
    # Declared header lengths stay authoritative over inferred values.
    declared=tmp_path/'declared.vcf'
    declared.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=14>\n'
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                        'chr1\t3\tlong\tACGTACGTAC\tA\t50\tPASS\t.\n')
    assert inspect_sequence(declared)['chromsizes']==[('chr1',14)]


def test_vcf_metadata_scan_stops_at_chrom_when_lengths_declared(tmp_path,monkeypatch):
    """Regression: a VCF whose header declares contig lengths must not scan
    the whole body during metadata inspection."""
    import poreckit_gui.sequenceio as module
    header=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000000>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    def lines(path):
        yield from header.splitlines(keepends=True)
        yield 'chr1\t1\tv\tA\tG\t50\tPASS\t.\n'
        raise AssertionError('metadata scan read past #CHROM into the body')
    monkeypatch.setattr(module,'text_lines',lines)
    meta=module.inspect_sequence(tmp_path/'big.vcf')
    assert meta['chromsizes']==[('chr1',10000000)]


def test_fixing_fai_invalidates_cached_reference_bases(tmp_path):
    """Regression: fasta_index() cached on FASTA identity only, so a repaired
    .fai kept serving stale bases from the old index."""
    import os
    fasta=tmp_path/'ref.fa'
    # Byte-exact writes: text mode would translate \n to \r\n on Windows and
    # break the FAI line geometry the offsets refer to.
    fasta.write_bytes(b'>chr1\nACGTACGT\n')
    fai=Path(str(fasta)+'.fai')
    fai.write_bytes(b'chr1\t8\t7\t8\t9\n')  # wrong offset reads CGTA
    assert reference_sequence(fasta,Region('chr1',0,4))=='CGTA'
    fai.write_bytes(b'chr1\t8\t6\t8\t9\n')  # corrected offset reads ACGT
    stamp=fai.stat().st_mtime_ns
    os.utime(fai,ns=(stamp,stamp+1_000_000_000))
    assert reference_sequence(fasta,Region('chr1',0,4))=='ACGT'
    assert inspect_sequence(fasta)['chromsizes']==[('chr1',8)]


def test_signal_track_uses_custom_tabix_index(tmp_path,monkeypatch):
    """Regression: query_signal() dropped track['index'], so a custom TBI/CSI
    never reached tabix_lines()."""
    import poreckit_gui.sequenceio as module
    from poreckit_gui.tracks import query_signal
    p=tmp_path/'signal.bg.gz'
    row=b'chr1\t0\t14\t2\n'
    p.write_bytes(bgzf(row)+bgzf(b''))
    names=b'chr1\0'
    index=b'TBI\x01'+struct.pack('<8i',1,2,1,2,0,35,0,len(names))+names
    index+=struct.pack('<iIiQQi',1,0,1,0,len(row),0)
    custom=str(p)+'.custom.tbi'
    Path(custom).write_bytes(gzip.compress(index))
    assert not Path(str(p)+'.tbi').exists()
    values=query_signal(str(p),'chr1',0,14,14,custom)
    assert float(values[0])==2.0 and float(np.nanmean(values))==2.0
    seen=[]
    original=module.tabix_lines
    def spy(path,region,index=None):
        seen.append(index)
        yield from original(path,region,index)
    monkeypatch.setattr(module,'tabix_lines',spy)
    track=dict(path=str(p),kind='signal',name='s',index=custom)
    result=query_sequence(dict(region='chr1:0-14',tracks=[track]))
    assert seen and seen[0]==custom
    assert result['tracks'][0]['data']['values'][0]==2.0


def test_cancelling_worker_terminates_external_tools(tmp_path):
    """Regression: force-killing the worker once left its external tool
    children (started like samtools) running. The worker now runs in its own
    process group, receives SIGTERM on cancel, and its finally blocks run."""
    import os
    import sys
    import textwrap
    from PySide6.QtCore import QProcess
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.jobs import JobProcess
    if os.name != 'posix':
        pytest.skip('process-group cancellation requires POSIX')
    app=QApplication.instance() or QApplication([])
    # The grandchild sits in its own session; only the worker's finally block
    # can stop it — so the marker proves the graceful path executed.
    script=tmp_path/'worker_tool.py'
    script.write_text(textwrap.dedent('''
        import os, signal, subprocess, sys, time
        os.setsid()
        child=subprocess.Popen(["sleep","30"],start_new_session=True)
        open(sys.argv[1],"w").write(f"{child.pid} {os.getpgrp()}")
        signal.signal(signal.SIGTERM,lambda signum,frame:(_ for _ in ()).throw(SystemExit(0)))
        try:
            time.sleep(600)
        finally:
            child.terminate()
            child.wait()
            open(sys.argv[1],"a").write(" cleaned")
    '''))
    pidfile=tmp_path/'tool.pid'
    job,result=tmp_path/'job.json',tmp_path/'result.json'
    job.write_text('{}'); result.write_text('')
    process=JobProcess()
    done=[]
    process.finished.connect(lambda code,status:done.append(code))
    try:
        process.start(sys.executable,[str(script),str(pidfile),str(job),str(result)])
        deadline=time.monotonic()+15
        while not pidfile.exists() and time.monotonic()<deadline:
            app.processEvents(); time.sleep(.01)
        tool_pid,worker_group=pidfile.read_text().split()
        worker_pid=process.processId()
        assert int(worker_group)!=os.getpgrp(),'worker did not get its own process group'
        process.kill()
        deadline=time.monotonic()+15
        while not done and time.monotonic()<deadline:
            app.processEvents(); time.sleep(.01)
        assert done
        assert process.child.state()==QProcess.ProcessState.NotRunning
        deadline=time.monotonic()+15
        while 'cleaned' not in pidfile.read_text() and time.monotonic()<deadline:
            app.processEvents(); time.sleep(.01)
        assert 'cleaned' in pidfile.read_text(),'worker finally block never ran'
        deadline=time.monotonic()+15
        while tool_pid and Path(f'/proc/{tool_pid}').exists() and time.monotonic()<deadline:
            app.processEvents(); time.sleep(.01)
        assert not Path(f'/proc/{tool_pid}').exists(),'external tool survived cancel'
    finally:
        process.shutdown()


@samtools_required
def test_generated_bam_index_is_rebuilt_after_session_restore(tmp_path):
    """Regression: saving dropped the transient generated index and reopening
    failed with 'BAM needs a BAI/CSI index'."""
    import json as jsonlib
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    app=QApplication.instance() or QApplication([])
    b=Browser()
    seq=b.sequence_browser
    p=write_bam(tmp_path)
    Path(str(p)+'.bai').unlink()
    try:
        b.load_file(str(p))
        wait(app,lambda:not seq.busy() and seq.result is not None)
        seq.tracks.setCurrentRow(0)
        seq.build_index()
        wait(app,lambda:not seq.busy() and seq.result['tracks'][0].get('data') is not None)
        assert seq.result['tracks'][0]['data']['accepted']==2
        data=seq.session_data(tmp_path)
        track=data['tracks'][0]
        assert track.get('index_rebuild') is True and 'index' not in track
        seq.restore(jsonlib.loads(jsonlib.dumps(data)),tmp_path)
        wait(app,lambda:not seq.busy() and seq.result is not None
             and seq.result['tracks'][0].get('data',{}).get('accepted')==2)
        item=seq.tracks.item(0)
        assert item.data(_qt_user_role()).get('index')
    finally:
        b.close()


def _qt_user_role():
    from PySide6.QtCore import Qt
    return Qt.ItemDataRole.UserRole


def test_corrupt_or_truncated_bgzf_raises_error_naming_the_file(tmp_path):
    """Regression: corrupt/truncated BGZF blocks must raise a ValueError that
    identifies the file, not a raw zlib/gzip exception."""
    from poreckit_gui.sequenceio import BGZF
    data = bgzf(b'chr1\t0\t100\t1.5\n') + bgzf(b'chr1\t100\t200\t2.5\n')
    good = tmp_path / 'signal.bg.gz'
    good.write_bytes(data)
    f = BGZF(str(good))
    try:
        assert f.readline() == b'chr1\t0\t100\t1.5\n'
        assert f.readline() == b'chr1\t100\t200\t2.5\n'
        assert f.readline() == b''
    finally:
        f.close()
    truncated = tmp_path / 'cut.bg.gz'
    truncated.write_bytes(bgzf(b'chr1\t0\t100\t1.5\n')[:-10])
    with pytest.raises(ValueError) as exc:
        BGZF(str(truncated)).readline()
    assert truncated.name in str(exc.value)
    corrupt = tmp_path / 'bad.bg.gz'
    damaged = bytearray(data)
    damaged[30] ^= 0xFF  # flip a byte inside the first block's compressed data
    corrupt.write_bytes(bytes(damaged))
    with pytest.raises(ValueError) as exc:
        BGZF(str(corrupt)).readline()
    assert corrupt.name in str(exc.value)


def test_index_with_absurd_counts_is_rejected(tmp_path):
    """Regression: corrupt index headers implying huge allocations must be
    rejected with a clear ValueError instead of huge reads."""
    from poreckit_gui.sequenceio import chunks, version
    p = write_bam(tmp_path)
    absurd_refs = tmp_path / 'refs.bai'
    absurd_refs.write_bytes(b'BAI\x01' + struct.pack('<i', 20_000_000))
    with pytest.raises(ValueError) as exc:
        chunks(absurd_refs, 0, 0, 100)
    assert absurd_refs.name in str(exc.value) and 'reference' in str(exc.value)
    # One ref, one bin, absurd chunk count: nref, n_bins, bin_key, n_chunks.
    absurd_bins = tmp_path / 'bins.bai'
    absurd_bins.write_bytes(b'BAI\x01' + struct.pack('<iiii', 1, 1, 0, 11_000_000))
    with pytest.raises(ValueError) as exc:
        chunks(absurd_bins, 0, 0, 100)
    assert absurd_bins.name in str(exc.value) and 'bins or chunks' in str(exc.value)
    # A real index with modest counts still parses.
    assert chunks(Path(str(p)+'.bai'), 0, 0, 2000)


def write_bam_with_raw_tags(path, tags):
    """BAM fixture identical to write_bam's first record, but with caller-
    supplied auxiliary tag bytes so malformed tags can be exercised."""
    text = b'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n'
    header = (b'BAM\x01' + struct.pack('<i', len(text)) + text
              + struct.pack('<ii', 1, 5) + b'chr1\0' + struct.pack('<i', 2000))
    name = b'read1\0'
    cigar = [(0, 14)]
    code = '=ACMGRSVTWYHKDBN'
    sequence = 'ACGTAACGTACGTA'
    encoded = bytes((code.index(sequence[i]) << 4) | (code.index(sequence[i+1]) if i+1 < len(sequence) else 0)
                    for i in range(0, len(sequence), 2))
    payload = struct.pack('<iiIIiiii', 0, 100, (4681 << 16) | (60 << 8) | len(name),
                          (0 << 16) | len(cigar), len(sequence), -1, -1, 0)
    payload += name + b''.join(struct.pack('<I', (n << 4) | op) for op, n in cigar) \
        + encoded + bytes([35]*len(sequence)) + tags
    data = header + struct.pack('<i', len(payload)) + payload
    path.write_bytes(bgzf(data) + bgzf(b''))
    index = b'BAI\x01' + struct.pack('<iiIiQQi', 1, 1, 0, 1, len(header), len(data), 0)
    Path(str(path)+'.bai').write_bytes(index)
    return path


def test_malformed_bam_tag_yields_per_track_error(tmp_path):
    """Regression: one BAM with a truncated auxiliary tag must produce a
    per-track error (never an uncaught IndexError) while other tracks in the
    same query still return data."""
    bad = write_bam_with_raw_tags(tmp_path/'bad.bam', b'XYA')  # 'A' tag without its value
    (tmp_path / 'good').mkdir()
    good = write_bam(tmp_path / 'good')
    bad_track = inspect_sequence(bad)
    good_track = inspect_sequence(good)
    result = query_sequence(dict(region='chr1:100-120',
                                 tracks=[bad_track, good_track], bins=20))
    assert 'error' in result['tracks'][0]
    assert bad.name in result['tracks'][0]['error']
    assert 'BAM record' in result['tracks'][0]['error']
    assert result['tracks'][1]['data']['accepted'] > 0
