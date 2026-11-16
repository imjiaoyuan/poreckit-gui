import io
import re
import struct
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlsplit

import numpy as np
import pytest

from poreckit_gui import remote
from poreckit_gui.model import Region
from poreckit_gui.sequenceio import inspect_sequence, alignment_query, reference_sequence, query_sequence
from test_sequence import write_bam


class Response(io.BytesIO):
    def __init__(self, data, headers, status=200):
        super().__init__(data)
        self.headers, self.status = headers, status


@pytest.fixture
def server(monkeypatch):
    files, calls = {}, []
    def request(url, headers=None, method=None):
        key = urlsplit(url).path
        if key not in files:
            raise HTTPError(url,404,'Missing',{},None)
        data = files[key]
        calls.append((key,method,dict(headers or {})))
        meta = {'Content-Length': str(len(data)), 'ETag': '"version-1"'}
        if method == 'HEAD':
            return Response(b'',meta)
        if headers and 'Range' in headers:
            a,b = map(int,re.fullmatch(r'bytes=(\d+)-(\d+)', headers['Range']).groups())
            return Response(data[a:b+1],{**meta,'Content-Length':str(b-a+1),'Content-Range':f'bytes {a}-{b}/{len(data)}'},206)
        return Response(data,meta)
    monkeypatch.setattr(remote,'_request',request)
    remote.clear_cache()
    yield files,calls
    remote.clear_cache()


def test_remote_ranges_seek_cache_and_signed_sidecars(server):
    files,calls=server
    files['/data.bam']=bytes(range(256))*5000
    url='https://data.test/data.bam?token=secret'
    with remote.RangeFile(url) as f:
        f.seek(remote.BLOCK_SIZE-5)
        assert f.read(20)==files['/data.bam'][remote.BLOCK_SIZE-5:remote.BLOCK_SIZE+15]
        f.seek(-12,2)
        assert f.read()==files['/data.bam'][-12:]
    count=len(calls)
    with remote.RangeFile(url) as f:
        f.seek(remote.BLOCK_SIZE-5)
        assert f.read(20)==files['/data.bam'][remote.BLOCK_SIZE-5:remote.BLOCK_SIZE+15]
    assert len(calls)==count
    assert all(c[2].get('If-Match')=='"version-1"' for c in calls if c[1]!='HEAD')
    assert remote.sidecar(url,'.bai')=='https://data.test/data.bam.bai?token=secret'
    assert remote.sidecar('/tmp/sample#1.bam','.bai')=='/tmp/sample#1.bam.bai'


def test_indexed_remote_bam_and_fasta_match_local_queries(server,tmp_path):
    files,_=server
    bam=write_bam(tmp_path)
    files['/reads.bam']=bam.read_bytes()
    files['/reads.bam.bai']=Path(str(bam)+'.bai').read_bytes()
    url='https://data.test/reads.bam?token=signed'
    track=inspect_sequence(url)
    assert track['index']=='https://data.test/reads.bam.bai?token=signed'
    result=alignment_query(track,Region('1',100,120),20)
    local=alignment_query(inspect_sequence(bam),Region('1',100,120),20)
    assert result==local
    files['/ref.fa']=b'>chr1\nACGTACGT\nAACCGGTT\n'
    files['/ref.fa.fai']=b'chr1\t16\t6\t8\t9\n'
    ref='https://data.test/ref.fa'
    assert inspect_sequence(ref)['chromsizes']==[('chr1',16)]
    assert reference_sequence(ref,Region('1',5,13))=='CGTAACCG'


def bigwig_fixture():
    tree=struct.pack('<IIIIQQ',0x78CA8C91,1,4,8,1,0)+struct.pack('<BBH4sII',1,0,1,b'chr1',0,1000)
    offset=64+len(tree)
    data=struct.pack('<IIIIIBBH',0,100,300,0,0,1,0,2)+struct.pack('<IIfIIf',100,200,7.,200,300,-3.)
    index=offset+len(data)
    spatial=struct.pack('<IIQIIIIQII',0x2468ACE0,1,1,0,100,0,300,index,1,0)
    spatial+=struct.pack('<BBHIIIIQQ',1,0,1,0,100,0,300,offset,len(data))
    head=struct.pack('<IHHQQQHHQQIQ',0x888FFC26,4,0,64,offset,index,0,0,0,0,0,0)
    return head+tree+data+spatial


def test_remote_bigwig_preserves_missing_negative_and_partial_bins(server,tmp_path):
    from poreckit_gui.tracks import query_signal
    files,_=server
    raw=bigwig_fixture()
    files['/signal.bw']=raw
    local=tmp_path/'signal.bw'
    local.write_bytes(raw)
    actual=query_signal('https://data.test/signal.bw','1',50,350,7)
    expected=query_signal(str(local),'chr1',50,350,7)
    np.testing.assert_allclose(actual,[np.nan,7,7,2,-3,-3,np.nan],equal_nan=True)
    np.testing.assert_allclose(actual,expected,equal_nan=True,rtol=1e-6)
    assert actual[-2]<0 and np.isnan(actual[0])
    data=query_sequence(dict(region='chr1:50-350',tracks=[inspect_sequence('https://data.test/signal.bw')],bins=7))
    assert 'error' not in data['tracks'][0]
    assert data['tracks'][0]['data']['ymin']<0


def test_rejects_ignored_ranges_and_wrong_content_range(monkeypatch):
    def no_range(url, headers=None, method=None):
        return Response(b'',{'Content-Length':str(remote.BLOCK_SIZE*4)},200)
    monkeypatch.setattr(remote,'_request',no_range)
    remote.clear_cache()
    with remote.RangeFile('https://data.test/large.bam') as f:
        with pytest.raises(ValueError,match='Range'):
            f.read(1)
    def wrong(url,headers=None,method=None):
        return Response(b'abc',{'Content-Length':'3','Content-Range':'bytes 1-3/3'},200 if method=='HEAD' else 206)
    monkeypatch.setattr(remote,'_request',wrong)
    remote.clear_cache()
    with remote.RangeFile('https://data.test/a') as f:
        with pytest.raises(ValueError,match='wrong HTTP'):
            f.read(1)
    remote.clear_cache()


def test_remote_url_session_roundtrip_preserves_sources(server,tmp_path):
    from PySide6.QtWidgets import QApplication
    from poreckit_gui.app import Browser
    app=QApplication.instance() or QApplication([])
    files,_=server
    bam=write_bam(tmp_path)
    files['/reads.bam']=bam.read_bytes()
    files['/custom-index']=Path(str(bam)+'.bai').read_bytes()
    url='https://data.test/reads.bam?token=signed'
    index='https://data.test/custom-index?token=signed'
    b=Browser()
    seq=b.sequence_browser
    def request(job,callback):
        if job['action']=='sequence_inspect':
            callback(dict(tracks=[inspect_sequence(p,job.get('chromsizes'),job.get('indexes',{}).get(p)) for p in job['paths']]))
        else:
            callback(query_sequence(job))
    seq.request=request
    try:
        seq.add_files([url],indexes={url:index})
        assert not seq.result['tracks'][0].get('error')
        assert seq.result['tracks'][0]['data']['accepted']==2
        state=seq.session_data(tmp_path)
        assert state['tracks'][0]['path']==url
        assert state['tracks'][0]['index']==index
        seq.restore(state,tmp_path)
        assert seq.result['tracks'][0]['track']['path']==url
        assert seq.result['tracks'][0]['data']['accepted']==2
    finally:
        b.close()


def test_remote_vcf_tbi_preserves_sample_names_and_interval_edges(server,tmp_path):
    import gzip
    from test_sequence import bgzf
    from poreckit_gui.sequenceio import interval_records
    files,_=server
    header=b'##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n'
    row=b'chr1\t11\tv\tA\tG\t50\tPASS\t.\tGT\t0/1\n'
    files['/calls.vcf.gz']=bgzf(header+row)+bgzf(b'')
    names=b'chr1\0'
    index=b'TBI\x01'+struct.pack('<8i',1,2,1,2,0,35,0,len(names))+names
    index+=struct.pack('<iIiQQi',1,0,1,len(header),len(header)+len(row),0)
    files['/calls.vcf.gz.tbi']=gzip.compress(index)
    url='https://data.test/calls.vcf.gz?token=signed'
    track=inspect_sequence(url)
    records=interval_records(track,Region('chr1',10,11))['records']
    assert records[0]['samples']=={'Sample':{'GT':'0/1'}}
    assert interval_records(track,Region('chr1',11,12))['records']==[]


def test_local_bedgraph_matches_bigwig_weighted_bins(tmp_path):
    from poreckit_gui.tracks import query_signal
    path=tmp_path/'signal.bg'
    path.write_text('track type=bedGraph\nchr1\t100\t200\t7\nchr1\t200\t300\t-3\n')
    actual=query_signal(str(path),'1',50,350,7)
    np.testing.assert_allclose(actual,[np.nan,7,7,2,-3,-3,np.nan],equal_nan=True)


def test_same_size_metadata_without_validators_eventually_invalidates(monkeypatch):
    """Regression: remote sources without ETag/Last-Modified but the same size
    kept metadata caches valid forever; identity now carries a generation."""
    class Response:
        def __init__(self, headers): self.headers = headers
        def __enter__(self): return self
        def __exit__(self, *args): return False
    def request(url, headers=None, method=None):
        return Response({'Content-Length': '5'})
    monkeypatch.setattr(remote, '_request', request)
    remote.clear_cache()
    url = 'https://data.test/stream.bin'
    monkeypatch.setattr(remote, '_generation', lambda: 1)
    first = remote.identity(url)
    monkeypatch.setattr(remote, '_generation', lambda: 2)
    second = remote.identity(url)
    assert first[:4] == second[:4] == (url, 5, '', '')
    assert first != second and second[-1] == 2
    monkeypatch.undo()
    remote.clear_cache()
