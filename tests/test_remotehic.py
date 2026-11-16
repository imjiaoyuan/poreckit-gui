"""Remote .hic over HTTP byte ranges: parity with the local poreckit reader."""
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np
import poreckit
import pytest
from hicfixture import build, dense_block, dense_block_v8

from poreckit_gui.remotehic import RemoteHicReader


def write_v9_hic(path, **kwargs):
    """Minimal v9 .hic built from bytes (chr1/chr2, 1bp bins, VC vectors)."""
    path.write_bytes(build(**kwargs))
    return path


class _Handler(BaseHTTPRequestHandler):
    files = {}
    requests = []

    def log_message(self, *args):
        pass

    def _serve(self, include_body=True):
        path = self.path.split('?')[0]
        data = self.files.get(path)
        if data is None:
            self.send_error(404)
            return
        self.requests.append((path, self.command, dict(self.headers)))
        base = {'Content-Length': str(len(data)), 'ETag': '"hic-1"',
                'Accept-Ranges': 'bytes'}
        if self.command == 'HEAD':
            self.send_response(200)
            for key, value in base.items():
                self.send_header(key, value)
            self.end_headers()
            return
        range_header = self.headers.get('Range')
        if range_header:
            start, stop = (int(part) for part in
                           range_header.removeprefix('bytes=').split('-'))
            stop = min(stop, len(data) - 1)
            chunk = data[start:stop + 1]
            self.send_response(206)
            self.send_header('Content-Range', f'bytes {start}-{stop}/{len(data)}')
            self.send_header('Content-Length', str(len(chunk)))
            self.end_headers()
            if include_body:
                self.wfile.write(chunk)
            return
        self.send_response(200)
        for key, value in base.items():
            self.send_header(key, value)
        self.end_headers()
        if include_body:
            self.wfile.write(data)

    do_HEAD = _serve

    def do_GET(self):
        self._serve()


@pytest.fixture
def hic_server(tmp_path):
    write_v9_hic(tmp_path / 'map.hic')
    handler = type('Handler', (_Handler,), {'files': {}, 'requests': []})
    handler.files['/map.hic'] = (tmp_path / 'map.hic').read_bytes()
    server = ThreadingHTTPServer(('127.0.0.1', 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    handler.requests.clear()
    yield f'http://127.0.0.1:{server.server_address[1]}/map.hic', handler
    server.shutdown()
    thread.join()


def _local(url, tmp_path):
    name = url.rsplit('/', 1)[-1]
    return poreckit.make_reader(str(tmp_path / name))


@pytest.mark.parametrize('version', [9, 8])
def test_remote_hic_matches_local_reader(tmp_path, version):
    dense = dense_block if version > 8 else dense_block_v8
    body = build(version=version, w1=(2.0, 2.0), w2=(1.0, float('nan')),
                 blk_intra=dense([5, 0, 0, 5], 2), blk_inter=dense([3, 0, 4, 7], 2))
    local_path = tmp_path / f'map{version}.hic'
    local_path.write_bytes(body)
    handler = type('Handler', (_Handler,), {'files': {}, 'requests': []})
    handler.files['/map.hic'] = body
    server = ThreadingHTTPServer(('127.0.0.1', 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        url = f'http://127.0.0.1:{server.server_address[1]}/map.hic'
        remote = RemoteHicReader(url)
        local = poreckit.make_reader(str(local_path))
        assert remote.chromsizes() == local.chromsizes()
        assert remote.resolutions() == local.resolutions()
        assert remote.meta() == local.meta()
        np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1), local.matrix('chr1', 0, 2, 1))
        np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1, 'VC'),
                                   local.matrix('chr1', 0, 2, 1, 'VC'))
        np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1, 'raw'),
                                   local.matrix('chr1', 0, 2, 1, 'raw'))
        np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1, 'VC', chrom2='chr2'),
                                   local.matrix('chr1', 0, 2, 1, 'VC', chrom2='chr2'),
                                   equal_nan=True)
        np.testing.assert_allclose(remote.matrix('chr2', 0, 2, 1, chrom2='chr1'),
                                   local.matrix('chr2', 0, 2, 1, chrom2='chr1'))
        with pytest.raises(ValueError, match='not found'):
            remote.matrix('chr1', 0, 2, 1, 'KR')
        with pytest.raises(ValueError, match='not in .hic'):
            remote.matrix('chrX', 0, 2, 1)
        with pytest.raises(ValueError, match='same-chromosome'):
            remote.matrix('chr1', 0, 2, 1, chrom2='chr1')
        assert isinstance(remote.matrix('chr1', 0, 2, 1), np.ndarray)
    finally:
        server.shutdown()


def test_remote_hic_sparse_block_and_subwindows(tmp_path):
    """v9 type-1 sparse block (upper triangle + mirroring) and diagonal
    sub-window queries must match the local reader exactly."""
    import struct as _struct
    import threading
    import zlib

    from hicfixture import build
    body = bytearray()
    body += _struct.pack('<i', 4)  # nRecords
    body += _struct.pack('<i', 0) + _struct.pack('<i', 0)  # bin offsets
    body += bytes([1, 1, 0, 1])  # useShort f32 counts, X wide i32, Y u16, sparse
    body += _struct.pack('<H', 2) + _struct.pack('<H', 0)  # row dy=0
    body += _struct.pack('<I', 2) + _struct.pack('<I', 0) + _struct.pack('<f', 3.0)
    body += _struct.pack('<I', 1) + _struct.pack('<f', 4.0)
    body += _struct.pack('<H', 1)  # row dy=1
    body += _struct.pack('<I', 1) + _struct.pack('<I', 1) + _struct.pack('<f', 7.0)
    data = build(blk_intra=zlib.compress(bytes(body)))
    path = tmp_path / 'sparse.hic'
    path.write_bytes(data)
    handler = type('Handler', (_Handler,), {'files': {'/sparse.hic': data}, 'requests': []})
    server = ThreadingHTTPServer(('127.0.0.1', 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        url = f'http://127.0.0.1:{server.server_address[1]}/sparse.hic'
        remote = RemoteHicReader(url)
        local = poreckit.make_reader(str(path))
        for region in ((0, 2), (0, 1), (1, 2)):
            np.testing.assert_allclose(remote.matrix('chr1', *region, 1),
                                       local.matrix('chr1', *region, 1))
            np.testing.assert_allclose(remote.matrix('chr1', *region, 1, 'VC'),
                                       local.matrix('chr1', *region, 1, 'VC'))
    finally:
        server.shutdown()


def test_remote_hic_range_requests_and_bounded_caching(hic_server, tmp_path):
    url, handler = hic_server
    remote = RemoteHicReader(url)
    local = _local(url, tmp_path)
    assert remote.chromsizes() == local.chromsizes() == [('chr1', 2), ('chr2', 2)]
    assert remote.resolutions() == [1]
    np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1), [[5, 0], [0, 5]])
    # every ranged response must be validated against the server ETag
    assert handler.requests
    assert all(call[2].get('If-Match') == '"hic-1"' for call in handler.requests
               if call[1] == 'GET')
    assert all(call[2].get('Accept-Encoding') == 'identity'
               for call in handler.requests)
    # served ranges never exceed the bounded 256 KiB block size
    for call in handler.requests:
        if call[1] == 'GET' and 'Range' in call[2]:
            start, stop = (int(p) for p in
                           call[2]['Range'].removeprefix('bytes=').split('-'))
            assert stop - start < 256 * 1024
    count = len(handler.requests)
    remote.matrix('chr1', 0, 2, 1)  # cached blocks and index: no new requests
    assert len(handler.requests) == count
    # inter-chromosomal query
    np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1, chrom2='chr2'),
                               [[3, 4], [0, 7]])
    # a chrom pair without a master-index entry answers with zeros
    np.testing.assert_allclose(remote.matrix('chr1', 0, 2, 1, chrom2='chr2', norm='NONE'),
                               [[3, 4], [0, 7]])


def test_remote_hic_worker_inspect(tmp_path):
    """worker.inspect_file dispatches remote .hic URLs to the range reader."""
    from poreckit_gui import worker
    path = write_v9_hic(tmp_path / 'inspect.hic')
    handler = type('Handler', (_Handler,), {'files': {}, 'requests': []})
    handler.files['/inspect.hic'] = path.read_bytes()
    server = ThreadingHTTPServer(('127.0.0.1', 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        url = f'http://127.0.0.1:{server.server_address[1]}/inspect.hic'
        meta = worker.inspect_file(url)
        assert meta['resolutions'] == [1]
        assert meta['chromsizes'] == [('chr1', 2), ('chr2', 2)]
        assert meta['assembly'] == 'hg19'
        local_meta = worker.inspect_file(str(path))
        meta.pop('input'), local_meta.pop('input')
        assert meta == local_meta
    finally:
        server.shutdown()


def test_remote_hic_open_routing(tmp_path):
    """File → Open routes a .hic URL into the contact-map inspect job and
    still rejects remote track formats on the main map."""
    from PySide6.QtWidgets import QApplication

    from poreckit_gui.app import Browser
    _app = QApplication.instance() or QApplication([])
    browser = Browser()
    jobs, errors = [], []
    browser.start_job = lambda job, message: jobs.append(job)
    browser.error = errors.append
    try:
        browser.load_file('https://data.test/map.hic?token=signed')
        assert jobs and jobs[0]['action'] == 'inspect'
        assert jobs[0]['input'] == 'https://data.test/map.hic?token=signed'
        browser.load_file('https://data.test/data.xyz')
        assert errors and 'remote' in errors[-1]
        assert len(jobs) == 1
    finally:
        browser.close()



def test_render_cache_remote_matches_local_render(hic_server, tmp_path):
    """RenderCache.tiled_rectangle over an HTTP .hic URL must produce the
    same matrix as the local-file render of the same bytes."""
    from poreckit_gui.cache import RenderCache
    from poreckit_gui.model import Region
    url, _handler = hic_server
    local_path = tmp_path / 'map.hic'
    x = Region('chr1', 0, 2)
    y = Region('chr2', 0, 2)
    remote_cache, local_cache = RenderCache(max_readers=4), RenderCache(max_readers=4)
    try:
        remote = remote_cache.tiled_rectangle(url, x, y, 1)
        local = local_cache.tiled_rectangle(str(local_path), x, y, 1)
        np.testing.assert_allclose(remote, local)
        remote = remote_cache.tiled_rectangle(url, x, x, 1)
        local = local_cache.tiled_rectangle(str(local_path), x, x, 1)
        np.testing.assert_allclose(remote, local)
        # matrix cache is bounded and versioned on the URL identity
        assert remote_cache.version(url)[0] == url
    finally:
        remote_cache.close()
        local_cache.close()


def test_render_cache_remote_keys_differ_and_invalidate(hic_server, tmp_path):
    """Two URLs yield distinct cache keys; a changed remote identity
    (server-side change) must invalidate the cached version key."""
    import poreckit_gui.remote as remote_module
    from poreckit_gui.cache import RenderCache
    url, handler = hic_server
    cache = RenderCache()
    try:
        base = cache.version(url)
        # a second URL (even for identical bytes) is a different cache entry
        other = f'{url.rsplit("/", 1)[0]}/other.hic'
        handler.files['/other.hic'] = handler.files['/map.hic']
        assert cache.version(other) != base
        # embedding identity(): a server-side change (new ETag) invalidates
        original = remote_module.identity
        remote_module.identity = lambda p: original(p)[:-1] + (999,) if p == url else original(p)
        try:
            assert cache.version(url) != base
        finally:
            remote_module.identity = original
    finally:
        cache.close()
