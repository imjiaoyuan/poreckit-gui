"""Seekable HTTP Range sources with bounded, version-aware byte caching."""
from __future__ import annotations

import io
import re
import time
from collections import OrderedDict
from urllib.error import HTTPError
from urllib.parse import urlsplit, urlunsplit
from urllib.request import Request, urlopen

BLOCK_SIZE = 256 * 1024
CACHE_BYTES = 64 * 1024 * 1024
_cache = OrderedDict()
_metadata = {}
_cache_size = 0


def is_remote(path):
    return urlsplit(str(path)).scheme.lower() in {'http', 'https'}


def source_path(path):
    return urlsplit(str(path)).path if is_remote(path) else str(path)


def sidecar(path, suffix, replace=False):
    if not is_remote(path):
        from pathlib import Path
        return str(Path(path).with_suffix(suffix)) if replace else str(path)+suffix
    parts = urlsplit(str(path))
    name = parts.path
    if replace:
        name = name.rsplit('.', 1)[0]
    return urlunsplit(parts._replace(path=name + suffix))


def _request(url, headers=None, method=None):
    return urlopen(Request(url, headers={'Accept-Encoding': 'identity',
                                       'User-Agent': 'poreckit-gui', **(headers or {})}, method=method), timeout=30)


def metadata(url):
    previous = _metadata.get(url)
    if previous and time.monotonic() - previous[0] < 10:
        return previous[1]
    try:
        response = _request(url, method='HEAD')
    except HTTPError as exc:
        if exc.code not in {405, 501, 403}:
            raise
        response = _request(url, {'Range': 'bytes=0-0'})
    with response:
        headers = response.headers
        content_range = headers.get('Content-Range', '')
        match = re.fullmatch(r'bytes \d+-\d+/(\d+)', content_range)
        size = int(match[1]) if match else int(headers.get('Content-Length', -1))
        if size < 0:
            raise ValueError('Remote server must provide a file size for indexed access')
        result = (size, headers.get('ETag', ''), headers.get('Last-Modified', ''))
    if len(_metadata) >= 128:
        _metadata.pop(next(iter(_metadata)))
    _metadata[url] = (time.monotonic(), result)
    return result


def exists(path):
    if not is_remote(path):
        from pathlib import Path
        return Path(path).is_file()
    try:
        metadata(str(path))
        return True
    except HTTPError as exc:
        if exc.code in {404, 410}:
            return False
        raise


def clear_cache():
    global _cache_size
    _cache.clear()
    _metadata.clear()
    _cache_size = 0


class RangeFile(io.RawIOBase):
    def __init__(self, url):
        super().__init__()
        self.url = str(url)
        self.length, self.etag, self.modified = metadata(self.url)
        self.position = 0
        self.epoch = _generation() if not (self.etag or self.modified) else 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=io.SEEK_SET):
        self._checkClosed()
        position = offset + (self.position if whence == io.SEEK_CUR else self.length if whence == io.SEEK_END else 0)
        if whence not in {io.SEEK_SET, io.SEEK_CUR, io.SEEK_END} or position < 0:
            raise ValueError('Invalid remote seek')
        self.position = position
        return position

    def _block(self, start):
        global _cache_size
        key = self.url, self.length, self.etag, self.modified, self.epoch, start
        if key in _cache:
            _cache.move_to_end(key)
            return _cache[key]
        stop = min(self.length, start + BLOCK_SIZE) - 1
        headers = {'Range': f'bytes={start}-{stop}'}
        if self.etag and not self.etag.startswith('W/'):
            headers['If-Match'] = self.etag
        elif self.modified:
            headers['If-Unmodified-Since'] = self.modified
        with _request(self.url, headers) as response:
            if response.headers.get('Content-Encoding', 'identity') not in {'identity', ''}:
                raise ValueError('Remote indexed data must not use HTTP content compression')
            status = response.status
            if status == 206:
                expected = f'bytes {start}-{stop}/{self.length}'
                if response.headers.get('Content-Range') != expected:
                    raise ValueError('Server returned the wrong HTTP byte range')
                content = response.read(stop-start+2)
                if len(content) != stop-start+1:
                    raise ValueError('Truncated HTTP range response')
            elif status == 200 and start == 0 and self.length <= BLOCK_SIZE:
                content = response.read(self.length+1)
                if len(content) != self.length:
                    raise ValueError('Remote file length changed during reading')
            else:
                raise ValueError('Server does not support HTTP Range access for this file')
            if self.etag and response.headers.get('ETag', self.etag) != self.etag:
                raise ValueError('Remote file changed during reading; reopen the source')
        while _cache and _cache_size+len(content) > CACHE_BYTES:
            _cache_size -= len(_cache.popitem(last=False)[1])
        _cache[key] = content
        _cache_size += len(content)
        return content

    def read(self, size=-1):
        self._checkClosed()
        if size is None or size < 0:
            size = self.length-self.position
        size = max(0, min(size, self.length-self.position))
        if size > 256*1024*1024:
            raise ValueError('Read a large remote source in bounded chunks')
        result = bytearray()
        while len(result) < size:
            start = self.position//BLOCK_SIZE*BLOCK_SIZE
            data = self._block(start)
            within = self.position-start
            count = min(size-len(result), len(data)-within)
            if count <= 0:
                raise ValueError('Invalid HTTP range block')
            result.extend(data[within:within+count])
            self.position += count
        return bytes(result)

    def readinto(self, target):
        content = self.read(len(target))
        target[:len(content)] = content
        return len(content)


def open_binary(path):
    return RangeFile(path) if is_remote(path) else open(path, 'rb')


def _generation():
    """Ticks every ~10s so same-size responses without validators eventually
    invalidate caches keyed on metadata alone."""
    return int(time.monotonic() // 10)


def identity(path):
    if is_remote(path):
        url = str(path)
        size, etag, modified = metadata(url)
        generation = _generation() if not (etag or modified) else 0
        return (url, size, etag, modified, generation)
    from pathlib import Path
    p = Path(path).resolve()
    stat = p.stat()
    return str(p), stat.st_size, stat.st_mtime_ns
