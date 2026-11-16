"""Cancellable data queries and portable export for native genomic views."""
import json
import os
import sys
import tempfile
from pathlib import Path

from PySide6.QtCore import QProcess, QProcessEnvironment, QRectF, Qt
from PySide6.QtGui import QImage, QPainter, QPdfWriter, QPageSize
from PySide6.QtSvg import QSvgGenerator
from PySide6.QtWidgets import QWidget, QFileDialog

from .jobs import JobProcess


class DataView(QWidget):
    def __init__(self, browser):
        super().__init__()
        self.browser = browser
        self.temp = tempfile.TemporaryDirectory(prefix='poreckit-data-')
        self.process = JobProcess(self)
        self.process.readyReadStandardOutput.connect(self.read_log)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.failed)
        self.job = None
        self.pending = []
        self.completing = False
        self.cancelled = False
        self.log = ''
        self.progress_buffer = ''

    def busy(self):
        return self.job is not None

    def request(self, job, callback):
        if self.busy() or self.completing:
            # Coalesce navigation, but never discard a file import or export.
            if self.pending and self.pending[-1][0]['action'] == job['action'] and job['action'] in {'sequence_query', 'assembly_preview'}:
                self.pending[-1] = job, callback
            else:
                self.pending.append((job, callback))
            self.status.setText('Query queued · Cancel stops the worker')
            return
        self.log = ''
        self.progress_buffer = ''
        self.cancelled = False
        self.job = job, callback
        folder = Path(self.temp.name)
        job_file, self.result_path = folder/'job.json', folder/'result.json'
        self.result_path.unlink(missing_ok=True)
        job_file.write_text(json.dumps(job), encoding='utf8')
        env = QProcessEnvironment.systemEnvironment()
        env.insert('PYTHONUNBUFFERED','1')
        env.insert('MPLCONFIGDIR',str(folder/'matplotlib'))
        self.process.setProcessEnvironment(env)
        args = ['--worker'] if getattr(sys,'frozen',False) else ['-m','poreckit_gui.worker']
        self.status.setText('Reading data…')
        self.process.start(sys.executable,[*args,str(job_file),str(self.result_path)])

    def read_log(self):
        value=bytes(self.process.readAllStandardOutput()).decode('utf8','replace')
        self.log=(self.log+value)[-131072:]
        self.progress_buffer+=value
        while '\n' in self.progress_buffer:
            line,self.progress_buffer=self.progress_buffer.split('\n',1)
            if line.startswith('@poreckit-progress:'):
                try:
                    self.status.setText(json.loads(line.split(':',1)[1])['message'])
                except (ValueError,KeyError,TypeError):
                    pass

    def failed(self,error):
        if error == QProcess.ProcessError.FailedToStart:
            self.log += self.process.errorString()
            self.finished(-1,None)

    def finished(self,code,_):
        self.read_log()
        previous,self.job = self.job,None
        if previous is None:
            return
        self.completing = True
        try:
            if self.cancelled:
                self.status.setText('Cancelled · Previous results retained')
            elif code:
                self.status.setText('Query failed · '+self.log.strip().split('\n')[-1])
                self.status.setToolTip(self.log)
            elif not (previous[0]['action'] in {'sequence_query', 'assembly_preview'} and self.pending):
                result=json.loads(self.result_path.read_text(encoding='utf8'))
                previous[1](result)
        except (ValueError,OSError,KeyError,TypeError) as exc:
            self.status.setText(str(exc))
        finally:
            self.completing = False
        if self.pending:
            self.request(*self.pending.pop(0))

    def cancel(self):
        self.pending=[]
        if self.busy():
            self.cancelled = True
            self.process.kill()
            self.status.setText('Cancelling…')

    def shutdown(self):
        self.pending=[]
        self.process.shutdown()
        from .sessionio import cleanup_temp; cleanup_temp(self.temp.name)

    def export(self):
        path,_=QFileDialog.getSaveFileName(self,'Export view','view.svg','SVG (*.svg);;PNG (*.png);;PDF (*.pdf)')
        if path:
            try:
                export_scene(self.canvas.scene(),path)
                self.status.setText('Exported '+path)
            except (ValueError,OSError) as exc:
                self.status.setText(str(exc))


def export_scene(scene,path):
    from PySide6.QtCore import QSize, QSizeF
    bounds=scene.sceneRect()
    size=QSize(max(1,int(bounds.width())),max(1,int(bounds.height())))
    suffix=Path(path).suffix.lower()
    if suffix=='.svg':
        device=QSvgGenerator()
        device.setFileName(str(path))
        device.setSize(size)
        device.setViewBox(bounds)
    elif suffix=='.pdf':
        device=QPdfWriter(str(path))
        device.setResolution(96)
        device.setPageSize(QPageSize(QSizeF(size.width()*.75, size.height()*.75),QPageSize.Unit.Point))
    elif suffix=='.png':
        if size.width()*size.height()>100_000_000:
            raise ValueError('Image is too large; export SVG or reduce visible reads')
        device=QImage(size,QImage.Format.Format_ARGB32)
        device.fill(Qt.GlobalColor.white)
    else:
        raise ValueError('Choose .svg, .png or .pdf')
    painter=QPainter(device)
    try:
        if not painter.isActive():
            raise OSError('Cannot open export destination')
        if suffix=='.pdf':
            target=QRectF(0,0,device.width(),device.height())
        else:
            target=QRectF(0,0,size.width(),size.height())
        scene.render(painter,target,bounds,Qt.AspectRatioMode.KeepAspectRatio)
    finally:
        painter.end()
    if suffix=='.png' and not device.save(str(path)):
        raise OSError('Unable to write PNG')
