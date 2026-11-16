"""Independent Qt signal, alignment, variant and reference browser."""
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGraphicsEllipseItem,
    QGraphicsRectItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from .dataview import DataView
from .model import Region
from .remote import is_remote, source_path
from .sequenceio import file_kind

BASE_COLORS={'A':'#279c50','C':'#3376ce','G':'#bd9415','T':'#cf4545','N':'#888888','=':'#888888'}


class TrackOptions(QDialog):
    def __init__(self,track,parent=None):
        super().__init__(parent)
        self.setWindowTitle('Track settings')
        self.track=track
        layout=QVBoxLayout(self)
        form=QFormLayout()
        layout.addLayout(form)
        self.fields: dict[str, QWidget] = {}
        def field(name,label,widget):
            self.fields[name]=widget
            form.addRow(label,widget)
            return widget
        field('name','Name',QLineEdit(track.get('name','')))
        field('color','Color',QLineEdit(track.get('color','#397ca8')))
        if track['kind']=='signal':
            for key,label,values in [('mode','Style',['fill','line','bar','points']),('scale','Y axis',['auto','fixed','shared'])]:
                box=QComboBox()
                box.addItems(values)
                box.setCurrentText(track.get(key,values[0]))
                field(key,label,box)
            field('group','Shared axis group',QLineEdit(track.get('group','signals')))
            for key,label,default in [('ymin','Minimum',0),('ymax','Maximum',1)]:
                spin=QDoubleSpinBox()
                spin.setRange(-1e15,1e15)
                spin.setDecimals(5)
                spin.setValue(float(track.get(key,default)))
                field(key,label,spin)
        elif track['kind']=='alignment':
            for key,label,default,maximum in [('mapq','Minimum mapping quality',0,255),('baseq','Minimum base quality (mismatch colors)',0,93),('max_reads','Displayed reads',1000,10000),('record_limit','Query record limit (0 = all)',0,2147483647)]:
                spin=QSpinBox()  # type: ignore[assignment]  # only setRange/setValue are shared with QDoubleSpinBox
                spin.setRange(0 if key!='max_reads' else 1,maximum)
                spin.setValue(int(track.get(key,default)))
                field(key,label,spin)
            for key,label,default in [('duplicates','Include duplicates',False),('secondary','Include secondary',False),('supplementary','Include supplementary',True),('pairs','Connect read pairs',False),('junction_arcs','Show splice junctions',True)]:
                check=QCheckBox()
                check.setChecked(track.get(key,default))
                field(key,label,check)
            for key,label,values in [('sort','Sort reads',['position','name','strand','mapq']),('read_color','Color reads by',['strand','mapping quality'])]:
                box=QComboBox()
                box.addItems(values)
                box.setCurrentText(track.get(key,values[0]))
                field(key,label,box)
        elif track['kind']=='variants':
            check=QCheckBox()
            check.setChecked(track.get('pass_only',False))
            field('pass_only','PASS / unfiltered only',check)
        if track['kind'] in {'alignment','variants'}:
            field('index','Index file / URL',QLineEdit(track.get('index') or ''))
        buttons=QDialogButtonBox(QDialogButtonBox.StandardButton.Ok|QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def accept(self):
        values: dict[str, str | int | float | bool] = {}
        for key,widget in self.fields.items():
            value: str | bool | int | float
            if isinstance(widget,QLineEdit):
                value=widget.text()
            elif isinstance(widget,QComboBox):
                value=widget.currentText()
            elif isinstance(widget,QCheckBox):
                value=widget.isChecked()
            else:
                assert isinstance(widget,(QDoubleSpinBox,QSpinBox))  # only spinboxes remain
                value=widget.value()
            values[key]=value
        if not QColor(values['color']).isValid():
            QMessageBox.warning(self,'Track settings','Enter a valid color')
            return
        if values.get('scale')=='fixed' and float(values['ymin'])>=float(values['ymax']):
            QMessageBox.warning(self,'Track settings','Minimum must be below maximum')
            return
        self.track.update(values)
        super().accept()


class TrackCanvas(QGraphicsView):
    result: dict
    navigate=Signal(float,float)
    pan=Signal(float)
    selected=Signal(object)
    hovered=Signal(str)

    def __init__(self):
        super().__init__()
        self.setScene(QGraphicsScene(self))
        self.setAlignment(Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignTop)
        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setMouseTracking(True)
        self.result=None  # type: ignore[assignment]  # cleared until a query result arrives
        self.drag=None
        self.margin=145
        self.plot_width=800
        self.resize_timer=QTimer(self)
        self.resize_timer.setSingleShot(True)
        self.resize_timer.setInterval(120)
        self.resize_timer.timeout.connect(self.redraw)

    def resizeEvent(self,event):
        super().resizeEvent(event)
        self.resize_timer.start()

    def x(self,pos):  # type: ignore[override]  # coordinate helper, not QWidget.x()
        r=self.region
        return self.margin+self.plot_width*(max(r.start,min(r.end,pos))-r.start)/(r.end-r.start)

    def text(self,text,x,y,color='#222222',size=9):
        item=self.scene().addSimpleText(str(text))
        font=item.font()
        font.setPointSize(size)
        item.setFont(font)
        item.setBrush(QColor(color))
        item.setPos(x,y)
        return item

    def rect(self,start,end,y,height,color,record=None):  # type: ignore[override]  # drawing helper, not QWidget.rect()
        item=self.scene().addRect(self.x(start),y,max(1,self.x(end)-self.x(start)),height,QPen(Qt.PenStyle.NoPen),QColor(color))
        if record is not None:
            item.setData(0,record)
            item.setToolTip(json.dumps(record,ensure_ascii=False,indent=2)[:12000])
        return item

    def show_result(self,result):
        self.result=result
        self.redraw()

    def zoom_by(self,factor):
        self.navigate.emit(1/factor,.5)

    def fit(self):
        self.resetTransform()
        self.redraw()

    def actual_size(self):
        self.resetTransform()

    def redraw(self):
        if not self.result:
            return
        self.scene().clear()
        self.region=Region.parse(self.result['region'])
        self.plot_width=max(400,self.viewport().width()-self.margin-30)
        r=self.region
        span=r.end-r.start
        self.text(r.text(),5,0)
        for i in range(6):
            pos=r.start+i*span/5
            x=self.x(pos)
            self.text(f'{int(pos):,}',x-15,20,size=8)
            self.scene().addLine(x,37,x,44,QPen(QColor('#888888')))
        y=48
        seq=self.result.get('sequence','')
        if self.result.get('reference'):
            self.text('Reference',5,y)
            if seq:
                for i,base in enumerate(seq):
                    self.rect(r.start+i,r.start+i+1,y,18,BASE_COLORS.get(base,'#888888'),dict(position=r.start+i,base=base))
                    if self.plot_width/span>10:
                        self.text(base,self.x(r.start+i)+1,y,'#ffffff',8)
            else:
                self.text(self.result.get('reference_error') or 'Zoom to ≤20,000 bases to inspect reference sequence',self.margin,y)
            y+=30
        for entry in self.result['tracks']:
            track=entry['track']
            self.scene().addLine(0,y,self.margin+self.plot_width,y,QPen(QColor('#cccccc')))
            label=self.text(track.get('name',Path(track['path']).name)[:21],5,y+4)
            label.setToolTip(track['path'])
            if entry.get('error'):
                self.text(entry['error'][:130],self.margin,y+5,'#a82828')
                y+=35
                continue
            data=entry['data']
            if track['kind']=='signal':
                y=self.signal(track,data,y)
            elif track['kind']=='alignment':
                y=self.alignments(track,data,y,seq)
            else:
                lanes: list = []
                hidden=0
                for record in sorted(data['records'],key=lambda v:v['start']):
                    row=next((i for i,end in enumerate(lanes) if end<=record['start']),len(lanes))
                    if row>=80:
                        hidden+=1
                        continue
                    if row==len(lanes):
                        lanes.append(0)
                    lanes[row]=record['end']+span*0.004
                    top=y+24+row*18
                    color=track.get('color','#397ca8')
                    if track['kind']=='variants':
                        alt=record['alt']
                        color='#a04fc4' if len(alt)!=len(record['ref']) or alt.startswith('<') else BASE_COLORS.get(alt,'#397ca8')
                    if record.get('blocks'):
                        self.scene().addLine(self.x(record['start']),top+6,self.x(record['end']),top+6,QPen(QColor(color)))
                        for a,b in record['blocks']:
                            self.rect(a,b,top+2,8,color,record)
                        for a,b in record.get('cds',[]):
                            self.rect(a,b,top,12,color,record)
                        direction='‹' if record.get('strand')=='-' else '›'
                        for x in range(int(self.x(record['start']))+20,int(self.x(record['end']))-5,45):
                            self.text(direction,x,top-3,color,9)
                    else:
                        self.rect(record['start'],record['end'],top,12,color,record)
                    if self.x(record['end'])-self.x(record['start'])>50:
                        self.text(record.get('name','')[:45],self.x(record['start']),top-2,size=7)
                y+=max(55,25+len(lanes)*18)
                if data.get('truncated') or hidden:
                    self.text(f"Display limited · {hidden} overlapping features hidden",self.margin,y,'#a65c00')
                    y+=24
        self.scene().setSceneRect(0,0,self.margin+self.plot_width+20,max(y+10,200))

    def signal(self,track,data,y):
        values=data['values']
        if not any(v is not None for v in values):
            self.text('No signal in this region',self.margin,y+40,'#777777')
        lo,hi=data['ymin'],data['ymax']
        if hi<=lo:
            hi=lo+1
        height=90
        self.text(f'{hi:g}',105,y+18,size=7)
        self.text(f'{lo:g}',105,y+height+10,size=7)
        def sy(value):
            return y+22+height*(1-(max(lo,min(hi,value))-lo)/(hi-lo))
        zero=sy(0)
        self.scene().addLine(self.margin,zero,self.margin+self.plot_width,zero,QPen(QColor('#aaaaaa')))
        path=QPainterPath()
        previous=False
        color=QColor(track.get('color','#397ca8'))
        mode=track.get('mode','fill')
        for i,value in enumerate(values):
            if value is None:
                previous=False
                continue
            a=self.region.start+i*(self.region.end-self.region.start)/len(values)
            b=self.region.start+(i+1)*(self.region.end-self.region.start)/len(values)
            x=self.x((a+b)/2)
            if mode in {'fill','bar'}:
                item: QGraphicsRectItem | QGraphicsEllipseItem = self.scene().addRect(self.x(a),min(zero,sy(value)),max(.3,self.x(b)-self.x(a)),max(.3,abs(zero-sy(value))),QPen(Qt.PenStyle.NoPen),color)
            elif mode=='points':
                item=self.scene().addEllipse(x-1,sy(value)-1,2,2,QPen(color),color)
            else:
                if previous:
                    path.lineTo(x,sy(value))
                else:
                    path.moveTo(x,sy(value))
                previous=True
                continue
            item.setData(0,dict(start=int(a),end=int(b),value=value,track=track['name']))
            item.setToolTip(f'{int(a):,}–{int(b):,}: {value:g}')
        if mode=='line':
            self.scene().addPath(path,QPen(color,1.2)).setToolTip(track['name']+' · line signal')
        return y+height+40

    def alignments(self,track,data,y,sequence):
        self.text(f"{data['accepted']:,} reads · mean depth",self.margin,y+2,size=8)
        values=data['coverage']
        covmax=max(values,default=1) or 1
        covpath=QPainterPath()
        covpath.moveTo(self.margin,y+70)
        for i,value in enumerate(values):
            covpath.lineTo(self.margin+self.plot_width*(i+.5)/len(values),y+70-40*value/covmax)
        covpath.lineTo(self.margin+self.plot_width,y+70)
        self.scene().addPath(covpath,QPen(QColor('#78909c')),QColor('#90a4ae'))
        self.text(f'{covmax:g}',105,y+27,size=7)
        y+=78
        if track.get('junction_arcs',True) and data.get('junctions'):
            for junction in data['junctions']:
                if junction['end']<=self.region.start or junction['start']>=self.region.end:
                    continue
                a,b=self.x(junction['start']),self.x(junction['end'])
                path=QPainterPath()
                path.moveTo(a,y+30)
                path.cubicTo(a,y-8,b,y-8,b,y+30)
                item=self.scene().addPath(path,QPen(QColor('#6f55a3'),1.2))
                item.setData(0,junction)
                item.setToolTip(f"Splice junction {junction['start']:,}–{junction['end']:,}: {junction['count']} reads")
                self.text(junction['count'],(a+b)/2,y-2,'#6f55a3',7)
            y+=40
        lanes: list = []
        pairs: dict = {}
        hidden=0
        for record in data['records']:
            row=next((i for i,end in enumerate(lanes) if end<=record['start']),len(lanes))
            if row>=200:
                hidden+=1
                continue
            if row==len(lanes):
                lanes.append(0)
            lanes[row]=record['end']+(self.region.end-self.region.start)*.002
            top=y+row*16
            if track.get('read_color')=='mapping quality':
                grey=max(40,210-min(record['mapq'],60)*2)
                color=QColor(grey,grey,grey).name()
            else:
                color='#4e91be' if record['flag'] & 16 else '#c97878'
            self.scene().addLine(self.x(record['start']),top+5,self.x(record['end']),top+5,QPen(QColor(color)))
            for a,b,q in record['blocks']:
                self.rect(a,b,top,10,color,record)
                if sequence and self.plot_width/(self.region.end-self.region.start)>=2:
                    lo,hi=max(a,self.region.start),min(b,self.region.end)
                    for pos in range(lo,hi):
                        qi=q+pos-a
                        if qi>=len(record['sequence']):
                            continue
                        base=record['sequence'][qi]
                        quality=record['qualities'][qi] if qi<len(record['qualities']) else 255
                        if quality < track.get('baseq',0):
                            continue
                        if base not in {'=',sequence[pos-self.region.start]}:
                            self.rect(pos,pos+1,top,10,BASE_COLORS.get(base,'#999999'),record)
                            if self.plot_width/(self.region.end-self.region.start)>10:
                                self.text(base,self.x(pos),top-2,'#ffffff',7)
            for op,pos,n,_q in record['edits']:
                if self.region.start<=pos<self.region.end:
                    if op in 'IS':
                        self.scene().addLine(self.x(pos),top-2,self.x(pos),top+12,QPen(QColor('#a531bd' if op=='I' else '#e2a22c'),2)).setToolTip(f'{op}: {n} bases')
                    elif op=='D':
                        self.scene().addLine(self.x(pos),top+5,self.x(pos+n),top+5,QPen(QColor('#111111'),1.4)).setToolTip(f'Deletion: {n} bases')
            if track.get('pairs'):
                name=record['name']
                if name in pairs:
                    px,py=pairs[name]
                    self.scene().addLine(px,py,self.x(record['start']),top+5,QPen(QColor('#999999'),.6,Qt.PenStyle.DashLine))
                else:
                    pairs[name]=self.x(record['end']),top+5
        y+=max(25,len(lanes)*16)+4
        if data.get('truncated') or data.get('reads_hidden') or hidden:
            message=f"Read display limited · {data.get('reads_hidden',0)+hidden:,} reads hidden"
            if data.get('truncated'):
                message+=' · coverage is partial (query record limit)'
            self.text(message,self.margin,y,'#a65c00',8)
            y+=25
        return y

    def wheelEvent(self,event):
        if self.result and not event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            fraction=max(0,min(1,(event.position().x()-self.margin)/self.plot_width))
            self.navigate.emit(.5 if event.angleDelta().y()>0 else 2.,fraction)
            event.accept()
        else:
            super().wheelEvent(event)

    def mousePressEvent(self,event):
        if event.button()==Qt.MouseButton.LeftButton:
            self.drag=event.position()
        super().mousePressEvent(event)

    def mouseMoveEvent(self,event):
        if self.result:
            p=self.mapToScene(event.position().toPoint())
            fraction=max(0,min(1,(p.x()-self.margin)/self.plot_width))
            pos=self.region.start+int(fraction*(self.region.end-self.region.start))
            messages=[f'{self.region.chrom}:{pos:,} (zero-based)']
            for entry in self.result['tracks']:
                if entry['track']['kind']=='signal' and 'data' in entry:
                    values=entry['data']['values']
                    value=values[min(len(values)-1,int(fraction*len(values)))] if values else None
                    messages.append(entry['track']['name']+': '+('missing' if value is None else f'{value:g}'))
            self.hovered.emit(' · '.join(messages))
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self,event):
        if self.drag is not None and self.result:
            dx=event.position().x()-self.drag.x()
            if abs(dx)>5:
                self.pan.emit(-dx/self.plot_width)
            else:
                for item in self.items(event.position().toPoint()):
                    record=item.data(0)
                    if record is not None:
                        self.selected.emit(record)
                        break
        self.drag=None
        super().mouseReleaseEvent(event)


class SequenceBrowser(DataView):
    canvas: TrackCanvas
    hover: QLabel
    result: dict
    def __init__(self,browser):
        super().__init__(browser)
        self.chromsizes={}
        self.reference=''
        self.region=None
        self.history=[]
        self.history_index=-1
        self.result=None  # type: ignore[assignment]  # cleared until a query result arrives
        self.loading=False
        self.current_session=None
        layout=QVBoxLayout(self)
        nav=QHBoxLayout()
        self.chrom=QComboBox()
        self.chrom.currentTextChanged.connect(self.chrom_changed)
        self.locus=QLineEdit()
        self.locus.setPlaceholderText('chr1:0-100,000 · zero-based, half-open')
        self.locus.returnPressed.connect(self.go)
        for label,callback in [('Open tracks…',self.open_files),('Open URL…',self.open_url),('Reference…',self.open_reference),('Back',lambda:self.travel(-1)),('Forward',lambda:self.travel(1))]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            nav.addWidget(b)
        nav.addWidget(self.chrom)
        nav.addWidget(self.locus,1)
        for label,callback in [('Go',self.go),('−',lambda:self.zoom(2)),('+',lambda:self.zoom(.5)),('←',lambda:self.move(-.5)),('→',lambda:self.move(.5)),('Whole chr',self.home)]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            nav.addWidget(b)
        layout.addLayout(nav)
        tools=QHBoxLayout()
        for label,callback in [('Use map locus',self.use_map),('Show locus in map',self.show_map),('Find gene…',self.find_gene),('Export…',self.export),('Save session…',self.save_session),('Open session…',self.load_session),('Cancel',self.cancel)]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            tools.addWidget(b)
        self.link=QCheckBox('Follow map navigation')
        tools.addWidget(self.link)
        self.track_link=QCheckBox('Link track windows')
        tools.addWidget(self.track_link)
        tools.addStretch()
        layout.addLayout(tools)
        splitter=QSplitter()
        left=QVBoxLayout()
        from PySide6.QtWidgets import QWidget
        panel=QWidget()
        panel.setLayout(left)
        left.addWidget(QLabel('Tracks · Drag to reorder; double-click for settings'))
        self.tracks=QListWidget()
        self.tracks.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.tracks.itemChanged.connect(lambda _:self.query())
        self.tracks.model().rowsMoved.connect(lambda *args:self.query())
        self.tracks.itemDoubleClicked.connect(self.settings)
        left.addWidget(self.tracks,1)
        for label,callback in [('Settings…',self.settings),('Remove',self.remove),('Build BAM index',self.build_index),('Duplicate in new window',self.new_window)]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            left.addWidget(b)
        panel.setMaximumWidth(290)
        splitter.addWidget(panel)
        self.canvas=TrackCanvas()
        self.canvas.navigate.connect(self.zoom)
        self.canvas.pan.connect(self.move)
        self.canvas.selected.connect(self.details)
        self.canvas.hovered.connect(lambda text:self.hover.setText(text))
        splitter.addWidget(self.canvas)
        splitter.setStretchFactor(1,1)
        layout.addWidget(splitter,1)
        self.hover=QLabel('Open bigWig, bedGraph, BAM/CRAM/SAM, VCF/BCF, FASTA, BED or GTF/GFF')
        self.hover.setWordWrap(True)
        layout.addWidget(self.hover)
        self.detail=QPlainTextEdit()
        self.detail.setReadOnly(True)
        self.detail.setMaximumHeight(130)
        self.detail.setPlaceholderText('Click a read, variant or feature for details. Sequences and genotype fields can be copied here.')
        layout.addWidget(self.detail)
        self.status=QLabel('Ready · No contact map required')
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.setAcceptDrops(True)

    def open_files(self):
        paths,_=QFileDialog.getOpenFileNames(self,'Open genomic tracks','',
            'Genomic tracks (*.bw *.bigwig *.bigWig *.bg *.bedgraph *.bedgraph.gz *.bam *.cram *.sam *.vcf *.vcf.gz *.bcf *.bed *.bed.gz *.gtf *.gtf.gz *.gff *.gff3 *.fa *.fasta *.fna *.fa.gz);;All files (*)')
        if paths:
            self.add_files(paths)

    def reset(self):
        self.cancel()
        self.loading=True
        self.tracks.clear()
        self.loading=False
        self.reference=''
        self.chromsizes={}
        self.region=None
        self.result=None  # type: ignore[assignment]  # cleared until a query result arrives
        self.history=[]
        self.history_index=-1
        self.canvas.scene().clear()
        self.canvas.result=None  # type: ignore[assignment]  # cleared until a query result arrives
        self.sync_chroms()
        self.locus.clear()

    def open_reference(self):
        path,_=QFileDialog.getOpenFileName(self,'Reference genome','','FASTA (*.fa *.fasta *.fna *.fa.gz *.fasta.gz)')
        if path:
            self.add_files([path])

    def open_url(self):
        dialog=QDialog(self)
        dialog.setWindowTitle('Open remote track')
        layout=QVBoxLayout(dialog)
        form=QFormLayout()
        url=QLineEdit()
        url.setPlaceholderText('https://…/sample.bw, .bam, .vcf.gz or .fa')
        index=QLineEdit()
        index.setPlaceholderText('Optional custom BAI/CSI/TBI URL')
        form.addRow('Data URL',url)
        form.addRow('Index URL',index)
        layout.addLayout(form)
        buttons=QDialogButtonBox(QDialogButtonBox.StandardButton.Open|QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.resize(600,150)
        if dialog.exec():
            source=url.text().strip()
            side=index.text().strip()
            if not is_remote(source) or side and not is_remote(side):
                self.status.setText('Enter HTTP or HTTPS URLs')
                return
            self.add_files([source],indexes={source:side} if side else {})

    def add_files(self,paths,indexes=None):
        self.request(dict(action='sequence_inspect',paths=[str(p) if is_remote(p) else str(Path(p).resolve()) for p in paths],
                          indexes=indexes or {},chromsizes=list(self.chromsizes.items())),self.add_inspected)

    def add_inspected(self,result):
        for track in result['tracks']:
            if track.get('error'):
                self.status.setText(track['path']+': '+track['error'])
                continue
            if self.chromsizes and track['kind']!='reference':
                # Unknown contigs can coexist, but conflicting lengths must not
                # silently shorten an explicitly loaded reference genome.
                for chrom,length in track['chromsizes']:
                    if chrom not in self.chromsizes:
                        self.chromsizes[chrom]=length
                    elif not self.reference:
                        self.chromsizes[chrom]=max(length,self.chromsizes[chrom])
            else:
                self.chromsizes=dict(track['chromsizes'])
            if track['kind']=='reference':
                self.reference=track['path']
            else:
                self.append(track)
        self.sync_chroms()
        if not self.region and self.chromsizes:
            name=next(iter(self.chromsizes))
            self.region=Region(name,0,min(self.chromsizes[name],100000))
        if self.region:
            if self.region.chrom not in self.chromsizes:
                name=next(iter(self.chromsizes))
                self.region=Region(name,0,min(self.chromsizes[name],100000))
            self.navigate(self.region)

    def append(self,track):
        item=QListWidgetItem(track.get('name',Path(track['path']).name))
        item.setData(Qt.ItemDataRole.UserRole,dict(track))
        item.setFlags(item.flags()|Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(Qt.CheckState.Checked if track.get('visible',True) else Qt.CheckState.Unchecked)
        self.loading=True
        self.tracks.addItem(item)
        self.loading=False

    def track_data(self):
        return [{**self.tracks.item(i).data(Qt.ItemDataRole.UserRole),'visible':self.tracks.item(i).checkState()==Qt.CheckState.Checked} for i in range(self.tracks.count())]

    def sync_chroms(self):
        self.chrom.blockSignals(True)
        self.chrom.clear()
        self.chrom.addItems(list(self.chromsizes))
        if self.region:
            self.chrom.setCurrentText(self.region.chrom)
        self.chrom.blockSignals(False)

    def navigate(self,region,remember=True,linked=False):
        if region.chrom not in self.chromsizes:
            from .sequenceio import chromosome
            try:
                name=chromosome(self.chromsizes,region.chrom)
                region=Region(name,region.start,region.end)
            except ValueError as exc:
                self.status.setText(str(exc))
                return
        length=self.chromsizes[region.chrom]
        region=region.window(region.start,region.end-region.start,length)
        self.region=region
        self.locus.setText(region.text())
        self.chrom.blockSignals(True)
        self.chrom.setCurrentText(region.chrom)
        self.chrom.blockSignals(False)
        if remember and (self.history_index<0 or self.history[self.history_index]!=region):
            self.history=self.history[:self.history_index+1]+[region]
            self.history_index=len(self.history)-1
        self.query()
        if self.track_link.isChecked() and not linked:
            from .sequenceio import chromosome
            from .windows import WindowController
            for controller in list(WindowController.instances):
                target=controller.browser.sequence_browser
                if target is self or not target.track_link.isChecked():
                    continue
                try:
                    name=chromosome(target.chromsizes,region.chrom)
                    if target.chromsizes[name]!=self.chromsizes[region.chrom]:
                        target.status.setText('Track navigation not linked: chromosome lengths differ')
                        continue
                    target.navigate(Region(name,region.start,region.end),linked=True)
                except ValueError:
                    target.status.setText('Track navigation not linked: chromosome is absent')

    def go(self):
        try:
            text=self.locus.text().strip()
            if text in self.chromsizes:
                self.navigate(Region(text,0,self.chromsizes[text]))
            else:
                self.navigate(Region.parse(text))
        except ValueError as exc:
            self.status.setText(str(exc))

    def chrom_changed(self,name):
        if name and self.region and name in self.chromsizes:
            self.navigate(Region(name,0,min(self.chromsizes[name],self.region.end-self.region.start)))

    def zoom(self,factor,fraction=.5):
        if self.region:
            r=self.region
            span=r.end-r.start
            self.navigate(r.window(r.start+span*fraction-span*factor*fraction,span*factor,self.chromsizes.get(r.chrom)))

    def move(self,fraction):  # type: ignore[override]  # pan helper, not QWidget.move()
        if self.region:
            self.navigate(self.region.pan(fraction,self.chromsizes.get(self.region.chrom)))

    def home(self):
        if self.region:
            self.navigate(Region(self.region.chrom,0,self.chromsizes[self.region.chrom]))

    def travel(self,delta):
        target=self.history_index+delta
        if 0<=target<len(self.history):
            self.history_index=target
            self.navigate(self.history[target],False)

    def query(self):
        if self.loading or not self.region:
            return
        job=dict(action='sequence_query',region=self.region.text(),reference=self.reference,
                 tracks=self.track_data(),bins=max(100,self.canvas.viewport().width()))
        self.request(job,self.ready)

    def ready(self,result):
        self.result=result
        self.canvas.show_result(result)
        errors=sum('error' in t for t in result['tracks'])
        self.status.setText(f"{result['region']} · {len(result['tracks'])} tracks · {errors} track errors"+(' · Reference: '+Path(self.reference).name if self.reference else ''))
        self.status.setToolTip('\n'.join(t['error'] for t in result['tracks'] if 'error' in t))

    def settings(self,*args):
        item=self.tracks.currentItem()
        if item:
            track=dict(item.data(Qt.ItemDataRole.UserRole))
            if TrackOptions(track,self).exec():
                self.loading=True
                item.setData(Qt.ItemDataRole.UserRole,track)
                item.setText(track['name'])
                self.loading=False
                self.query()

    def remove(self):
        item=self.tracks.takeItem(self.tracks.currentRow())
        if item:
            self.query()

    def details(self,record):
        self.detail.setPlainText(json.dumps(record,ensure_ascii=False,indent=2))

    def new_window(self):
        import shutil
        import uuid
        new=self.browser.windows.create()
        if new is None:
            return
        destination=new.sequence_browser
        def persistent(path):
            if not path:
                return ''
            if is_remote(path):
                return path
            p=Path(path).resolve()
            roots=(Path(self.temp.name).resolve(),Path(self.browser.temp.name).resolve())
            if not any(p.is_relative_to(root) for root in roots):
                return str(p)
            target=Path(destination.temp.name)/(uuid.uuid4().hex+'-'+p.name)
            shutil.copyfile(p,target)
            for suffix in ('.fai','.gzi','.tbi','.bai','.csi'):
                sidecar=Path(str(p)+suffix)
                if sidecar.is_file():
                    shutil.copyfile(sidecar,str(target)+suffix)
            return str(target)
        tracks=self.track_data()
        for track in tracks:
            track['path']=persistent(track['path'])
            if track.get('index'):
                track['index']=persistent(track['index'])
        data=dict(type='poreckit-sequence',version=1,reference=persistent(self.reference),
                  chromsizes=list(self.chromsizes.items()),region=self.region.text() if self.region else '',tracks=tracks)
        destination.restore(data,Path.cwd())
        new.tabs.setCurrentWidget(destination)
        return new

    def use_map(self):
        if self.browser.region:
            self.navigate(self.browser.region)

    def show_map(self):
        if self.region and self.browser.input_path:
            self.browser.locus.setText(self.region.text())
            self.browser.go()
            self.browser.tabs.setCurrentIndex(0)

    def map_changed(self,region):
        if self.link.isChecked() and self.chromsizes:
            self.navigate(region)

    def find_gene(self):
        from PySide6.QtWidgets import QInputDialog
        query,ok=QInputDialog.getText(self,'Find gene','Gene name, ID or alias in loaded GTF/GFF/BED tracks')
        if ok and query.strip():
            paths=[t['path'] for t in self.track_data() if t['kind'] in {'genes','intervals'}]
            self.request(dict(action='gene_search',paths=paths,query=query,chromsizes=list(self.chromsizes.items())),self.gene_results)

    def gene_results(self,result):
        from PySide6.QtWidgets import QInputDialog
        matches=result.get('results',[])
        if not matches:
            self.status.setText('No matching genes in loaded annotations')
            return
        labels=[f"{v['name']} · {v['chrom']}:{v['start']}-{v['end']}" for v in matches]
        label,ok=QInputDialog.getItem(self,'Find gene','Matches',labels,0,False)
        if ok:
            v=matches[labels.index(label)]
            r=Region(v['chrom'],v['start'],v['end'])
            self.navigate(r.window(r.start-1000,r.end-r.start+2000,self.chromsizes.get(r.chrom)))

    def build_index(self):
        item=self.tracks.currentItem()
        if not item:
            return
        self.rebuild_index(item.data(Qt.ItemDataRole.UserRole), lambda result:self.index_ready(item,result))

    def rebuild_index(self,track,done=None):
        """Generate a transient .bai for a local BAM track in this tab's
        temporary directory and invoke `done(result)` on completion."""
        from .sequenceio import file_kind
        if is_remote(track['path']):
            self.status.setText('Remote BAM needs a published BAI/CSI index; edit Index URL in Settings')
            return
        if file_kind(track['path'])!='alignment' or not source_path(track['path']).lower().endswith('.bam'):
            self.status.setText('Select a BAM track; CRAM uses its source CRAI index')
            return
        import hashlib

        from .sequenceio import version
        name=hashlib.sha256(repr(version(track['path'])).encode()).hexdigest()+'.bai'
        output=str(Path(self.temp.name)/name)
        self.request(dict(action='sequence_index',path=track['path'],output=output),
                     done if done else lambda result:self.index_ready(self.tracks.currentItem(),result))

    def index_ready(self,item,result):
        if self.tracks.row(item)<0:
            return
        track=item.data(Qt.ItemDataRole.UserRole)
        track['index']=result['index']
        self.loading=True
        item.setData(Qt.ItemDataRole.UserRole,track)
        self.loading=False
        self.query()

    def track_index_ready(self,item,result):
        track=item.data(Qt.ItemDataRole.UserRole)
        track['index']=result['index']
        self.loading=True
        item.setData(Qt.ItemDataRole.UserRole,track)
        self.loading=False
        self.query()

    def session_data(self,base):
        from .sessionio import to_relative
        def relative(path):
            return to_relative(path, base)
        tracks=self.track_data()
        for path in [self.reference,*[t['path'] for t in tracks]]:
            resolved=Path(path).resolve() if path else path
            if path and (resolved.is_relative_to(Path(self.browser.temp.name).resolve()) or resolved.is_relative_to(Path(self.temp.name).resolve())):
                raise ValueError('Synthetic demo files are temporary. Load your own data before saving a session.')
        for track in tracks:
            track['path']=relative(track['path'])
            # Transient generated indexes are rebuilt after reopening.
            if track.get('index'):
                if is_remote(track['index']):
                    continue
                index=Path(track['index'])
                if index.resolve().is_relative_to(Path(self.temp.name).resolve()):
                    track.pop('index',None)
                    track['index_rebuild']=True
                else:
                    track['index']=relative(str(index))
        return dict(type='poreckit-sequence',version=1,reference=relative(self.reference),
                    chromsizes=list(self.chromsizes.items()),region=self.region.text() if self.region else '',tracks=tracks)

    def save_session(self):
        path,_=QFileDialog.getSaveFileName(self,'Save track session','tracks.json','JSON (*.json)')
        if path:
            if not Path(path).suffix:
                path+='.json'
            try:
                from .sessionio import save_json
                save_json(path, self.session_data(Path(path).parent))
                self.current_session=path
                self.status.setText('Session saved: '+path)
            except (OSError,ValueError) as exc:
                self.status.setText(str(exc))

    def restore(self,data,base):
        if data.get('type')!='poreckit-sequence' or data.get('version')!=1:
            raise ValueError('Unsupported track session')
        from .sessionio import to_absolute
        def absolute(path):
            p=to_absolute(path,base)
            if not is_remote(path) and not Path(p).is_file():
                raise ValueError('File not found: '+str(p))
            return p
        import copy
        tracks=copy.deepcopy(data['tracks'])
        reference=absolute(data['reference']) if data.get('reference') else ''
        for track in tracks:
            track['path']=absolute(track['path'])
            if track.get('index'):
                track['index']=absolute(track['index'])
            if file_kind(track['path'])!=track['kind']:
                raise ValueError('Track type does not match file format')
        sizes=dict(data['chromsizes'])
        if any(not isinstance(n,int) or n<=0 for n in sizes.values()):
            raise ValueError('Invalid chromosome lengths')
        region=Region.parse(data['region']) if data.get('region') else None
        if region and (region.chrom not in sizes or region.end>sizes[region.chrom]):
            raise ValueError('Session locus exceeds chromosome lengths')
        self.cancel()
        self.loading=True
        self.tracks.clear()
        self.loading=False
        self.reference,self.chromsizes=reference,sizes
        self.history=[]
        self.history_index=-1
        for track in tracks:
            self.append(track)
        for row,track in enumerate(tracks):
            # Transient generated indexes are not saved; rebuild them before
            # the first query needs them. The item keeps its own dict copy.
            if track.pop('index_rebuild',None):
                item=self.tracks.item(row)
                self.rebuild_index(item.data(Qt.ItemDataRole.UserRole),
                                   lambda result,item=item:self.track_index_ready(item,result))
        self.region=region
        self.sync_chroms()
        if region:
            self.navigate(region)

    def load_session(self):
        path,_=QFileDialog.getOpenFileName(self,'Open track session','','JSON (*.json)')
        if path:
            try:
                from .sessionio import load_json
                self.restore(load_json(path),Path(path).parent)
                self.current_session=path
            except (OSError,ValueError,KeyError,TypeError) as exc:
                self.status.setText(str(exc))

    def dragEnterEvent(self,event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self,event):
        paths=[u.toLocalFile() if u.isLocalFile() else u.toString() for u in event.mimeData().urls() if u.isLocalFile() or is_remote(u.toString())]
        if paths:
            self.add_files(paths)
            event.acceptProposedAction()
