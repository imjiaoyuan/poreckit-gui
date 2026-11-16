"""Native scaffold editing with a live, lifted contact map."""
import json
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,QCheckBox,QComboBox,QFileDialog,QHBoxLayout,QInputDialog,
    QLabel,QLineEdit,QMessageBox,QPushButton,QSpinBox,QSplitter,QTreeWidget,
    QTreeWidgetItem,QVBoxLayout,
)

from .assembly import Assembly
from .dataview import DataView
from .model import Region
from .widgets import ImageView


class ScaffoldTree(QTreeWidget):
    def __init__(self,editor):
        super().__init__()
        self.editor=editor
        self.setColumnCount(4)
        self.setHeaderLabels(['Scaffold / fragment','Direction','Length','Source'])
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.setDefaultDropAction(Qt.DropAction.MoveAction)

    def dropEvent(self,event):
        super().dropEvent(event)
        if self.editor.document:
            groups=[]
            for i in range(self.topLevelItemCount()):
                root=self.topLevelItem(i)
                groups.append([root.child(j).data(0,Qt.ItemDataRole.UserRole) for j in range(root.childCount())])
            try:
                self.editor.document.edit(lambda:setattr(self.editor.document,'groups',groups))
                self.editor.refresh()
            except (ValueError,TypeError):
                self.editor.refresh()


class AssemblyEditor(DataView):
    def __init__(self,browser):
        super().__init__(browser)
        self.document=None
        self.input_path=''
        self.reference=''
        self.session_path=None
        self.preview_result=None
        self.preview_document=None
        layout=QVBoxLayout(self)
        row=QHBoxLayout()
        for label,callback in [('Use loaded map',self.from_map),('Import .assembly…',self.import_assembly),('Open project…',self.load_session),('Save project…',self.save_session),('Export…',self.export_menu),('Cancel',self.cancel)]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            row.addWidget(b)
        row.addStretch()
        layout.addLayout(row)
        self.source=QLabel('Load a contact map in Genome Browser, then initialize or import an assembly.')
        self.source.setWordWrap(True)
        layout.addWidget(self.source)
        row=QHBoxLayout()
        for label,callback in [('Invert',self.invert),('Move up',lambda:self.move(-1)),('Move down',lambda:self.move(1)),('Move to scaffold…',self.move_to),('Cut fragment…',self.cut),('Split scaffold',self.split),('Join next scaffold',self.join),('Undo',self.undo),('Redo',self.redo)]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            row.addWidget(b)
        layout.addLayout(row)
        row=QHBoxLayout()
        self.locus=QLineEdit()
        self.locus.setPlaceholderText('assembly:0-1,000,000 · zero-based, half-open')
        self.locus.returnPressed.connect(self.preview)
        self.resolution=QComboBox()
        self.resolution.currentIndexChanged.connect(lambda _:self.preview())
        self.log_scale=QCheckBox('Log colors')
        self.log_scale.setChecked(True)
        self.log_scale.toggled.connect(lambda _:self.preview())
        row.addWidget(QLabel('Preview locus'))
        row.addWidget(self.locus,1)
        row.addWidget(QLabel('Assembly bp / bin'))
        row.addWidget(self.resolution)
        row.addWidget(self.log_scale)
        for label,callback in [('Update map',self.preview),('Whole assembly',self.whole),('−',lambda:self.zoom(2)),('+',lambda:self.zoom(.5))]:
            b=QPushButton(label)
            b.clicked.connect(callback)
            row.addWidget(b)
        layout.addLayout(row)
        split=QSplitter()
        self.tree=ScaffoldTree(self)
        self.tree.itemDoubleClicked.connect(self.focus)
        split.addWidget(self.tree)
        self.canvas=ImageView()
        self.canvas.data_navigation=True
        self.canvas.dataZoom.connect(lambda factor,x,y:self.zoom(factor,x))
        self.canvas.dataPan.connect(lambda dx,dy:self.pan(dx))
        self.canvas.binHovered.connect(self.hover)
        split.addWidget(self.canvas)
        split.setSizes([370,850])
        split.setStretchFactor(1,1)
        layout.addWidget(split,1)
        self.hover_label=QLabel('Drag fragments between scaffolds; double-click a fragment to focus the map.')
        self.hover_label.setWordWrap(True)
        layout.addWidget(self.hover_label)
        self.status=QLabel('Ready')
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

    def allow_replace(self):
        if self.document and self.document.dirty:
            answer=QMessageBox.question(self,'Unsaved assembly','Discard unsaved assembly edits?',QMessageBox.StandardButton.Discard|QMessageBox.StandardButton.Cancel)
            return answer==QMessageBox.StandardButton.Discard
        return True

    def from_map(self):
        b=self.browser
        if not b.input_path or not b.metadata.get('chromsizes'):
            self.status.setText('Load a contact map with chromosome lengths first')
            return
        if not self.allow_replace():
            return
        self.input_path=b.input_path
        self.document=Assembly.from_chromsizes(b.metadata['chromsizes'])
        self.configure(b.metadata['resolutions'])

    def import_assembly(self):
        b=self.browser
        if not b.input_path:
            self.status.setText('Load the original contact map first')
            return
        if not self.allow_replace():
            return
        path,_=QFileDialog.getOpenFileName(self,'Import .assembly','','Assembly (*.assembly)')
        if not path:
            return
        mode,ok=QInputDialog.getItem(self,'Source coordinates','How does the original map store the assembly?',
                                  ['Concatenated assembly chromosome','One chromosome per original fragment'],0,False)
        if not ok:
            return
        chrom='assembly'
        if mode.startswith('Concatenated'):
            names=[name for name,_ in b.metadata.get('chromsizes',[])]
            chrom,ok=QInputDialog.getItem(self,'Source chromosome','Original concatenated chromosome',names,0,False)
            if not ok:
                return
        scale,ok=QInputDialog.getInt(self,'Coordinate scale','Original assembly bp per contact-map bp (the common default is 1 + total bp // 2,100,000,000)',1,1,1000000)
        if not ok:
            return
        try:
            document=Assembly.read(path,chrom,scale,by_name=mode.startswith('One chromosome'))
            self.input_path=b.input_path
            self.document=document
            self.configure(b.metadata['resolutions'])
        except (OSError,ValueError) as exc:
            self.status.setText(str(exc))

    def configure(self,levels):
        scale=next(iter(self.document.fragments.values()))['scale']
        choices=sorted(levels)
        if not choices:
            self.status.setText('Assembly editing needs a contact map with stored resolutions')
            return
        self.resolution.blockSignals(True)
        self.resolution.clear()
        if str(self.input_path).lower().endswith((".cool", ".mcool")):
            import math
            base=choices[-1]
            target=max(1,math.ceil(self.document.length/scale/600))
            choices=sorted(set([*choices,math.ceil(target/base)*base]))
        self.resolution.addItems([str(level*scale) for level in choices])
        target=max(1,self.document.length//600)
        level=next((v*scale for v in choices if v*scale>=target),choices[-1]*scale)
        self.resolution.setCurrentText(str(level))
        self.resolution.blockSignals(False)
        self.locus.setText(Region('assembly',0,self.document.length).text())
        self.refresh()

    def refresh(self):
        if not self.document:
            return
        selected={abs(item.data(0,Qt.ItemDataRole.UserRole)) for item in self.tree.selectedItems() if item.parent() is not None}
        self.tree.clear()
        for i,group in enumerate(self.document.groups):
            length=sum(self.document.fragments[abs(v)]['length'] for v in group)
            root=QTreeWidgetItem([f'scaffold_{i+1}','',f'{length:,}',''])
            root.setFlags(root.flags() & ~Qt.ItemFlag.ItemIsDragEnabled)
            self.tree.addTopLevelItem(root)
            for signed in group:
                f=self.document.fragments[abs(signed)]
                item=QTreeWidgetItem([f['name'],'−' if signed<0 else '+',f"{f['length']:,}",f"{f['source_name']}:{f['source_start']:,}"])
                item.setData(0,Qt.ItemDataRole.UserRole,signed)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsDropEnabled)
                root.addChild(item)
                item.setSelected(abs(signed) in selected)
            root.setExpanded(True)
        self.source.setText(f"{Path(self.input_path).name} · {len(self.document.groups)} scaffolds · {len(self.document.fragments)} fragments · {self.document.length:,} bp"+(' · unsaved edits' if self.document.dirty else ''))
        self.preview()

    def selection(self):
        selected=self.tree.selectedItems()
        if not selected:
            raise ValueError('Select one or more contiguous fragments')
        if any(item.parent() is None for item in selected):
            raise ValueError('Select fragment rows inside a scaffold')
        parents={self.tree.indexOfTopLevelItem(item.parent()) for item in selected}
        if len(parents)!=1:
            raise ValueError('Select fragments within one scaffold')
        group=parents.pop()
        positions=sorted(item.parent().indexOfChild(item) for item in selected)
        if positions!=list(range(positions[0],positions[-1]+1)):
            raise ValueError('Select contiguous fragments')
        return group,positions[0],positions[-1]

    def perform(self,operation):
        if not self.document:
            return
        try:
            operation()
            self.refresh()
        except (ValueError,IndexError) as exc:
            self.status.setText(str(exc))

    def invert(self):
        self.perform(lambda:self.document.invert(*self.selection()))

    def move(self,delta):
        def operation():
            group,first,last=self.selection()
            target=max(0,min(len(self.document.groups[group]),first-1 if delta<0 else last+2))
            self.document.move(group,first,last,group,target)
        self.perform(operation)

    def move_to(self):
        def operation():
            group,first,last=self.selection()
            dest,ok=QInputDialog.getInt(self,'Move fragments','Target scaffold number (append at end)',1,1,len(self.document.groups))
            if ok:
                self.document.move(group,first,last,dest-1,len(self.document.groups[dest-1]))
        self.perform(operation)

    def cut(self):
        def operation():
            group,first,last=self.selection()
            if first!=last:
                raise ValueError('Select one fragment to cut')
            f=self.document.fragments[abs(self.document.groups[group][first])]
            # Text input supports coordinates above Qt's signed 32-bit spinbox range.
            value,ok=QInputDialog.getText(self,'Cut fragment','Cut offset from the left in the current orientation (bp)')
            if ok:
                self.document.split_fragment(group,first,int(value.replace(',','')))
        self.perform(operation)

    def split(self):
        self.perform(lambda:self.document.split_group(self.selection()[0],self.selection()[1]))

    def join(self):
        self.perform(lambda:self.document.join(self.selection()[0]))

    def undo(self):
        self.perform(self.document.undo if self.document else lambda:None)

    def redo(self):
        self.perform(self.document.redo if self.document else lambda:None)

    def focus(self,item,*_):
        if not self.document or item.parent() is None:
            return
        signed=item.data(0,Qt.ItemDataRole.UserRole)
        fragment=next(f for f in self.document.layout() if f['id']==abs(signed))
        span=fragment['length']*2
        self.locus.setText(Region('assembly',0,self.document.length).window(fragment['start']-span/4,span,self.document.length).text())
        self.preview()

    def whole(self):
        if self.document:
            self.locus.setText(Region('assembly',0,self.document.length).text())
            self.preview()

    def zoom(self,factor,fraction=.5):
        if self.document:
            try:
                r=Region.parse(self.locus.text())
                span=r.end-r.start
                r=r.window(r.start+span*fraction-span*factor*fraction,span*factor,self.document.length)
                self.locus.setText(r.text())
                self.preview()
            except ValueError as exc:
                self.status.setText(str(exc))

    def pan(self,fraction):
        if self.document:
            r=Region.parse(self.locus.text()).pan(fraction,self.document.length)
            self.locus.setText(r.text())
            self.preview()

    def preview(self):
        if not self.document or not self.resolution.currentText():
            return
        try:
            region=Region.parse(self.locus.text())
            if region.chrom!='assembly' or region.end>self.document.length:
                raise ValueError('Preview locus must lie within the edited assembly')
            job=dict(action='assembly_preview',assembly=self.document.data(),input=self.input_path,
                     start=region.start,end=region.end,resolution=int(self.resolution.currentText()),
                     log=self.log_scale.isChecked(),output=str(Path(self.temp.name)/'assembly.svg'))
            self.request(job,self.ready)
        except ValueError as exc:
            self.status.setText(str(exc))

    def ready(self,result):
        self.preview_result=result
        self.preview_document=Assembly.restore(result['assembly'])
        self.canvas.show_image(result['output'])
        self.canvas.set_map_data(result['map_data'])
        self.status.setText('Preview updated · cuts within source bins use fractional overlap; source data is preserved')

    def hover(self,record):
        if record and self.preview_document:
            try:
                a=self.preview_document.lift(record['x'][1])
                b=self.preview_document.lift(record['y'][1])
                self.hover_label.setText(f"X {a['fragment']} → {a['chrom']}:{a['position']:,} · Y {b['fragment']} → {b['chrom']}:{b['position']:,} · contacts {record['value']:g}")
            except ValueError:
                pass

    def export_menu(self):
        if not self.document:
            return
        kind,ok=QInputDialog.getItem(self,'Export assembly','Format',
                                   ['.assembly','AGP','Corrected FASTA','Corrected contact map (.cool)','Map image'],0,False)
        if not ok:
            return
        if kind=='Map image':
            self.export()
            return
        suffix={'.assembly':'.assembly','AGP':'.agp','Corrected FASTA':'.fa','Corrected contact map (.cool)':'.cool'}[kind]
        path,_=QFileDialog.getSaveFileName(self,'Export assembly','corrected'+suffix,'Output (*'+suffix+')')
        if not path:
            return
        if not Path(path).suffix:
            path+=suffix
        try:
            if Path(path).resolve()==Path(self.input_path).resolve():
                raise ValueError('Choose a new file; the source contact map must be preserved')
            if kind=='.assembly':
                self.document.write(path)
                self.status.setText('Exported '+path+' · Save project also retains source coordinate mapping')
            elif kind=='AGP':
                self.document.write_agp(path)
                self.status.setText('Exported '+path)
            else:
                job=dict(action='assembly_export',assembly=self.document.data(),input=self.input_path,output=path,
                         format='fasta' if kind=='Corrected FASTA' else 'cool',resolution=int(self.resolution.currentText()))
                if kind=='Corrected FASTA':
                    reference,_=QFileDialog.getOpenFileName(self,'Original contig FASTA','','FASTA (*.fa *.fasta *.fna)')
                    if not reference:
                        return
                    if Path(path).resolve()==Path(reference).resolve():
                        raise ValueError('Choose a destination different from the reference FASTA')
                    job['reference']=reference
                self.request(job,lambda result:self.status.setText('Exported '+result['output']))
        except (OSError,ValueError) as exc:
            self.status.setText(str(exc))

    def export(self):
        if not self.preview_result:
            return
        path,_=QFileDialog.getSaveFileName(self,'Export assembly map','assembly.png','PNG (*.png);;SVG (*.svg)')
        if path:
            import shutil
            if path.lower().endswith('.svg'):
                shutil.copyfile(self.preview_result['output'],path)
            elif not self.canvas.scene().sceneRect().isEmpty():
                from .dataview import export_scene
                export_scene(self.canvas.scene(),path)

    def save_session(self):
        if not self.document:
            return
        path,_=QFileDialog.getSaveFileName(self,'Save assembly project',self.session_path or 'assembly-project.json','JSON (*.json)')
        if path:
            if not Path(path).suffix:
                path+='.json'
            try:
                from .sessionio import save_json, to_relative
                source=Path(self.input_path).resolve() if self.input_path else None
                roots=(Path(self.browser.temp.name).resolve(),Path(self.temp.name).resolve())
                if source and any(source.is_relative_to(root) for root in roots):
                    raise ValueError('Synthetic demo files are temporary. Load your own data before saving a project.')
                data=dict(type='poreckit-assembly',version=1,input=to_relative(self.input_path,Path(path).parent),
                          assembly=self.document.data(),region=self.locus.text(),resolution=int(self.resolution.currentText()))
                save_json(path,data)
                self.session_path=path
                self.document.mark_saved()
                self.source.setText(self.source.text().replace(' · unsaved edits',''))
            except (OSError,ValueError) as exc:
                self.status.setText(str(exc))

    def load_session(self):
        if not self.allow_replace():
            return
        path,_=QFileDialog.getOpenFileName(self,'Open assembly project','','JSON (*.json)')
        if not path:
            return
        try:
            from .sessionio import load_json, to_absolute
            data=load_json(path)
            if data.get('type')!='poreckit-assembly' or data.get('version')!=1:
                raise ValueError('Unsupported assembly project')
            doc=Assembly.restore(data['assembly'])
            source=Path(to_absolute(data['input'],Path(path).parent))
            if not source.is_file():
                raise ValueError('Original contact map is missing')
            self.request(dict(action='inspect',input=str(source)),lambda meta:self.accept_project(meta,data,doc,str(source),path))
        except (OSError,ValueError,KeyError,TypeError) as exc:
            self.status.setText(str(exc))

    def restore_controls(self,meta,data):
        self.configure(meta['resolutions'])
        self.resolution.setCurrentText(str(data['resolution']))
        self.locus.setText(data['region'])
        self.preview()

    def accept_project(self,meta,data,doc,source,path):
        self.document,self.input_path=doc,source
        self.session_path=path
        self.restore_controls(meta,data)
