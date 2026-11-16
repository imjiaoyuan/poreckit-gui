"""Editable annotation-layer tools connected to the genome viewport."""
import json
from pathlib import Path

from PySide6.QtCore import QObject, Qt
from PySide6.QtGui import QColor, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QColorDialog,
    QDialog,
    QDialogButtonBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .annotations import AnnotationModel, read_bedpe, write_bedpe


class AnnotationController(QObject):
    def __init__(self, browser, file_menu, view_menu):
        super().__init__(browser)
        self.browser = browser
        self.model = AnnotationModel(self.refresh)
        self.snapshot_revision = -1
        self.snapshot_path = None
        self.ui_revision = 0
        self.dock = QDockWidget("2D annotation layers", browser)
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.layers = QListWidget()
        self.layers.setMaximumHeight(110)
        self.layers.currentRowChanged.connect(self.show_features)
        self.layers.itemChanged.connect(self.visibility_changed)
        layout.addWidget(self.layers)
        buttons = QHBoxLayout()
        for name, callback in [("New", self.new_layer), ("Load…", self.load), ("Export…", self.export_layer),
                               ("Color…", self.color), ("Remove", self.remove_layer)]:
            button = QPushButton(name)
            button.clicked.connect(callback)
            buttons.addWidget(button)
        layout.addLayout(buttons)
        row = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setPlaceholderText("Filter annotations…")
        self.search.textChanged.connect(self.show_features)
        row.addWidget(self.search)
        self.mirror = QCheckBox("Mirror across diagonal")
        self.mirror.toggled.connect(self.mirror_changed)
        row.addWidget(self.mirror)
        layout.addLayout(row)
        self.table = QTableWidget(0, 8)
        self.columns = ("xchrom", "x0", "x1", "ychrom", "y0", "y1", "name", "score")
        self.table.setHorizontalHeaderLabels(["X chr", "X start", "X end", "Y chr", "Y start", "Y end", "Name", "Score"])
        self.table.itemChanged.connect(self.cell_changed)
        self.table.itemDoubleClicked.connect(self.centre_feature)
        layout.addWidget(self.table)
        row = QHBoxLayout()
        for text, callback in [("Edit…", self.edit_selected_dialog), ("Delete selected", self.delete_selected),
                               ("Undo", self.model.undo), ("Redo", self.model.redo)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            row.addWidget(button)
        delete_key = QShortcut(QKeySequence.StandardKey.Delete, self.table)
        delete_key.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        delete_key.activated.connect(self.delete_key_pressed)
        self.status = QLabel("Draw on the matrix to add an annotation")
        row.addWidget(self.status, 1)
        layout.addLayout(row)
        self.dock.setWidget(panel)
        browser.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.dock)
        self.dock.hide()
        self.initial_size = False
        self.dock.visibilityChanged.connect(self.size_dock)
        view_menu.addAction(self.dock.toggleViewAction())
        browser.action(file_menu, "Load 2D annotations…", self.load)
        browser.action(view_menu, "Undo edit", lambda: self.active_history(False), "Ctrl+Z")
        browser.action(view_menu, "Redo edit", lambda: self.active_history(True), "Ctrl+Shift+Z")
        self.draw = browser.action(view_menu, "Draw 2D annotation", self.toggle_draw)
        self.draw.setCheckable(True)
        self.draw.setToolTip("Drag a rectangle on the matrix. Ctrl-click selects; Ctrl-drag moves; Ctrl+Shift-drag resizes.")
        draw_button = QToolButton()
        draw_button.setDefaultAction(self.draw)
        buttons.insertWidget(0, draw_button)
        browser.image.annotationSelected.connect(self.draw_rectangle)
        browser.image.featureSelected.connect(self.select_feature)
        browser.image.featureEdited.connect(self.edit_feature)

    def active_history(self, redo):
        editor = getattr(self.browser, "assembly_editor", None)
        if editor is not None and self.browser.tabs.currentWidget() is editor:
            (editor.redo if redo else editor.undo)()
        else:
            (self.model.redo if redo else self.model.undo)()

    def size_dock(self, visible):
        if visible and not self.initial_size:
            from PySide6.QtCore import QTimer
            self.initial_size = True
            QTimer.singleShot(0, lambda: self.browser.resizeDocks([self.dock], [280], Qt.Orientation.Vertical))

    def current(self):
        item = self.layers.currentItem()
        return next((v for v in self.model.layers if v["id"] == item.data(Qt.ItemDataRole.UserRole)), None) if item is not None else None

    def select(self, key):
        for row in range(self.layers.count()):
            if self.layers.item(row).data(Qt.ItemDataRole.UserRole) == key:
                self.layers.setCurrentRow(row)
                return

    def refresh(self):
        previous = self.current()
        selected = previous["id"] if previous else None
        self.layers.blockSignals(True)
        self.layers.clear()
        for layer in self.model.layers:
            item = QListWidgetItem(layer["name"])
            item.setData(Qt.ItemDataRole.UserRole, layer["id"])
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if layer.get("visible", True) else Qt.CheckState.Unchecked)
            self.layers.addItem(item)
        self.layers.blockSignals(False)
        if selected and any(v["id"] == selected for v in self.model.layers):
            self.select(selected)
        elif self.layers.count():
            self.layers.setCurrentRow(0)
        self.show_features()
        if self.ui_revision != self.model.revision:
            self.ui_revision = self.model.revision
            self.browser.mark_dirty()

    def show_features(self, *args):
        layer = self.current()
        needle = self.search.text().lower()
        values = [v for v in layer["features"] if not needle or needle in (v["id"] + " " + " ".join(str(v.get(k, "")) for k in self.columns)).lower()] if layer else []
        self.table.blockSignals(True)
        self.table.setRowCount(min(len(values), 1000))
        for row, value in enumerate(values[:1000]):
            for col, name in enumerate(self.columns):
                item = QTableWidgetItem(str(value.get(name, "")))
                item.setData(Qt.ItemDataRole.UserRole, value["id"])
                self.table.setItem(row, col, item)
        self.table.blockSignals(False)
        self.mirror.blockSignals(True)
        self.mirror.setChecked(bool(layer and layer.get("mirror", True)))
        self.mirror.blockSignals(False)
        self.status.setText(f"{len(values)} annotations" + (" · first 1000 shown" if len(values)>1000 else ""))

    def visibility_changed(self, item):
        self.model.set_layer(item.data(Qt.ItemDataRole.UserRole), visible=item.checkState() == Qt.CheckState.Checked)

    def mirror_changed(self, checked):
        if self.current():
            self.model.set_layer(self.current()["id"], mirror=checked)

    def cell_changed(self, item):
        layer = self.current()
        if layer:
            try:
                self.model.update(layer["id"], item.data(Qt.ItemDataRole.UserRole), **{self.columns[item.column()]: item.text()})
            except ValueError as exc:
                self.browser.error(exc)
                self.show_features()

    def new_layer(self):
        name, ok = QInputDialog.getText(self.browser, "New annotation layer", "Name")
        if ok:
            self.select(self.model.add_layer(name or "Annotations"))
            self.dock.show()

    def remove_layer(self):
        if self.current():
            self.model.remove_layer(self.current()["id"])

    def delete_selected(self):
        layer = self.current()
        ids = {item.data(Qt.ItemDataRole.UserRole) for item in self.table.selectedItems()}
        if layer:
            for key in ids:
                self.model.remove(layer["id"], key)

    def delete_key_pressed(self):
        if self.table.state() != QAbstractItemView.State.EditingState:
            self.delete_selected()

    def selected_feature(self):
        layer = self.current()
        items = self.table.selectedItems()
        if not layer or not items:
            return layer, None
        key = items[0].data(Qt.ItemDataRole.UserRole)
        return layer, next((v for v in layer["features"] if v["id"] == key), None)

    def edit_selected_dialog(self):
        layer, value = self.selected_feature()
        if layer is None or value is None:
            return
        dialog = FeatureDialog(self.browser, layer, value)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            self.model.update(layer["id"], value["id"], **dialog.values())
        except ValueError as exc:
            self.browser.error(exc)

    def color(self):
        if not self.current():
            return
        color = QColorDialog.getColor(parent=self.browser, options=QColorDialog.ColorDialogOption.DontUseNativeDialog)
        if color.isValid():
            self.model.set_layer(self.current()["id"], color=color.name())

    def load(self):
        path, _ = QFileDialog.getOpenFileName(self.browser, "Load 2D annotations", "", "BEDPE (*.bedpe *.txt *.tsv);;All files (*)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            try:
                key = self.model.add_layer(Path(path).stem, read_bedpe(path))
                self.model.mark_saved(key)
                self.select(key)
                self.dock.show()
            except (ValueError, OSError) as exc:
                self.browser.error(exc)

    def export_layer(self):
        layer = self.current()
        if not layer:
            return
        path, _ = QFileDialog.getSaveFileName(self.browser, "Export annotation layer", "annotations.bedpe", "BEDPE (*.bedpe)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            try:
                write_bedpe(path, layer)
                self.model.mark_saved(layer["id"])
            except OSError as exc:
                self.browser.error(exc)

    def toggle_draw(self, checked):
        if checked and self.browser.image.map_data is None:
            self.browser.error("Render a square map before drawing annotations")
            self.draw.setChecked(False)
            return
        self.browser.image.annotation_mode = checked
        if checked:
            self.dock.show()
            self.browser.statusBar().showMessage("Drag to add · Ctrl-click to select · Ctrl-drag to move · Ctrl+Shift-drag to resize")

    def select_feature(self, hit):
        if self.hit_value(hit) is None:
            return
        self.dock.show()
        self.select(hit["layer"])
        self.search.setText(hit["id"])
        self.table.selectRow(0)

    def hit_value(self, hit):
        layer = next((v for v in self.model.layers if v["id"] == hit["layer"]), None)
        return next((v for v in layer["features"] if v["id"] == hit["id"]), None) if layer else None

    def edit_feature(self, hit, dx, dy, resize):
        value = self.hit_value(hit)
        if value is None:
            return
        viewport = self.browser.current_viewport()
        dx, dy = round(dx*(viewport.x.end-viewport.x.start)), round(dy*(viewport.y.end-viewport.y.start))
        if hit.get("reflected"):
            dx, dy = dy, dx
        def length(chrom):
            return next((size for name, size in self.browser.metadata.get("chromsizes", [])
                         if name.removeprefix("chr") == chrom.removeprefix("chr")), None)
        lx, ly = length(value["xchrom"]), length(value["ychrom"])
        if resize:
            x1, y1 = value["x1"]+dx, value["y1"]+dy
            changes = dict(x1=max(value["x0"]+1, min(x1, lx) if lx else x1),
                           y1=max(value["y0"]+1, min(y1, ly) if ly else y1))
        else:
            if lx:
                dx = min(dx, lx-value["x1"])
            if ly:
                dy = min(dy, ly-value["y1"])
            dx, dy = max(-value["x0"], dx), max(-value["y0"], dy)
            changes = dict(x0=value["x0"]+dx, x1=value["x1"]+dx, y0=value["y0"]+dy, y1=value["y1"]+dy)
        self.model.update(hit["layer"], hit["id"], **changes)

    def draw_rectangle(self, x0, x1, y0, y1):
        viewport = self.browser.current_viewport().select(x0, x1, y0, y1)
        if self.current() is None:
            self.model.add_layer()
        layer = self.current()
        name, ok = QInputDialog.getText(self.browser, "New annotation", "Name (optional)")
        if ok:
            self.model.add(layer["id"], dict(xchrom=viewport.x.chrom, x0=viewport.x.start, x1=viewport.x.end,
                                            ychrom=viewport.y.chrom, y0=viewport.y.start, y1=viewport.y.end, name=name))

    def centre_feature(self, item):
        layer = self.current()
        if layer:
            value = next(v for v in layer["features"] if v["id"] == item.data(Qt.ItemDataRole.UserRole))
            from .model import Region
            from .viewport import Viewport
            self.browser.navigate_viewport(Viewport(Region(value["xchrom"], value["x0"], value["x1"]),
                                                     Region(value["ychrom"], value["y0"], value["y1"])))

    def view_options(self):
        if not self.model.layers:
            return {}
        if self.snapshot_revision != self.model.revision:
            path = Path(self.browser.temp.name) / f"layers-{self.model.revision}.json"
            path.write_text(json.dumps(self.model.layers, ensure_ascii=False), encoding="utf-8")
            self.snapshot_revision = self.model.revision
            self.snapshot_path = str(path)
        return {"layers_file": self.snapshot_path}


class FeatureDialog(QDialog):
    """Edit the coordinates and attributes of one loop/TAD annotation."""

    def __init__(self, parent, layer, value):
        super().__init__(parent)
        self.setWindowTitle("Edit annotation")
        self.feature = value
        form = QFormLayout(self)
        self.name = QLineEdit(str(value.get("name", "")))
        self.score = QLineEdit(str(value.get("score", "")))
        self.x0 = QLineEdit(str(value["x0"]))
        self.x1 = QLineEdit(str(value["x1"]))
        self.y0 = QLineEdit(str(value["y0"]))
        self.y1 = QLineEdit(str(value["y1"]))
        form.addRow("Name", self.name)
        form.addRow("Score", self.score)
        form.addRow(f"{value['xchrom']} start", self.x0)
        form.addRow(f"{value['xchrom']} end", self.x1)
        form.addRow(f"{value['ychrom']} start", self.y0)
        form.addRow(f"{value['ychrom']} end", self.y1)
        color_row = QHBoxLayout()
        label = QLabel(value.get("color") or layer["color"])
        def pick():
            chosen = QColorDialog.getColor(parent=self, initial=QColor(label.text()),
                                           options=QColorDialog.ColorDialogOption.DontUseNativeDialog)
            if chosen.isValid():
                label.setText(chosen.name())
        choose = QPushButton("Color…")
        choose.clicked.connect(pick)
        color_row.addWidget(choose)
        color_row.addWidget(label)
        color_row.addStretch(1)
        form.addRow("Color", color_row)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)
        self.color_label = label

    def values(self):
        parsed = QColor(self.color_label.text())
        return dict(name=self.name.text(), score=self.score.text(),
                    x0=self.x0.text(), x1=self.x1.text(), y0=self.y0.text(), y1=self.y1.text(),
                    color=parsed.name() if parsed.isValid() else "")
