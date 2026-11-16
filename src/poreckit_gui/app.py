"""Native browser with locus navigation, track lists and a shared view."""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from PySide6.QtCore import QLocale, QProcess, QProcessEnvironment, QSettings, Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QComboBox, QLineEdit, QCheckBox, QFormLayout,
    QListWidget, QListWidgetItem, QAbstractItemView, QSplitter, QToolBar,
    QFileDialog, QMessageBox, QDockWidget, QPlainTextEdit, QProgressBar,
    QGroupBox, QStackedWidget, QTabWidget, QDoubleSpinBox, QScrollArea, QSizePolicy,
)

from .model import Region, validate_view
from .widgets import ImageView, ChromosomeOverview
from .jobs import JobProcess
from .viewport import Viewport
from .contactview import DISPLAY_MODES, CONTROL_MODES, SIGNED_MODES

TRACK_NAMES = {"contacts": "Contact map", "genes": "Genes", "signal": "Signal", "tads": "TAD", "loops": "Loops"}


class Browser(QMainWindow):
    def __init__(self, settings=None):
        super().__init__()
        self.setWindowTitle("poreckit · Genome Browser")
        self.resize(1360, 900)
        self.setMinimumSize(1000, 660)
        self.settings = settings if settings is not None else QSettings("poreckit", "GenomeBrowser")
        self.temp = tempfile.TemporaryDirectory(prefix="poreckit-gui-")
        self.metadata = {}
        self.input_path = ""
        self.region = None
        self.control_path = ""
        self.control_metadata = {}
        self.pending_viewport = None
        self.pending_linked = False
        self.hovered_bin = None
        self.annotations = None
        self.history = []
        self.history_index = -1
        self.job = None
        self.serial = 0
        self.pending_session = None
        self.last_view = None
        self.last_files = {}
        self.cancelled = False
        self.loading_controls = True
        self.pending_display_update = False
        self.display_timer = QTimer(self)
        self.display_timer.setSingleShot(True)
        self.display_timer.setInterval(450)
        self.display_timer.timeout.connect(self.apply_display_changes)
        self.plot_workbench = None
        self.process = JobProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.read_log)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.process_error)
        self._build()
        self.loading_controls = False
        geometry = self.settings.value("geometry")
        if geometry:
            self.restoreGeometry(geometry)
        self.setAcceptDrops(True)

    def action(self, menu, title, callback, shortcut=None):
        action = QAction(title, self)
        action.triggered.connect(callback)
        if shortcut:
            action.setShortcut(shortcut)
        menu.addAction(action)
        return action

    def _build(self):
        file_menu = self.menuBar().addMenu("File")
        self.action(file_menu, "Load contact data…", self.open_file, QKeySequence.StandardKey.Open)
        self.action(file_menu, "Load control map…", self.open_control)
        self.action(file_menu, "Open genomic tracks…", self.open_sequence_files)
        self.action(file_menu, "Open remote track…", self.open_remote_track)
        self.action(file_menu, "Export matrix data…", self.export_matrix)
        annotation_menu = file_menu.addMenu("Add annotation track")
        self.annotation_menu = annotation_menu
        for kind, name in TRACK_NAMES.items():
            self.action(annotation_menu, name + "…", lambda checked=False, k=kind: self.add_track(k))
        file_menu.addSeparator()
        self.action(file_menu, "Open session…", self.load_session)
        self.action(file_menu, "Save session…", self.save_session, QKeySequence.StandardKey.Save)
        self.action(file_menu, "Export image…", self.export, "Ctrl+Shift+S")
        file_menu.addSeparator()
        self.action(file_menu, "Open synthetic demo", self.open_demo)
        self.action(file_menu, "Open synthetic track demo", self.open_track_demo)
        self.action(file_menu, "Quit", self.close, QKeySequence.StandardKey.Quit)

        toolbar = QToolBar("Data and navigation", self)
        self.genome_toolbar = toolbar
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        load = QPushButton("Load data…")
        load.clicked.connect(self.open_file)
        toolbar.addWidget(load)
        toolbar.addSeparator()
        self.control = QWidget()
        nav = QHBoxLayout(self.control)
        nav.setContentsMargins(4, 2, 4, 2)
        self.back = QPushButton("Back")
        self.back.setFixedWidth(60)
        self.back.setEnabled(False)
        self.back.clicked.connect(lambda: self.travel(-1))
        self.forward = QPushButton("Forward")
        self.forward.setFixedWidth(80)
        self.forward.setEnabled(False)
        self.forward.clicked.connect(lambda: self.travel(1))
        nav.addWidget(self.back)
        nav.addWidget(self.forward)
        self.chrom = QComboBox()
        self.chrom.setMinimumContentsLength(6)
        self.chrom.currentTextChanged.connect(self.chrom_changed)
        nav.addWidget(self.chrom)
        self.locus = QLineEdit("chr1:0-2,000,000")
        self.locus.setMinimumWidth(235)
        self.locus.setToolTip("Genomic region: zero-based, half-open")
        self.locus.returnPressed.connect(self.go)
        nav.addWidget(self.locus, 1)
        go = QPushButton("Go")
        go.setFixedWidth(50)
        go.clicked.connect(self.go)
        nav.addWidget(go)
        for text, callback in [("◀", lambda: self.pan(-0.5)), ("▶", lambda: self.pan(0.5)),
                               ("−", lambda: self.zoom(2)), ("+", lambda: self.zoom(0.5)),
                               ("Whole chr", self.home)]:
            button = QPushButton(text)
            button.setFixedWidth(80 if text == "Whole chr" else 32)
            button.setToolTip({"◀": "Pan left by half a window", "▶": "Pan right by half a window",
                               "−": "Zoom out", "+": "Zoom in"}.get(text, text))
            button.clicked.connect(callback)
            nav.addWidget(button)
        toolbar.addWidget(self.control)

        splitter = QSplitter()
        self.tabs = QTabWidget()
        self.tabs.addTab(splitter, "Genome Browser")
        self.setCentralWidget(self.tabs)
        self.sidebar = QWidget()
        side = QVBoxLayout(self.sidebar)
        self.dataset = QLabel("No data loaded")
        self.dataset.setWordWrap(True)
        self.dataset.setTextFormat(Qt.TextFormat.PlainText)
        side.addWidget(self.dataset)
        track_hint = QLabel("Tracks · Double-click to edit; drag to reorder")
        track_hint.setWordWrap(True)
        side.addWidget(track_hint)
        self.tracks = QListWidget()
        self.tracks.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.tracks.itemChanged.connect(lambda item: self.mark_dirty())
        self.tracks.itemDoubleClicked.connect(self.edit_track)
        self.tracks.model().rowsMoved.connect(self.mark_dirty)
        self.tracks.setMinimumHeight(150)
        side.addWidget(self.tracks, 1)
        row = QHBoxLayout()
        add = QPushButton("Add track…")
        add.clicked.connect(self.choose_track)
        remove = QPushButton("Remove")
        remove.clicked.connect(self.remove_track)
        edit = QPushButton("Settings…")
        edit.clicked.connect(self.edit_track)
        row.addWidget(add)
        row.addWidget(edit)
        row.addWidget(remove)
        side.addLayout(row)
        group = QGroupBox("Contact map display")
        form = QFormLayout(group)
        self.mode = QComboBox()
        self.mode.addItem("Triangle", "triangle")
        self.mode.addItem("Square", "square")
        self.mode.setCurrentIndex(1)
        form.addRow("View", self.mode)
        self.map_display = QComboBox()
        for value, label in DISPLAY_MODES.items():
            self.map_display.addItem(label, value)
        self.map_display.currentIndexChanged.connect(self.display_changed)
        form.addRow("Data", self.map_display)
        self.control_label = QLabel("No control map")
        self.control_label.setWordWrap(True)
        self.control_label.setTextFormat(Qt.TextFormat.PlainText)
        form.addRow("Control", self.control_label)
        self.control_norm = QComboBox()
        self.control_norm.addItems(["Same as observed", "raw", "balanced", "KR", "VC", "VC_SQRT"])
        form.addRow("Control norm", self.control_norm)
        self.pseudocount = QDoubleSpinBox()
        self.pseudocount.setRange(0, 1e9)
        self.pseudocount.setDecimals(4)
        self.pseudocount.setValue(1)
        form.addRow("Ratio pseudocount", self.pseudocount)
        self.resolution = QComboBox()
        self.resolution.setEditable(True)
        self.resolution.addItem("auto")
        self.resolution.setToolTip("Select auto or enter a resolution in bp. For mcool files, use a resolution stored in the file.")
        form.addRow("Resolution / bp", self.resolution)
        self.norm = QComboBox()
        self.norm.setEditable(True)
        self.norm.addItems(["raw", "balanced", "VC", "VC_SQRT", "KR"])
        form.addRow("Normalization", self.norm)
        self.tile_cache = QCheckBox("Reuse contact tiles when navigating")
        self.tile_cache.setChecked(True)
        self.tile_cache.setToolTip("Cache aligned 128-bin source tiles so overlapping regions reuse data reads.")
        form.addRow(self.tile_cache)
        self.cmap = QComboBox()
        self.cmap.addItems(["hic_whitered", "Reds", "inferno", "viridis", "magma", "YlOrRd", "RdBu_r"])
        self.cmap.setToolTip("White–red is suitable for observed contact counts. Diverging palettes are mainly used for differences or O/E maps.")
        form.addRow("Colormap", self.cmap)
        self.log = QCheckBox("log1p color scale")
        self.log.setChecked(False)
        form.addRow(self.log)
        self.color_range = QComboBox()
        self.color_range.addItem("Auto · 95th percentile", 95.0)
        self.color_range.addItem("Auto · 99th percentile", 99.0)
        self.color_range.addItem("Full data range", None)
        self.color_range.addItem("Manual maximum", "manual")
        self.color_range.setToolTip("Auto excludes the diagonal and saturates high contacts; it does not change the contact counts.")
        form.addRow("Color range", self.color_range)
        self.color_max = QDoubleSpinBox()
        self.color_max.setRange(0.000001, 1e12)
        self.color_max.setDecimals(6)
        self.color_max.setValue(50)
        self.color_max.setEnabled(False)
        self.color_max.setToolTip("Upper color limit in contact units, before any log1p transform")
        form.addRow("Color maximum", self.color_max)
        preset = QPushButton("Reset display settings")
        preset.setToolTip("Restore square view, white–red colors, linear scale and automatic color range. The file, region, resolution and normalization stay as selected.")
        preset.clicked.connect(self.reset_display_settings)
        form.addRow(preset)
        side.addWidget(group)
        self.auto_update = QCheckBox("Auto update")
        self.auto_update.setChecked(True)
        self.auto_update.toggled.connect(self.auto_update_changed)
        side.addWidget(self.auto_update)
        self.refresh = QPushButton("Update map")
        self.refresh.clicked.connect(self.render)
        side.addWidget(self.refresh)
        for combo in (self.mode, self.resolution, self.norm, self.cmap):
            combo.currentTextChanged.connect(self.mark_dirty)
        self.log.toggled.connect(self.mark_dirty)
        self.color_range.currentIndexChanged.connect(self.color_range_changed)
        self.color_max.valueChanged.connect(self.mark_dirty)
        self.control_norm.currentTextChanged.connect(self.mark_dirty)
        self.pseudocount.valueChanged.connect(self.mark_dirty)
        self.mode.currentTextChanged.connect(self.map_mode_changed)
        self.resolution.currentTextChanged.connect(self.update_normalizations)
        hint = QLabel("Coordinates: zero-based, half-open.\nDisplay settings apply automatically. Image zoom changes the viewport; top navigation changes the genomic region.")
        hint.setWordWrap(True)
        side.addWidget(hint)
        for combo in self.sidebar.findChildren(QComboBox):
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(8)
            combo.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
        self.map_display.currentTextChanged.connect(self.map_display.setToolTip)
        self.map_display.setToolTip(self.map_display.currentText())
        sidebar_scroll = QScrollArea()
        sidebar_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        sidebar_scroll.setWidgetResizable(True)
        sidebar_scroll.setWidget(self.sidebar)
        sidebar_scroll.setMinimumWidth(340)
        splitter.addWidget(sidebar_scroll)
        center = QWidget()
        content = QVBoxLayout(center)
        content.setContentsMargins(0, 0, 0, 0)
        axes_row = QHBoxLayout()
        self.independent_y = QCheckBox("Independent Y axis")
        self.independent_y.toggled.connect(self.y_axis_changed)
        axes_row.addWidget(self.independent_y)
        self.locus_y = QLineEdit("chr1:0-2,000,000")
        self.locus_y.setEnabled(False)
        self.locus_y.returnPressed.connect(self.go)
        axes_row.addWidget(self.locus_y, 1)
        swap = QPushButton("Swap X / Y")
        swap.clicked.connect(self.swap_axes)
        axes_row.addWidget(swap)
        content.addLayout(axes_row)
        self.overview = ChromosomeOverview()
        self.overview.navigate.connect(self.overview_go)
        content.addWidget(self.overview)
        self.view_label = QLabel("  Genome view")
        self.view_label.setTextFormat(Qt.TextFormat.PlainText)
        content.addWidget(self.view_label)
        self.stack = QStackedWidget()
        welcome = QWidget()
        welcome_layout = QVBoxLayout(welcome)
        welcome_layout.addStretch()
        title = QLabel("poreckit\nHi-C / Pore-C Genome Browser")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        font = title.font()
        font.setPointSize(22)
        title.setFont(font)
        welcome_layout.addWidget(title)
        subtitle = QLabel("Load .cool / .mcool / .hic / .h5 / .pairs / .parquet\nDrop a file to load it, or try the synthetic demo.")
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        welcome_layout.addWidget(subtitle)
        demo = QPushButton("Open synthetic demo")
        demo.clicked.connect(self.open_demo)
        welcome_layout.addWidget(demo, alignment=Qt.AlignmentFlag.AlignCenter)
        welcome_layout.addStretch()
        self.image = ImageView()
        self.image.dataZoom.connect(self.data_zoom)
        self.image.dataPan.connect(self.data_pan)
        self.image.dataSelection.connect(self.data_select)
        self.image.binHovered.connect(self.show_bin)
        image_tools = QHBoxLayout()
        for text, callback in [("Image zoom +", lambda: self.image.zoom_by(1.25)),
                               ("Image zoom −", lambda: self.image.zoom_by(0.8)),
                               ("Fit", self.image.fit), ("100%", self.image.actual_size),
                               ("Export…", self.export), ("Viewport…", self.export_viewport)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            image_tools.addWidget(button)
        image_tools.addStretch()
        self.navigate_data = QCheckBox("Navigate data")
        self.navigate_data.setChecked(True)
        self.navigate_data.setToolTip("Drag to query another region; wheel/double-click to query a zoomed region; Alt-drag selects a region. Ctrl-wheel zooms the image.")
        self.navigate_data.toggled.connect(lambda checked: setattr(self.image, "data_navigation", checked))
        image_tools.addWidget(self.navigate_data)
        content.addLayout(image_tools)
        self.stack.addWidget(welcome)
        self.stack.addWidget(self.image)
        content.addWidget(self.stack, 1)
        self.hover_label = QLabel("Hover over the map to inspect genomic bins and contact values")
        self.hover_label.setTextFormat(Qt.TextFormat.PlainText)
        self.hover_label.setWordWrap(True)
        content.addWidget(self.hover_label)
        splitter.addWidget(center)
        splitter.setSizes([290, 1070])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(1500)
        self.log_dock = QDockWidget("Run log", self)
        self.log_dock.setWidget(self.log_view)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.log_dock)
        self.log_dock.hide()
        view_menu = self.menuBar().addMenu("View")
        view_menu.addAction(self.log_dock.toggleViewAction())
        crosshair = self.action(view_menu, "Crosshair", lambda checked: setattr(self.image, "crosshair", checked))
        crosshair.setCheckable(True)
        self.action(view_menu, "Copy hovered bin", self.copy_bin, "Ctrl+Shift+C")
        grid = self.action(view_menu, "Gridlines", lambda checked: self.mark_dirty())
        grid.setCheckable(True)
        self.grid_action = grid
        self.action(view_menu, "Fit to window", lambda: self.active_image().fit(), "Ctrl+0")
        self.action(view_menu, "Image zoom in", lambda: self.active_image().zoom_by(1.25), "Ctrl+=")
        self.action(view_menu, "Image zoom out", lambda: self.active_image().zoom_by(0.8), "Ctrl+-")
        help_menu = self.menuBar().addMenu("Help")
        self.action(help_menu, "About and controls", self.about)
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(130)
        self.progress.setRange(0, 0)
        self.progress.hide()
        self.stop = QPushButton("Cancel")
        self.stop.clicked.connect(self.cancel)
        self.stop.hide()
        self.statusBar().addPermanentWidget(self.progress)
        self.statusBar().addPermanentWidget(self.stop)
        self.statusBar().showMessage("Ready · Load contact data to begin")
        from .workbench import PlotWorkbench
        self.plot_workbench = PlotWorkbench(self)
        self.tabs.addTab(self.plot_workbench, "Plot Workbench")
        self.tabs.currentChanged.connect(self.active_tab_changed)
        from .annotationui import AnnotationController
        self.annotations = AnnotationController(self, file_menu, view_menu)
        from .states import StateController
        self.bookmarks = StateController(self)
        from .genomeui import GenomeDialog
        self.genome_dialog = GenomeDialog(self)
        self.action(view_menu, "Whole-genome overview…", self.genome_dialog.open, "Ctrl+G")
        genome_button = QPushButton("Whole genome")
        genome_button.clicked.connect(self.genome_dialog.open)
        self.genome_toolbar.addWidget(genome_button)
        from .searchui import GeneSearchDialog
        self.gene_search = GeneSearchDialog(self)
        self.action(view_menu, "Find gene…", self.find_gene, "Ctrl+F")
        find_button = QPushButton("Find gene…")
        find_button.clicked.connect(self.gene_search.open_search)
        self.genome_toolbar.addWidget(find_button)
        self.locus.setToolTip("Enter a genomic region (zero-based, half-open), chromosome, or gene name")
        from .windows import WindowController
        self.windows = WindowController(self)
        from .sequenceview import SequenceBrowser
        self.sequence_browser = SequenceBrowser(self)
        self.tabs.addTab(self.sequence_browser, "Tracks / alignments")
        from .assemblyview import AssemblyEditor
        self.assembly_editor = AssemblyEditor(self)
        self.tabs.addTab(self.assembly_editor, "Assembly editor")

    def open_control(self):
        if self.busy() or not self.input_path:
            return
        path, _ = QFileDialog.getOpenFileName(self, "Load control map", "", "Contact maps (*.hic *.cool *.mcool *.h5);;All files (*)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            self.start_job({"action": "inspect", "input": str(Path(path).resolve()), "target": "control"}, "Reading control map…")

    def display_changed(self, *args):
        if self.loading_controls:
            return
        mode = self.map_display.currentData()
        if mode != "observed":
            self.mode.setCurrentIndex(self.mode.findData("square"))
            self.log.setChecked(False)
        if mode in SIGNED_MODES:
            self.cmap.setCurrentText("RdBu_r")
        self.mark_dirty()

    def map_mode_changed(self, *args):
        square = self.mode.currentData() == "square"
        self.independent_y.setEnabled(square)
        if not square:
            self.independent_y.setChecked(False)
            self.map_display.setCurrentIndex(0)
        if self.annotations:
            self.annotations.draw.setEnabled(square)
            if not square:
                self.annotations.draw.setChecked(False)

    def y_axis_changed(self, enabled):
        self.locus_y.setEnabled(enabled)
        if enabled and self.region and not self.loading_controls:
            self.locus_y.setText(self.locus.text())
        self.mark_dirty()

    def swap_axes(self):
        if self.busy():
            return
        viewport = self.current_viewport()
        self.navigate_viewport(Viewport(viewport.y, viewport.x))

    def current_viewport(self):
        if self.pending_viewport:
            return self.pending_viewport
        if self.last_view:
            return Viewport.from_view(self.last_view)
        region = Region.parse(self.locus.text())
        return Viewport(region, Region.parse(self.locus_y.text()) if self.independent_y.isChecked() else region)

    def navigate_viewport(self, viewport, remember=True, *, linked=False):
        self.pending_viewport = viewport
        self.pending_linked = linked
        if self.busy():
            return
        self.pending_viewport = None
        self.pending_linked = False
        self.independent_y.blockSignals(True)
        self.independent_y.setChecked(not viewport.diagonal)
        self.locus_y.setEnabled(not viewport.diagonal)
        self.independent_y.blockSignals(False)
        self.locus.setText(viewport.x.text())
        self.locus_y.setText(viewport.y.text())
        if not viewport.diagonal:
            self.mode.setCurrentIndex(self.mode.findData("square"))
        self.render(remember=remember, linked=linked)

    def data_zoom(self, factor, fx, fy):
        viewport = self.current_viewport()
        resolution = self.image.map_data.spec["binsize"] if self.image.map_data else 1
        if viewport.diagonal and abs(fx-fy) < resolution / (viewport.x.end-viewport.x.start):
            fy = fx
        self.navigate_viewport(viewport.zoom(factor, fx, fy, dict(self.metadata.get("chromsizes", []))))

    def data_pan(self, dx, dy):
        self.navigate_viewport(self.current_viewport().pan(dx, dy, dict(self.metadata.get("chromsizes", []))))

    def data_select(self, x0, x1, y0, y1):
        self.navigate_viewport(self.current_viewport().select(x0, x1, y0, y1, dict(self.metadata.get("chromsizes", []))))

    def show_bin(self, info):
        if not info:
            return
        self.hovered_bin = info
        def locus(value):
            return f"{value[0]}:{value[1]:,}-{value[2]:,}"
        text = f"X {locus(info['x'])} · Y {locus(info['y'])} · observed {info['observed']:g} · value {info['value']:g}"
        if "control" in info:
            text += f" · control {info['control']:g}"
        self.hover_label.setText(text)

    def copy_bin(self):
        if self.hovered_bin:
            QApplication.clipboard().setText(self.hover_label.text())

    def export_matrix(self):
        data = self.image.map_data
        if data is None:
            self.error("Render a square map before exporting matrix data")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export matrix and coordinate metadata", "matrix.npz", "NPZ (*.npz)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            import numpy as np
            arrays = {"observed": data.observed, "displayed": data.displayed,
                      "x": np.asarray(data.spec["x"], dtype=str), "y": np.asarray(data.spec["y"], dtype=str),
                      "binsize": np.asarray(data.spec["binsize"])}
            if data.control is not None:
                arrays["control"] = data.control
            np.savez(path if path.endswith(".npz") else path + ".npz", **arrays)

    def active_tab_changed(self, index):
        is_browser = self.tabs.currentIndex() == 0
        self.genome_toolbar.setVisible(is_browser)
        self.statusBar().setVisible(is_browser)

    def active_image(self):
        if self.tabs.currentWidget() is getattr(self, "assembly_editor", None):
            return self.assembly_editor.canvas
        if self.tabs.currentWidget() is getattr(self, "sequence_browser", None):
            return self.sequence_browser.canvas
        if self.tabs.currentWidget() is self.plot_workbench:
            return self.plot_workbench.image
        return self.image

    def find_gene(self):
        if self.tabs.currentWidget() is self.sequence_browser:
            self.sequence_browser.find_gene()
        else:
            self.gene_search.open_search()

    def error(self, message):
        QMessageBox.warning(self, "poreckit", str(message))

    def busy(self):
        return self.job is not None

    def mark_dirty(self, *args):
        if self.loading_controls or not self.last_view:
            return
        if self.auto_update.isChecked():
            self.pending_display_update = True
            self.statusBar().showMessage("Settings changed · Updating map…")
            self.display_timer.start()
        else:
            self.statusBar().showMessage("Settings changed · Click Update map to apply")

    def apply_display_changes(self):
        if self.busy():
            return
        if self.pending_display_update and self.auto_update.isChecked():
            self.render(remember=False)

    def auto_update_changed(self, enabled):
        if enabled:
            self.mark_dirty()
        else:
            self.display_timer.stop()
            self.pending_display_update = False

    def color_range_changed(self, *args):
        self.color_max.setEnabled(self.color_range.currentData() == "manual")
        self.mark_dirty()

    def reset_display_settings(self):
        self.mode.setCurrentIndex(self.mode.findData("square"))
        self.cmap.setCurrentText("hic_whitered")
        self.log.setChecked(False)
        self.color_range.setCurrentIndex(0)
        self.mark_dirty()

    def update_normalizations(self, *args):
        selected = self.norm.currentText()
        resolution = self.resolution.currentText()
        if resolution == "auto":
            available = sorted(self.metadata.get("resolutions", []))
            try:
                region = Region.parse(self.locus.text())
                span = region.end - region.start
            except ValueError:
                span = 2000000
            target = max(1, (span + 249) // 250)
            resolution = next((r for r in available if r >= target), available[-1] if available else target)
        try:
            resolution = int(resolution)
        except ValueError:
            return
        choices = ["raw"] + [name for name, levels in self.metadata.get("normalizations", {}).items() if resolution in levels]
        self.norm.blockSignals(True)
        self.norm.clear()
        self.norm.addItems(choices)
        self.norm.setEditable(False)
        self.norm.setCurrentText(selected if selected in choices else "raw")
        self.norm.blockSignals(False)

    def open_file(self):
        if self.busy():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Load contact data", self.settings.value("last_folder", ""),
            "Genomic data (*.cool *.mcool *.hic *.h5 *.pairs *.pairs.gz *.parquet *.bw *.bigwig *.bigWig *.bedgraph *.bg *.bam *.cram *.sam *.vcf *.vcf.gz *.bcf *.fa *.fasta *.bed *.gtf *.gff);;All files (*)",
            options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            self.load_file(path)

    def load_file(self, path):
        if self.busy():
            return
        from .sequenceio import file_kind
        try:
            file_kind(path)
        except ValueError:
            pass
        else:
            self.tabs.setCurrentWidget(self.sequence_browser)
            self.sequence_browser.add_files([path])
            return
        from .remote import is_remote
        if is_remote(path):
            self.error("Open a local contact map or a supported remote genomic track")
            return
        path = str(Path(path).resolve())
        self.tabs.setCurrentIndex(0)
        self.settings.setValue("last_folder", str(Path(path).parent))
        self.start_job({"action": "inspect", "input": path}, "Reading chromosomes and resolutions…")

    def open_sequence_files(self):
        self.tabs.setCurrentWidget(self.sequence_browser)
        self.sequence_browser.open_files()

    def open_remote_track(self):
        self.tabs.setCurrentWidget(self.sequence_browser)
        self.sequence_browser.open_url()

    def open_demo(self):
        if self.busy():
            return
        from .demo import create_demo
        self.load_file(create_demo(self.temp.name))

    def open_track_demo(self):
        from .demo import create_track_demo
        self.tabs.setCurrentWidget(self.sequence_browser)
        self.sequence_browser.reset()
        self.sequence_browser.add_files(create_track_demo(Path(self.temp.name) / "tracks"))

    def choose_track(self):
        self.annotation_menu.exec(self.cursor().pos())

    def add_track(self, kind):
        if self.tabs.currentWidget() is self.sequence_browser and kind != "contacts":
            self.sequence_browser.open_files()
            return
        if self.busy():
            return
        if not self.input_path:
            self.error("Load contact data before adding annotation tracks")
            return
        filters = {"contacts": "Contact maps (*.cool *.mcool *.hic *.h5 *.pairs *.pairs.gz *.parquet)",
                   "genes": "Genes (*.gtf *.gff *.gff3 *.gtf.gz *.gff3.gz)",
                   "signal": "Signal (*.bw *.bigwig *.bigWig *.bedgraph *.bg *.bedgraph.gz)",
                   "tads": "TAD (*.bed *.bed.gz)", "loops": "Loops (*.bedpe *.bedpe.gz)"}
        paths, _ = QFileDialog.getOpenFileNames(self, f"Add {TRACK_NAMES[kind]} track", "", filters[kind] + ";;All files (*)",
                                              options=QFileDialog.Option.DontUseNativeDialog)
        for path in paths:
            colors = {"tads": "#d4a500", "loops": "#299fca", "signal": "#397ca8", "genes": "#397ca8"}
            self.append_track({"kind": kind, "path": path, "visible": True, "placement": "auto",
                               "color": colors.get(kind, "#397ca8")})
        if paths:
            self.render()

    def append_track(self, track):
        item = QListWidgetItem(f"{TRACK_NAMES[track['kind']]} · {track.get('name') or Path(track['path']).name}")
        item.setData(Qt.ItemDataRole.UserRole, track)
        item.setToolTip(track["path"])
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        item.setCheckState(Qt.CheckState.Checked if track.get("visible", True) else Qt.CheckState.Unchecked)
        self.tracks.addItem(item)

    def edit_track(self, item=None):
        if self.busy():
            return
        item = item if isinstance(item, QListWidgetItem) else self.tracks.currentItem()
        if item is None:
            return
        from .tracksettings import TrackSettings
        dialog = TrackSettings(item.data(Qt.ItemDataRole.UserRole), self)
        if dialog.exec():
            track = dialog.track
            item.setData(Qt.ItemDataRole.UserRole, track)
            item.setText(f"{TRACK_NAMES[track['kind']]} · {track.get('name') or Path(track['path']).name}")
            self.mark_dirty()

    def remove_track(self):
        if not self.busy() and self.tracks.currentRow() >= 0:
            self.tracks.takeItem(self.tracks.currentRow())
            self.render()

    def track_data(self):
        tracks = []
        for index in range(self.tracks.count()):
            item = self.tracks.item(index)
            tracks.append({**item.data(Qt.ItemDataRole.UserRole),
                           "visible": item.checkState() == Qt.CheckState.Checked})
        return tracks

    def view_data(self):
        resolution = self.resolution.currentText().strip()
        region = Region.parse(self.locus.text())
        length = self.length(region.chrom)
        if self.metadata.get("chromsizes") and length is None:
            raise ValueError("Chromosome is not present in the current data file")
        if length is not None and region.end > length:
            raise ValueError(f"Region exceeds chromosome length: {length:,} bp")
        value = {"input": self.input_path, "region": region.text(),
                 "binsize": None if resolution == "auto" else int(resolution),
                 "mode": self.mode.currentData(), "norm": self.norm.currentText().strip(),
                 "cmap": self.cmap.currentText(), "log": self.log.isChecked(), "tracks": self.track_data()}
        if self.independent_y.isChecked():
            y = Region.parse(self.locus_y.text())
            y_length = self.length(y.chrom)
            if self.metadata.get("chromsizes") and y_length is None:
                raise ValueError("Y chromosome is not present in the data file")
            if y_length is not None and y.end > y_length:
                raise ValueError("Y region exceeds chromosome length")
            value["region_y"] = y.text()
        value.update(display_mode=self.map_display.currentData(), control=self.control_path,
                     pseudocount=self.pseudocount.value(), grid=self.grid_action.isChecked(),
                     tile_cache=self.tile_cache.isChecked())
        if self.annotations:
            value.update(self.annotations.view_options())
        if self.control_norm.currentText() != "Same as observed":
            value["control_norm"] = None if self.control_norm.currentText() == "raw" else self.control_norm.currentText()
        if value["display_mode"] in CONTROL_MODES and not self.control_path:
            raise ValueError("Load a control map before selecting this display")
        color_range = self.color_range.currentData()
        value["vmax"] = self.color_max.value() if color_range == "manual" else None
        value["vmax_percentile"] = None if color_range == "manual" else color_range
        if value["norm"] == "raw":
            value["norm"] = None
        validate_view(value)
        return value

    def length(self, chrom=None):
        return dict(self.metadata.get("chromsizes", [])).get(chrom or (self.region.chrom if self.region else ""))

    def go(self):
        if self.busy() or not self.input_path:
            return
        text = self.locus.text().strip()
        if ":" not in text:
            length = self.length(text)
            if length:
                self.resolution.setCurrentText("auto")
                self.navigate(Region(text, 0, length))
            else:
                self.gene_search.search_text(text)
                if self.region:
                    self.locus.setText(self.region.text())
            return
        self.render()

    def navigate(self, region, remember=True):
        if self.busy() or not self.input_path:
            return
        y = self.current_viewport().y if self.independent_y.isChecked() else region
        self.navigate_viewport(Viewport(region, y), remember)

    def pan(self, fraction):
        if self.region:
            self.navigate(self.region.pan(fraction, self.length()))

    def zoom(self, factor):
        if self.region:
            self.navigate_viewport(self.current_viewport().zoom(factor, lengths=dict(self.metadata.get("chromsizes", []))))

    def home(self):
        if self.region and self.length():
            self.navigate(Region(self.region.chrom, 0, self.length()))

    def overview_go(self, fraction):
        if self.region and self.length():
            span = self.region.end - self.region.start
            self.navigate(self.region.window(self.length() * fraction - span / 2, span, self.length()))

    def chrom_changed(self, chrom):
        length = self.length(chrom)
        if length and self.region and not self.busy():
            self.navigate(Region(chrom, 0, min(length, self.region.end - self.region.start)))

    def travel(self, delta):
        index = self.history_index + delta
        if not self.busy() and 0 <= index < len(self.history):
            self.history_index = index
            state = self.history[index]
            self.navigate_viewport(state if isinstance(state, Viewport) else Viewport(state, state), remember=False)

    def render(self, checked=False, *, remember=True, output=None, linked=False):
        if self.busy():
            return
        self.display_timer.stop()
        self.pending_display_update = False
        try:
            view = self.view_data()
        except (ValueError, KeyError) as exc:
            self.error(exc)
            return
        path = output or str(Path(self.temp.name) / f"preview-{self.serial + 1}.png")
        self.start_job({"action": "render", "view": view, "output": path,
                        "export": output is not None, "remember": remember, "linked": linked}, "Rendering contact map…")

    def start_job(self, job, message):
        self.serial += 1
        self.job = job
        self.cancelled = False
        self.result_path = Path(self.temp.name) / f"result-{self.serial}.json"
        job_path = Path(self.temp.name) / f"job-{self.serial}.json"
        job_path.write_text(json.dumps(job, ensure_ascii=False), encoding="utf-8")
        self.sidebar.setEnabled(False)
        self.control.setEnabled(False)
        self.overview.setEnabled(False)
        self.progress.show()
        self.stop.show()
        self.statusBar().showMessage(message)
        self.log_view.appendPlainText(f"\n[{job['action']}] {job.get('input') or job['view']['region']}")
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "1")
        env.insert("MPLBACKEND", "Agg")
        env.insert("MPLCONFIGDIR", str(Path(self.temp.name) / "matplotlib"))
        self.process.setProcessEnvironment(env)
        if getattr(sys, "frozen", False):
            args = ["--worker", str(job_path), str(self.result_path)]
        else:
            args = ["-m", "poreckit_gui.worker", str(job_path), str(self.result_path)]
        self.process.start(sys.executable, args)

    def read_log(self):
        text = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if text:
            self.log_view.appendPlainText(text.rstrip())

    def process_error(self, error):
        if error == QProcess.ProcessError.FailedToStart and self.job:
            self.log_view.appendPlainText(self.process.errorString())
            self.finished(-1, QProcess.ExitStatus.CrashExit)

    def finished(self, code, status):
        self.read_log()
        job, self.job = self.job, None
        self.sidebar.setEnabled(True)
        self.control.setEnabled(True)
        self.overview.setEnabled(True)
        self.progress.hide()
        self.stop.hide()
        if not job:
            return
        if self.cancelled:
            self.pending_session = None
            self.statusBar().showMessage("Task cancelled")
            return
        if code != 0 or not self.result_path.exists():
            self.pending_session = None
            self.log_dock.show()
            self.statusBar().showMessage("Task failed · See the run log")
            self.error("Data loading or rendering failed. See the run log for details.")
            return
        try:
            result = json.loads(self.result_path.read_text(encoding="utf-8"))
            if job["action"] == "inspect":
                if job.get("target") == "control":
                    self.control_path = result["input"]
                    self.control_metadata = result
                    self.control_label.setText(Path(self.control_path).name)
                    self.control_label.setToolTip(self.control_path)
                    self.mark_dirty()
                    self.statusBar().showMessage("Control map loaded · Select a data display")
                else:
                    self.loaded(result)
            elif job["action"] == "genome":
                self.genome_dialog.ready(result, job)
                self.statusBar().showMessage("Whole-genome overview ready")
            elif job["action"] == "gene_search":
                self.gene_search.ready(result)
                self.statusBar().showMessage(f'Gene search: {result["total"]} matches')
            elif job.get("export"):
                self.statusBar().showMessage(f"Exported: {result['output']}")
            else:
                self.last_files = result.get("files", {})
                self.image.show_image(self.last_files.get("svg") or result["output"])
                self.image.set_map_data(result.get("map_data"))
                self.stack.setCurrentWidget(self.image)
                self.last_view = job["view"]
                self.region = Region.parse(result["region"])
                viewport = Viewport(self.region, Region.parse(result.get("region_y") or result["region"]))
                self.last_view["region"] = viewport.x.text()
                if "region_y" in self.last_view or not viewport.diagonal:
                    self.last_view["region_y"] = viewport.y.text()
                self.locus.setText(self.region.text())
                self.locus_y.setText(viewport.y.text())
                self.chrom.blockSignals(True)
                self.chrom.setCurrentText(self.region.chrom)
                self.chrom.blockSignals(False)
                if job["remember"] and (not self.history or self.history[self.history_index] != viewport):
                    self.history = self.history[:self.history_index + 1] + [viewport]
                    self.history_index = len(self.history) - 1
                self.overview.update_region(self.region, self.length())
                self.back.setEnabled(self.history_index > 0)
                self.forward.setEnabled(self.history_index < len(self.history) - 1)
                display = result.get("display", {})
                settings = f"{DISPLAY_MODES.get(display.get('display_mode'), 'Observed')} · {display.get('mode', '')} · {display.get('norm', 'raw')} · {'log1p' if display.get('log') else 'linear'}"
                upper = display.get("color_max")
                if upper is not None:
                    settings += f" · color max {upper:g}"
                    if self.color_range.currentData() != "manual":
                        self.color_max.blockSignals(True)
                        self.color_max.setValue(upper)
                        self.color_max.blockSignals(False)
                self.view_label.setText(f"  {self.region.text()}   |   {result['binsize']:,} bp / bin   |   {settings}")
                if not viewport.diagonal:
                    self.view_label.setText(f"  X {viewport.x.text()} · Y {viewport.y.text()} | {result['binsize']:,} bp / bin | {settings}")
                self.statusBar().showMessage("Ready · Map updated")
                self.windows.completed(viewport, linked=job.get("linked", False) or self.pending_viewport is not None)
                self.sequence_browser.map_changed(self.region)
                self.reap_temporaries()
        except (ValueError, KeyError, OSError) as exc:
            self.error(exc)
        if self.pending_viewport:
            pending, self.pending_viewport = self.pending_viewport, None
            linked, self.pending_linked = self.pending_linked, False
            self.navigate_viewport(pending, linked=linked)
            return
        if self.pending_display_update and self.auto_update.isChecked():
            self.display_timer.start()

    def reap_temporaries(self):
        """Keep only the current and previous numbered preview artifacts;
        older successful runs and failed/cancelled leftovers are removed so
        navigating does not grow the temporary folder without bound."""
        import re
        keep = {self.serial, max(1, self.serial - 1)}
        for path in Path(self.temp.name).iterdir():
            match = re.fullmatch(r"(preview|job|result)-(\d+)(?:\..*)?", path.name)
            if match and int(match.group(2)) not in keep and path.is_file():
                try:
                    path.unlink()
                except OSError:
                    pass

    def loaded(self, result):
        self.display_timer.stop()
        self.pending_display_update = False
        self.loading_controls = True
        self.pending_viewport = None
        self.control_path = ""
        self.control_metadata = {}
        self.control_label.setText("No control map")
        self.independent_y.setChecked(False)
        self.map_display.setCurrentIndex(0)
        self.metadata = result
        self.input_path = result["input"]
        self.tracks.clear()
        self.stack.setCurrentIndex(0)
        self.image.clear_image()
        self.view_label.setText("  Data loaded · Awaiting rendering")
        self.history, self.history_index, self.last_view = [], -1, None
        self.last_files = {}
        self.dataset.setText(Path(self.input_path).name)
        self.setWindowTitle(f"poreckit · {Path(self.input_path).name}")
        self.dataset.setToolTip(self.input_path)
        self.chrom.blockSignals(True)
        self.chrom.clear()
        self.chrom.addItems([c for c, _ in result["chromsizes"]])
        self.chrom.setEditable(not bool(result["chromsizes"]))
        self.chrom.blockSignals(False)
        self.resolution.clear()
        self.resolution.addItems(["auto", *map(str, result["resolutions"])])
        self.resolution.setEditable(not bool(result["resolutions"]))
        self.update_normalizations()
        self.norm.setCurrentText("raw")
        if result["chromsizes"]:
            chrom, length = result["chromsizes"][0]
            # Single-resolution files need a window wide enough to show their bins.
            smallest = min(result["resolutions"], default=10000)
            region = Region(chrom, 0, min(length, max(2000000, smallest * 100)))
        else:
            region = Region("chr1", 0, 2000000)
        self.locus.setText(region.text())
        self.locus_y.setText(region.text())
        self.region = region
        has_session = self.pending_session is not None
        # The session was fully validated by validate_session() before the
        # file load; still guarantee the flag never stays stuck on a failure.
        try:
            if has_session:
                session, self.pending_session = self.pending_session, None
                self.gene_search.paths = list(session.get("gene_annotations", []))
                if "bookmarks" in session:
                    self.bookmarks.restore(session["bookmarks"])
                self.locus.setText(session["region"])
                if session.get("region_y"):
                    self.independent_y.setChecked(True)
                    self.locus_y.setText(session["region_y"])
                self.control_path = session.get("control", "")
                if self.annotations:
                    self.annotations.model.restore(session.get("layers", []))
                self.control_label.setText(Path(self.control_path).name if self.control_path else "No control map")
                self.map_display.setCurrentIndex(max(0, self.map_display.findData(session.get("display_mode", "observed"))))
                self.pseudocount.setValue(float(session.get("pseudocount", 1)))
                self.control_norm.setCurrentText(session.get("control_norm") or ("raw" if "control_norm" in session else "Same as observed"))
                self.grid_action.setChecked(session.get("grid", False))
                self.tile_cache.setChecked(session.get("tile_cache", True))
                self.mode.setCurrentIndex(1 if session["mode"] == "square" else 0)
                self.resolution.setCurrentText(str(session["binsize"]) if session["binsize"] else "auto")
                self.norm.setCurrentText(session.get("norm") or "raw")
                self.cmap.setCurrentText(session["cmap"])
                self.log.setChecked(session["log"])
                if session.get("vmax") is not None:
                    self.color_range.setCurrentIndex(self.color_range.findData("manual"))
                    self.color_max.setValue(float(session["vmax"]))
                else:
                    percentile = session.get("vmax_percentile", 95.0)
                    index = self.color_range.findData(percentile)
                    self.color_range.setCurrentIndex(index if index >= 0 else 0)
                for track in session["tracks"]:
                    self.append_track(track)
        finally:
            self.loading_controls = False
        if result["chromsizes"] or has_session:
            self.render()
        else:
            self.statusBar().showMessage("File loaded · No chromosome lengths available. Enter a locus and click Go.")

    def export(self):
        if self.tabs.currentWidget() is self.assembly_editor:
            self.assembly_editor.export()
            return
        if self.tabs.currentWidget() is self.sequence_browser:
            self.sequence_browser.export()
            return
        if self.tabs.currentWidget() is self.plot_workbench:
            self.plot_workbench.export()
            return
        if self.busy() or not self.last_view:
            return
        path, selected = QFileDialog.getSaveFileName(self, "Export full figure", "poreckit.png",
                                                     "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)",
                                                     options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            if not Path(path).suffix:
                path += {"PNG (*.png)": ".png", "PDF (*.pdf)": ".pdf", "SVG (*.svg)": ".svg"}[selected]
            if Path(path).suffix.lower() not in {".png", ".pdf", ".svg"}:
                self.error("Select PNG, PDF or SVG format")
                return
            extension = Path(path).suffix.lower().lstrip(".")
            if extension in self.last_files:
                try:
                    shutil.copyfile(self.last_files[extension], path)
                    self.statusBar().showMessage(f"Exported: {path}")
                except OSError as exc:
                    self.error(exc)
            else:
                self.start_job({"action": "render", "view": self.last_view, "output": path,
                                "export": True, "remember": False}, "Exporting figure…")

    def export_viewport(self):
        if self.tabs.currentWidget() is self.plot_workbench:
            self.plot_workbench.export_viewport()
            return
        if not self.last_view:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export visible viewport", "viewport.png", "PNG (*.png)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path and not self.image.viewport().grab().save(path if path.lower().endswith(".png") else path + ".png"):
            self.error("Unable to save the viewport")

    def session_data(self, include_bookmarks=False):
        import copy
        view = self.view_data()
        view["version"] = 1
        view.pop("layers_file", None)
        if self.annotations:
            view["layers"] = copy.deepcopy(self.annotations.model.layers)
        if include_bookmarks:
            view["bookmarks"] = copy.deepcopy(self.bookmarks.entries)
        view["gene_annotations"] = list(self.gene_search.paths)
        return view

    def validate_session(self, data):
        if data.get("version") != 1:
            raise ValueError("Unsupported session version")
        validate_view(data)
        # loaded() indexes these directly; reject early, before any state changes.
        for key in ("binsize", "cmap", "log", "tracks"):
            if key not in data:
                raise ValueError(f"Session is missing: {key}")
        if data.get("pseudocount", 1) is not None:
            try:
                if float(data.get("pseudocount", 1)) < 0:
                    raise ValueError("Pseudocount must be nonnegative")
            except (TypeError, ValueError) as exc:
                raise ValueError("Invalid pseudocount") from exc
        if data.get("vmax") is not None:
            try:
                if float(data["vmax"]) <= 0:
                    raise ValueError("Color maximum must be positive")
            except (TypeError, ValueError) as exc:
                raise ValueError("Invalid color maximum") from exc
        if data.get("control") and not Path(data["control"]).is_file():
            raise ValueError("Control map file does not exist")
        if not isinstance(data["log"], bool) or data["cmap"] not in [self.cmap.itemText(i) for i in range(self.cmap.count())]:
            raise ValueError("Invalid display settings")
        for track in data["tracks"]:
            if not isinstance(track, dict) or track.get("kind") not in TRACK_NAMES:
                raise ValueError("Unsupported track type")
            if not isinstance(track.get("path"), str) or not Path(track["path"]).is_file():
                raise ValueError(f"Track file not found: {track.get('path')}")
        for saved in data.get("bookmarks", []):
            self.validate_session(saved["view"])
        for annotation in data.get("gene_annotations", []):
            if not Path(annotation).is_file():
                raise ValueError(f"Gene annotation file not found: {annotation}")
        # Validate annotations before replacing any current state.
        from .annotations import AnnotationModel
        AnnotationModel().restore(data.get("layers", []))

    def apply_session(self, data):
        self.validate_session(data)
        self.pending_session = data
        if str(Path(data["input"]).resolve()) == self.input_path and self.metadata:
            self.tabs.setCurrentIndex(0)
            self.loaded(self.metadata)
        else:
            self.load_file(data["input"])

    def save_session(self):
        if self.tabs.currentWidget() is self.assembly_editor:
            self.assembly_editor.save_session()
            return
        if self.tabs.currentWidget() is self.sequence_browser:
            self.sequence_browser.save_session()
            return
        if self.tabs.currentWidget() is self.plot_workbench:
            self.plot_workbench.save_preset()
            return
        if self.busy() or not self.input_path:
            return
        try:
            view = self.session_data(include_bookmarks=True)
            if Path(view["input"]).resolve().is_relative_to(Path(self.temp.name).resolve()):
                raise ValueError("Synthetic demo data is temporary. Load your own data before saving a session.")
            path, _ = QFileDialog.getSaveFileName(self, "Save session (settings and file references; data is not copied)", "session.json", "JSON (*.json)",
                                                 options=QFileDialog.Option.DontUseNativeDialog)
            if not path:
                return
            if not Path(path).suffix:
                path += ".json"
            from .states import relative_references
            from .sessionio import save_json
            view = relative_references(view, Path(path).parent)
            save_json(path, view)
            if self.annotations:
                self.annotations.model.mark_saved()
            self.bookmarks.mark_saved()
            self.statusBar().showMessage(f"Session saved: {path}")
        except (ValueError, OSError) as exc:
            self.error(exc)

    def load_session(self):
        if self.tabs.currentWidget() is self.assembly_editor:
            self.assembly_editor.load_session()
            return
        if self.tabs.currentWidget() is self.sequence_browser:
            self.sequence_browser.load_session()
            return
        if self.tabs.currentWidget() is self.plot_workbench:
            self.plot_workbench.load_preset()
            return
        if self.busy():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open session", "", "JSON (*.json)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if not path:
            return
        try:
            from .states import absolute_references
            from .sessionio import load_json
            data = absolute_references(load_json(path), Path(path).parent)
            self.apply_session(data)
        except (ValueError, OSError, KeyError, TypeError, AttributeError) as exc:
            self.error(f"Unable to open session: {exc}")

    def cancel(self):
        self.display_timer.stop()
        self.pending_display_update = False
        self.pending_viewport = None
        if self.busy():
            self.cancelled = True
            self.process.kill()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls() and not self.busy():
            event.acceptProposedAction()

    def dropEvent(self, event):
        urls = event.mimeData().urls()
        if urls and urls[0].isLocalFile() and not self.busy():
            self.load_file(urls[0].toLocalFile())

    def about(self):
        QMessageBox.about(self, "About poreckit GUI",
            "poreckit Genome Browser\n\n"
            "Browse local contact maps with gene, signal, TAD and loop annotations.\n"
            "Coordinates are zero-based, half-open; + / − zoom, arrows pan.\n"
            "Export PNG / PDF / SVG; sessions save settings and data references.\n\n"
            "Synthetic demo data is for exploration and is deleted on exit.\n"
            "Plot Workbench exposes plotting commands, functions and all track types.\n"
            "Scroll to zoom any figure, drag to pan, and use Fit to reset.\n"
            "Export full figures or the visible viewport; save plot presets for reuse.")

    def closeEvent(self, event):
        if not self.assembly_editor.allow_replace():
            event.ignore()
            return
        self.display_timer.stop()
        if (self.annotations and self.annotations.model.dirty) or self.bookmarks.dirty:
            answer = QMessageBox.question(self, "Unsaved edits", "Save the session with annotations and bookmarks before closing?",
                QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel)
            if answer == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
            if answer == QMessageBox.StandardButton.Save:
                self.tabs.setCurrentIndex(0)
                self.save_session()
                if self.annotations.model.dirty or self.bookmarks.dirty:
                    event.ignore()
                    return
        if self.busy() or self.plot_workbench.busy() or self.sequence_browser.busy() or self.assembly_editor.busy():
            if QMessageBox.question(self, "Quit", "A task is running. Stop the task and quit?") != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.process.blockSignals(True)
            self.process.kill()
            self.process.waitForFinished(3000)
        self.plot_workbench.shutdown()
        self.sequence_browser.shutdown()
        self.assembly_editor.shutdown()
        self.bookmarks.timer.stop()
        self.process.shutdown()
        self.settings.setValue("geometry", self.saveGeometry())
        from .sessionio import cleanup_temp; cleanup_temp(self.temp.name)
        self.windows.close()
        event.accept()


def main():
    if sys.argv[1:2] == ["--self-test"]:
        if len(sys.argv) != 3:
            raise SystemExit("Usage: poreckit-gui --self-test REPORT.json")
        from .selftest import run
        raise SystemExit(run(sys.argv[2]))
    if sys.argv[1:2] == ["--worker"]:
        from .worker import main as worker_main, serve
        if sys.argv[2:] == ["--serve"]:
            serve()
            return
        raise SystemExit(worker_main(*sys.argv[2:]))
    QLocale.setDefault(QLocale("en_US"))
    app = QApplication(sys.argv)
    app.setApplicationName("poreckit")
    app.setOrganizationName("poreckit")
    window = Browser()
    window.show()
    raise SystemExit(app.exec())
