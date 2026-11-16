"""Whole-genome overview window and navigation into chromosome pairs."""
import shutil
from pathlib import Path

from PySide6.QtWidgets import QCheckBox, QComboBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout

from .model import Region
from .viewport import Viewport
from .widgets import ImageView


class GenomeDialog(QDialog):
    # Named to avoid clashing with QDialog.result(); holds the last
    # whole-genome job payload.
    overview: dict | None
    image: ImageView

    def __init__(self, browser):
        super().__init__(browser)
        self.browser = browser
        self.setWindowTitle("Whole-genome overview")
        self.resize(920, 850)
        self.setMinimumSize(650, 500)
        self.overview = None
        self.source = ""
        self.source_version = None
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel("Resolution / bp"))
        self.resolution = QComboBox()
        self.resolution.setEditable(True)
        self.resolution.addItem("auto")
        row.addWidget(self.resolution)
        self.log = QCheckBox("log1p color scale")
        self.log.setChecked(False)
        row.addWidget(self.log)
        update = QPushButton("Update overview")
        update.clicked.connect(self.update_overview)
        row.addWidget(update)
        for text, callback in [("Fit", lambda: self.image.fit()), ("+", lambda: self.image.zoom_by(1.25)),
                               ("−", lambda: self.image.zoom_by(.8)), ("Export…", self.export), ("Matrix…", self.export_matrix)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            row.addWidget(button)
        layout.addLayout(row)
        self.image = ImageView()
        self.image.data_navigation = False
        self.image.setToolTip("Scroll to zoom; drag to pan. Double-click a chromosome block to open it in the browser.")
        self.image.binActivated.connect(self.open_bin)
        self.image.binHovered.connect(self.hover)
        layout.addWidget(self.image, 1)
        self.label = QLabel("Raw contacts · Double-click a chromosome block to open it")
        layout.addWidget(self.label)

    def open(self):
        if not self.browser.input_path or self.browser.busy():
            return
        self.show()
        self.raise_()
        stat = Path(self.browser.input_path).stat()
        version = stat.st_mtime_ns, stat.st_size, stat.st_ino
        if self.source != self.browser.input_path or self.source_version != version or self.overview is None:
            self.overview = None
            self.image.clear_image()
            self.resolution.clear()
            self.resolution.addItems(["auto", *map(str, self.browser.metadata.get("resolutions", []))])
            self.update_overview()

    def update_overview(self):
        b = self.browser
        if b.busy() or not b.input_path:
            return
        try:
            text = self.resolution.currentText().strip()
            resolution = None if text == "auto" else int(text)
            if resolution is not None and resolution <= 0:
                raise ValueError("Overview resolution must be a positive integer")
            output = str(Path(b.temp.name) / f"genome-{b.serial+1}.svg")
            b.start_job(dict(action="genome", input=b.input_path, metadata=b.metadata,
                             resolution=resolution, output=output, log=self.log.isChecked()), "Building whole-genome overview…")
        except ValueError as exc:
            b.error(exc)

    def ready(self, result, job):
        self.overview, self.source = result, job["input"]
        stat = Path(self.source).stat()
        self.source_version = stat.st_mtime_ns, stat.st_size, stat.st_ino
        self.image.show_image(result["output"])
        self.image.set_map_data(result["map_data"])
        self.label.setText(f"Raw contacts · {result['binsize']:,} bp / bin · Double-click a chromosome block to open it")

    def hover(self, info):
        if info:
            x, y = info["x"], info["y"]
            self.label.setText(f"X {x[0]}:{x[1]:,}–{x[2]:,} · Y {y[0]}:{y[1]:,}–{y[2]:,} · Contact {info['value']:g}")

    def open_bin(self, info):
        b = self.browser
        if not info or b.busy() or self.source != b.input_path:
            return
        sizes = dict(b.metadata["chromsizes"])
        x, y = info["x"][0], info["y"][0]
        b.loading_controls = True
        b.map_display.setCurrentIndex(b.map_display.findData("observed"))
        b.mode.setCurrentIndex(b.mode.findData("square"))
        b.resolution.setCurrentText("auto")
        b.loading_controls = False
        self.hide()
        b.tabs.setCurrentIndex(0)
        b.navigate_viewport(Viewport(Region(x, 0, sizes[x]), Region(y, 0, sizes[y])))

    def export(self):
        b = self.browser
        if b.busy() or self.overview is None:
            return
        path, selected = QFileDialog.getSaveFileName(self, "Export overview figure", "genome.svg",
            "SVG (*.svg);;PNG (*.png);;PDF (*.pdf)", options=QFileDialog.Option.DontUseNativeDialog)
        if not path:
            return
        if not Path(path).suffix:
            path += ".png" if "PNG" in selected else ".pdf" if "PDF" in selected else ".svg"
        try:
            if Path(path).suffix.lower() == ".svg":
                shutil.copyfile(self.overview["output"], path)
            elif Path(path).suffix.lower() in {".png", ".pdf"}:
                b.start_job(dict(action="genome_export", input=self.source, output=path,
                                 result=self.overview, export=True), "Exporting whole-genome overview…")
            else:
                raise ValueError("Select SVG, PNG or PDF")
        except (ValueError, OSError) as exc:
            b.error(exc)

    def export_matrix(self):
        if self.overview is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export genome matrix and bin coordinates", "genome.npz", "NumPy (*.npz)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            try:
                import numpy as np
                data = self.image.map_data
                assert data is not None  # the overview render always attaches map data
                layout = self.overview["chromosomes"]
                resolution = self.overview["binsize"]
                chroms, starts, ends = [], [], []
                for c in layout:
                    for i in range(c["count"]):
                        chroms.append(c["name"])
                        starts.append(i*resolution)
                        ends.append(min((i+1)*resolution, c["length"]))
                np.savez(path, matrix=data.observed, chrom=np.array(chroms), start=starts, end=ends, binsize=resolution)
            except OSError as exc:
                self.browser.error(exc)
