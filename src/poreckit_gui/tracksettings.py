"""Per-track display options for the guided genome browser."""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QFormLayout, QLineEdit, QComboBox,
                              QDoubleSpinBox, QDialogButtonBox, QMessageBox)


class TrackSettings(QDialog):
    def __init__(self, track, parent=None):
        super().__init__(parent)
        self.track = dict(track)
        self.setWindowTitle("Track settings")
        layout = QVBoxLayout(self)
        form = QFormLayout()
        layout.addLayout(form)
        self.name = QLineEdit(track.get("name") or Path(track["path"]).name)
        form.addRow("Name", self.name)
        self.chrom = QLineEdit(track.get("chrom", ""))
        self.chrom.setPlaceholderText("Use main chromosome")
        self.chrom.setToolTip("Optional source name, such as chr1 instead of 1. Coordinates are not converted.")
        form.addRow("Source chromosome", self.chrom)
        self.chrom_y = QLineEdit(track.get("chrom_y", ""))
        self.chrom_y.setPlaceholderText("Use Y chromosome")
        form.addRow("Y source chromosome", self.chrom_y)
        self.mode = QComboBox()
        modes = {"signal": [("Filled", "fill"), ("Line", "line"), ("Bars", "bar"), ("Points", "points")],
                 "contacts": [("Triangle", "triangle"), ("Square", "square")],
                 "tads": [("Triangles", "triangle"), ("Boxes", "box")],
                 "loops": [("Anchor boxes / arcs", "box"), ("Arcs", "arc"), ("Bezier curves", "bezier"), ("Anchor rings / points", "point")],
                 "genes": [("Arrows", "arrow"), ("Basic", "basic")]}
        for label, value in modes[track["kind"]]:
            self.mode.addItem(label, value)
        self.mode.setCurrentIndex(max(0, self.mode.findData(track.get("style"))))
        form.addRow("Display mode", self.mode)
        self.placement = QComboBox()
        self.placement.addItem("Follow map layout", "auto")
        self.placement.addItem("Separate panel below", "panel")
        if track["kind"] in {"genes", "signal"}:
            self.placement.addItem("Above main map", "top")
            self.placement.addItem("Left of main map", "left")
            self.placement.addItem("Above and left of main map", "both")
        self.placement.setToolTip("Follow map layout places genes/signals above and left of square maps, and below triangle maps. TAD/loop tracks overlay square maps and stack below triangle maps.")
        if track["kind"] in {"tads", "loops"}:
            self.placement.addItem("Overlay main contact map", "overlay")
        default_placement = "overlay" if track["kind"] in {"tads", "loops"} else "panel"
        self.placement.setCurrentIndex(max(0, self.placement.findData(track.get("placement", default_placement))))
        form.addRow("Placement", self.placement)
        self.mode.setToolTip("In-map TADs use diagonal boxes on square maps. In-map loops use anchor boxes or rings; standalone triangle tracks use arcs or points.")
        self.color = QLineEdit(track.get("color", "#397ca8"))
        self.color.setPlaceholderText("Color name or #RRGGBB")
        if track["kind"] != "contacts":
            form.addRow("Color", self.color)
        self.cmap = QComboBox()
        self.cmap.addItem("Use main map colors", None)
        for name in ("hic_whitered", "Reds", "inferno", "viridis", "magma", "YlOrRd", "RdBu_r"):
            self.cmap.addItem(name, name)
        self.cmap.setCurrentIndex(max(0, self.cmap.findData(track.get("cmap"))))
        self.resolution = QLineEdit(str(track.get("binsize", "auto")))
        self.resolution.setPlaceholderText("auto or a resolution in bp")
        self.norm = QComboBox()
        self.norm.addItem("Use main normalization", None)
        for value in ("raw", "balanced", "KR", "VC", "VC_SQRT"):
            self.norm.addItem(value, value)
        self.norm.setCurrentIndex(max(0, self.norm.findData(track.get("norm"))))
        if track["kind"] == "contacts":
            form.addRow("Colormap", self.cmap)
            form.addRow("Resolution / bp", self.resolution)
            form.addRow("Normalization", self.norm)
        self.height = QDoubleSpinBox()
        self.height.setRange(0.1, 4)
        self.height.setSingleStep(0.1)
        self.height.setValue(float(track.get("height", 1.0 if track["kind"] == "contacts" else 0.3)))
        form.addRow("Relative panel height", self.height)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def accept(self):
        from matplotlib.colors import is_color_like
        if self.track["kind"] != "contacts" and not is_color_like(self.color.text()):
            QMessageBox.warning(self, "Track settings", "Enter a valid color name or #RRGGBB value")
            return
        self.track.update(name=self.name.text().strip() or Path(self.track["path"]).name,
                          style=self.mode.currentData(), placement=self.placement.currentData(),
                          color=self.color.text(), height=self.height.value())
        self.track["chrom"] = self.chrom.text().strip()
        self.track["chrom_y"] = self.chrom_y.text().strip()
        if self.track["kind"] == "contacts":
            value = self.resolution.text().strip()
            if value != "auto" and (not value.isdigit() or int(value) <= 0):
                QMessageBox.warning(self, "Track settings", "Resolution must be auto or a positive integer")
                return
            self.track.update(cmap=self.cmap.currentData(), binsize=value)
            if self.norm.currentData() is None:
                self.track.pop("norm", None)
            else:
                self.track["norm"] = self.norm.currentData()
        super().accept()
