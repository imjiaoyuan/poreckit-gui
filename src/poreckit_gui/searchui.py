"""Search annotation files and navigate either genomic axis."""
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton,
    QLabel, QComboBox, QSpinBox, QTableWidget, QTableWidgetItem, QFileDialog, QAbstractItemView)

from .model import Region
from .viewport import Viewport


class GeneSearchDialog(QDialog):
    def __init__(self, browser):
        super().__init__(browser)
        self.browser = browser
        self.paths = []
        self.results = []
        self.setWindowTitle("Find gene")
        self.resize(820, 500)
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        self.query = QLineEdit()
        self.query.setPlaceholderText("Gene name, ID or alias")
        self.query.returnPressed.connect(self.search)
        row.addWidget(self.query, 1)
        button = QPushButton("Search")
        button.clicked.connect(self.search)
        row.addWidget(button)
        button = QPushButton("Add annotation file…")
        button.clicked.connect(self.add_file)
        row.addWidget(button)
        layout.addLayout(row)
        self.sources = QLabel("Searches loaded gene tracks and added GTF / GFF / BED files")
        self.sources.setWordWrap(True)
        self.sources.setTextFormat(Qt.TextFormat.PlainText)
        layout.addWidget(self.sources)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Gene", "Chromosome", "Start", "End", "Annotation file"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.itemDoubleClicked.connect(lambda _: self.navigate())
        self.table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table, 1)
        row = QHBoxLayout()
        row.addWidget(QLabel("Navigate"))
        self.axis = QComboBox()
        self.axis.addItems(["Both axes", "X axis", "Y axis"])
        row.addWidget(self.axis)
        row.addWidget(QLabel("Flanking bp"))
        self.padding = QSpinBox()
        self.padding.setRange(0, 100000000)
        self.padding.setValue(50000)
        row.addWidget(self.padding)
        button = QPushButton("Go to selected")
        button.clicked.connect(self.navigate)
        row.addWidget(button)
        self.status = QLabel("Enter a name and click Search")
        row.addWidget(self.status, 1)
        layout.addLayout(row)

    def all_paths(self):
        return list(dict.fromkeys([t["path"] for t in self.browser.track_data() if t["kind"] == "genes"] + self.paths))

    def open_search(self):
        self.show()
        self.raise_()
        self.query.setFocus()
        self.sources.setText("Annotation files: " + (", ".join(Path(p).name for p in self.all_paths()) or "Add a gene track or an annotation file"))

    def search_text(self, text):
        self.open_search()
        self.query.setText(text)
        self.search()

    def add_file(self):
        paths, _ = QFileDialog.getOpenFileNames(self, "Add gene search annotations", "",
            "Annotations (*.gtf *.gff *.gff3 *.bed *.gtf.gz *.gff.gz *.gff3.gz *.bed.gz);;All files (*)",
            options=QFileDialog.Option.DontUseNativeDialog)
        self.paths = list(dict.fromkeys(self.paths+paths))
        self.open_search()
        if self.query.text().strip():
            self.search()

    def search(self):
        b = self.browser
        if b.busy():
            self.status.setText("Wait for the current task or cancel it")
            return
        if not b.input_path:
            self.status.setText("Load contact data first")
            return
        if not self.all_paths():
            self.status.setText("Add a gene track or an annotation file first")
            return
        if not self.query.text().strip():
            self.status.setText("Enter a gene name, ID or alias")
            return
        self.status.setText("Searching annotations…")
        b.start_job(dict(action="gene_search", input=b.input_path, query=self.query.text(), paths=self.all_paths(),
                         chromsizes=b.metadata.get("chromsizes", [])), "Searching gene annotations…")

    def ready(self, result):
        self.results = result["results"]
        self.table.setRowCount(len(self.results))
        for row, gene in enumerate(self.results):
            for col, value in enumerate([gene["name"], gene["chrom"], f'{gene["start"]:,}', f'{gene["end"]:,}', Path(gene["source"]).name]):
                item = QTableWidgetItem(str(value))
                item.setToolTip("Aliases: " + ", ".join(gene["aliases"]) + "\n" + gene["source"])
                self.table.setItem(row, col, item)
        if self.results:
            self.table.selectRow(0)
        self.status.setText(f'{result["total"]} matches' + (f" · first {len(self.results)} shown" if result["total"] > len(self.results) else ""))

    def navigate(self):
        b = self.browser
        row = self.table.currentRow()
        if row < 0 or b.busy():
            return
        gene = self.results[row]
        sizes = dict(b.metadata.get("chromsizes", []))
        length = sizes.get(gene["chrom"])
        if length is None:
            self.status.setText("This chromosome is absent from the current contact map")
            return
        minimum = min(b.metadata.get("resolutions", []), default=1000)*20
        span = max(gene["end"]-gene["start"]+2*self.padding.value(), minimum)
        region = Region(gene["chrom"], gene["start"], gene["end"]).window((gene["start"]+gene["end"]-span)/2, span, length)
        current = b.current_viewport()
        target = (Viewport(region, region), Viewport(region, current.y), Viewport(current.x, region))[self.axis.currentIndex()]
        self.hide()
        b.tabs.setCurrentIndex(0)
        b.resolution.setCurrentText("auto")
        b.navigate_viewport(target)
