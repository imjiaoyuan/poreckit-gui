"""Full plotting catalog, editable track stacks and a shared figure viewer."""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

from PySide6.QtCore import QProcess, QProcessEnvironment, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .catalog import Entry, catalog
from .forms import ParameterForm
from .jobs import JobProcess
from .model import Region
from .widgets import ImageView

DIALOG = QFileDialog.Option.DontUseNativeDialog


class PlotWorkbench(QWidget):
    result: dict
    # Widgets/attributes created or replaced after construction; every use
    # below happens after the corresponding setup has run.
    selected: Entry
    form: ParameterForm

    def __init__(self, browser):
        super().__init__()
        self.browser = browser
        self.entries = catalog()
        self.temp = tempfile.TemporaryDirectory(prefix="poreckit-plots-")
        self.job = None
        self.serial = 0
        self.result = None  # type: ignore[assignment]  # cleared until a plot result arrives
        self.form = None  # type: ignore[assignment]  # ParameterForm built when an entry is selected
        self.selected = None  # type: ignore[assignment]  # Entry set in selection_changed
        self.cancelled = False
        self.form_cache = {}
        self.process = JobProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.read_log)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.process_error)
        self._build()

    def _build(self):
        outer = QVBoxLayout(self)
        context_row = QHBoxLayout()
        context_row.addWidget(QLabel("Locus"))
        self.locus = QLineEdit("chr1:0-2,000,000")
        self.locus.setToolTip("Track-stack region; zero-based, half-open. Command/function region parameters are in the form.")
        context_row.addWidget(self.locus)
        apply_region = QPushButton("Apply locus")
        apply_region.setToolTip("Copy this locus into region parameters and all tracks before rendering")
        apply_region.clicked.connect(self.apply_locus)
        context_row.addWidget(apply_region)
        self.use_browser = QPushButton("Use browser data / region")
        self.use_browser.clicked.connect(self.prefill)
        context_row.addWidget(self.use_browser)
        context_row.addWidget(QLabel("Data folder"))
        self.cwd = QLineEdit(str(Path.cwd()))
        self.cwd.setToolTip("Relative input paths and structured data references resolve from this folder")
        context_row.addWidget(self.cwd)
        folder = QPushButton("Browse…")
        folder.clicked.connect(self.choose_folder)
        context_row.addWidget(folder)
        outer.addLayout(context_row)
        splitter = QSplitter()
        outer.addWidget(splitter, 1)
        navigation = QWidget()
        nav = QVBoxLayout(navigation)
        nav.setContentsMargins(0, 0, 0, 0)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Find a plot or track…")
        self.search.textChanged.connect(self.filter_catalog)
        nav.addWidget(self.search)
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setMinimumWidth(180)
        groups = {}
        for entry in self.entries.values():
            if entry.group not in groups:
                groups[entry.group] = QTreeWidgetItem(self.tree, [entry.group])
            item = QTreeWidgetItem(groups[entry.group], [entry.title])
            item.setData(0, Qt.ItemDataRole.UserRole, entry.key)
            item.setToolTip(0, entry.description)
        self.tree.itemSelectionChanged.connect(self.selection_changed)
        self.tree.expandToDepth(0)
        nav.addWidget(self.tree, 1)
        counts = {kind: sum(e.kind == kind for e in self.entries.values()) for kind in ("cli", "api", "track")}
        nav.addWidget(QLabel(f"{counts['cli']} commands · {counts['api']} functions\n{counts['track']} track types"))
        splitter.addWidget(navigation)

        self.settings_panel = QWidget()
        settings = QVBoxLayout(self.settings_panel)
        settings.setContentsMargins(0, 0, 0, 0)
        self.heading = QLabel("Choose a plot")
        self.heading.setWordWrap(True)
        self.heading.setTextFormat(Qt.TextFormat.PlainText)
        settings.addWidget(self.heading)
        self.description = QPlainTextEdit()
        self.description.setReadOnly(True)
        self.description.setMaximumHeight(95)
        self.description.setPlaceholderText("Select a command, function or track from the catalog.")
        settings.addWidget(self.description)
        self.forms = QStackedWidget()
        self.forms.addWidget(QWidget())
        settings.addWidget(self.forms, 1)
        self.stack_label = QLabel("Track stack · Check to show, drag to reorder")
        self.stack_label.setWordWrap(True)
        settings.addWidget(self.stack_label)
        self.track_stack = QListWidget()
        self.track_stack.setMaximumHeight(130)
        self.track_stack.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self.track_stack.itemDoubleClicked.connect(self.edit_track)
        settings.addWidget(self.track_stack)
        row = QHBoxLayout()
        self.add_track_button = QPushButton("Add / update track")
        self.add_track_button.clicked.connect(self.add_track)
        row.addWidget(self.add_track_button)
        remove = QPushButton("Remove")
        self.remove_track_button = remove
        remove.clicked.connect(self.remove_track)
        row.addWidget(remove)
        settings.addLayout(row)
        row = QHBoxLayout()
        save = QPushButton("Save preset…")
        save.clicked.connect(self.save_preset)
        load = QPushButton("Open preset…")
        load.clicked.connect(self.load_preset)
        row.addWidget(save)
        row.addWidget(load)
        settings.addLayout(row)
        self.run_button = QPushButton("Render plot")
        self.run_button.clicked.connect(self.run)
        settings.addWidget(self.run_button)
        splitter.addWidget(self.settings_panel)

        results = QWidget()
        display = QVBoxLayout(results)
        display.setContentsMargins(0, 0, 0, 0)
        self.figure_list = QComboBox()
        self.figure_list.setToolTip("Every figure from the selected operation is available here")
        self.figure_list.currentIndexChanged.connect(self.show_figure)
        display.addWidget(self.figure_list)
        tools = QHBoxLayout()
        self.image = ImageView()
        for text, callback in [("Zoom +", lambda: self.image.zoom_by(1.25)),
                               ("Zoom −", lambda: self.image.zoom_by(0.8)),
                               ("Fit", self.image.fit), ("100%", self.image.actual_size),
                               ("Export…", self.export), ("Viewport…", self.export_viewport),
                               ("Export all…", self.export_all)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            tools.addWidget(button)
        display.addLayout(tools)
        display.addWidget(self.image, 1)
        self.zoom_label = QLabel("Scroll to zoom · Drag to pan · Exports preserve the full figure")
        self.zoom_label.setWordWrap(True)
        self.image.zoomChanged.connect(lambda ratio: self.zoom_label.setText(f"Zoom {ratio:.0%} of fit · Drag to pan · Viewport exports the visible area"))
        display.addWidget(self.zoom_label)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(1500)
        self.log.setMaximumHeight(110)
        display.addWidget(self.log)
        splitter.addWidget(results)
        splitter.setSizes([220, 340, 800])
        splitter.setStretchFactor(2, 1)
        status = QHBoxLayout()
        self.status = QLabel("Ready · Choose a plot or build a track stack")
        self.status.setWordWrap(True)
        status.addWidget(self.status, 1)
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(120)
        self.progress.setRange(0, 0)
        self.progress.hide()
        status.addWidget(self.progress)
        self.stop = QPushButton("Cancel")
        self.stop.clicked.connect(self.cancel)
        self.stop.hide()
        status.addWidget(self.stop)
        outer.addLayout(status)
        self.editing_track = None
        self.select_entry("cli:heatmap")

    def apply_locus(self):
        if self.busy():
            return
        try:
            region = Region.parse(self.locus.text())
            values = {"chrom": region.chrom, "start": region.start, "end": region.end}
            if self.form:
                self.form.set_values(values)
            from .catalog import fields
            for index in range(self.track_stack.count()):
                item = self.track_stack.item(index)
                track = item.data(Qt.ItemDataRole.UserRole)
                names = {f["name"] for f in fields(self.entries[track["entry"]])}
                track["values"].update({k: v for k, v in values.items() if k in names})
                item.setData(Qt.ItemDataRole.UserRole, track)
            self.status.setText("Locus applied · Click Render to query this region")
        except ValueError as exc:
            self.error(exc)

    def error(self, message):
        QMessageBox.warning(self, "poreckit", str(message))

    def busy(self):
        return self.job is not None

    def choose_folder(self):
        path = QFileDialog.getExistingDirectory(self, "Data folder", self.cwd.text(), options=DIALOG)
        if path:
            self.cwd.setText(path)

    def filter_catalog(self, text):
        needle = text.lower().strip()
        for index in range(self.tree.topLevelItemCount()):
            group = self.tree.topLevelItem(index)
            assert group is not None  # index < topLevelItemCount()
            visible = False
            for child_index in range(group.childCount()):
                item = group.child(child_index)
                assert item is not None  # child_index < childCount()
                entry = self.entries[item.data(0, Qt.ItemDataRole.UserRole)]
                match = not needle or needle in (entry.title + " " + entry.group + " " + entry.description).lower()
                item.setHidden(not match)
                visible |= match
            group.setHidden(not visible)
            if needle and visible:
                group.setExpanded(True)

    def select_entry(self, key):
        if key not in self.entries:
            raise ValueError(f"Unknown plot entry: {key}")
        for index in range(self.tree.topLevelItemCount()):
            group = self.tree.topLevelItem(index)
            assert group is not None  # index < topLevelItemCount()
            for child_index in range(group.childCount()):
                item = group.child(child_index)
                assert item is not None  # child_index < childCount()
                if item.data(0, Qt.ItemDataRole.UserRole) == key:
                    self.tree.setCurrentItem(item)
                    group.setExpanded(True)
                    return

    def selection_changed(self):
        selected = self.tree.selectedItems()
        if not selected:
            return
        key = selected[0].data(0, Qt.ItemDataRole.UserRole)
        if not key:
            return
        self.selected = self.entries[key]
        if key not in self.form_cache:
            form = ParameterForm(self.selected)
            self.form_cache[key] = form
            self.forms.addWidget(form)
        self.form = self.form_cache[key]
        self.forms.setCurrentWidget(self.form)
        self.heading.setText(self.selected.title)
        help_text = self.selected.description
        if self.selected.kind in {"api", "track"}:
            help_text += ('\n\nData inputs accept JSON or NPY/NPZ files. Enable Array for inline numeric arrays. '
                          'Use {"$array": [...]} for arrays nested in result dictionaries, '
                          'or {"$file": "path.npz", "key": "matrix"} for a specific array. '
                          'Current matrix queries the Genome Browser input at the selected locus.')
        self.description.setPlainText(help_text)
        self.add_track_button.setEnabled(self.selected.kind == "track")
        for widget in (self.stack_label, self.track_stack, self.add_track_button, self.remove_track_button):
            widget.setVisible(self.selected.kind == "track")
        self.run_button.setText("Render track stack" if self.selected.kind == "track" else "Render plot")
        self.editing_track = None

    def context(self):
        context = dict(self.browser.last_view or {"binsize": None, "norm": None})
        if self.browser.input_path:
            context["input"] = self.browser.input_path
        context["region"] = Region.parse(self.locus.text()).text()
        return context

    def prefill(self):
        if self.browser.region:
            self.locus.setText(self.browser.locus.text())
        if not self.form:
            return
        try:
            context = self.context()
            region = Region.parse(context["region"])
            resolution = context.get("binsize")
            if resolution is None and self.browser.metadata.get("resolutions"):
                available = sorted(self.browser.metadata["resolutions"])
                target = max(1, (region.end - region.start + 249) // 250)
                resolution = next((n for n in available if n >= target), available[-1])
            resolution = resolution or max(1, (region.end - region.start + 249) // 250)
            values = {"chrom": region.chrom, "start": region.start, "end": region.end,
                      "bin_size": resolution, "binsize": resolution}
            if context.get("norm"):
                values["norm"] = context["norm"]
            if context.get("input"):
                for name in ("contacts", "pairs", "source", "input", "matrix1"):
                    if self.selected.kind == "cli":
                        values[name] = context["input"]
                for name in ("reader", "reader1"):
                    values[name] = {"$reader": context["input"]}
                for name in ("pairs_path", "bam_path"):
                    if name == "pairs_path" and context["input"].endswith((".pairs", ".pairs.gz", ".parquet")):
                        values[name] = context["input"]
                    elif name == "bam_path" and context["input"].endswith(".bam"):
                        values[name] = context["input"]
            self.form.set_values(values)
        except (ValueError, KeyError) as exc:
            self.error(exc)

    def add_track(self):
        if self.busy() or not self.selected or self.selected.kind != "track":
            return
        try:
            data = {"entry": self.selected.key, "values": self.form.values(), "visible": True}
            item = self.editing_track
            if item is None:
                item = QListWidgetItem()
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(Qt.CheckState.Checked)
                self.track_stack.addItem(item)
            item.setText(data["values"].get("name") or self.selected.title)
            item.setData(Qt.ItemDataRole.UserRole, data)
            self.editing_track = None
            self.status.setText("Track saved · Render the stack to apply changes")
        except (ValueError, TypeError) as exc:
            self.error(exc)

    def edit_track(self, item):
        data = item.data(Qt.ItemDataRole.UserRole)
        self.select_entry(data["entry"])
        # A new form avoids leaving values from a different instance of this type.
        self.reset_form(data["entry"], data["values"])
        self.editing_track = item

    def remove_track(self):
        self.editing_track = None
        self.track_stack.takeItem(self.track_stack.currentRow())

    def tracks(self):
        return [{**self.track_stack.item(i).data(Qt.ItemDataRole.UserRole),
                 "visible": self.track_stack.item(i).checkState() == Qt.CheckState.Checked}
                for i in range(self.track_stack.count())]

    def run(self):
        if self.busy() or not self.selected:
            return
        try:
            context = self.context()
            values = self.form.values()
            tracks = self.tracks()
            if self.selected.kind == "track" and not tracks:
                tracks = [{"entry": self.selected.key, "values": values, "visible": True}]
            cwd = Path(self.cwd.text()).expanduser().resolve()
            if not cwd.is_dir():
                raise ValueError("Select a valid data folder")
            self.serial += 1
            directory = Path(self.temp.name) / f"run-{self.serial}"
            directory.mkdir()
            job = {"action": "plot", "kind": self.selected.kind, "entry": self.selected.key,
                   "values": values, "tracks": tracks, "context": context, "cwd": str(cwd), "directory": str(directory)}
            self.start(job)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            self.error(exc)

    def start(self, job):
        self.job = job
        self.cancelled = False
        job_path = Path(self.temp.name) / "job.json"
        self.result_path = Path(self.temp.name) / "result.json"
        self.result_path.unlink(missing_ok=True)
        job_path.write_text(json.dumps(job, ensure_ascii=False), encoding="utf-8")
        self.settings_panel.setEnabled(False)
        self.tree.setEnabled(False)
        self.use_browser.setEnabled(False)
        self.progress.show()
        self.stop.show()
        self.status.setText("Rendering · Previous results remain available")
        self.log.appendPlainText(f"\nRender: {job['entry']}")
        env = QProcessEnvironment.systemEnvironment()
        env.insert("MPLBACKEND", "Agg")
        env.insert("MPLCONFIGDIR", str(Path(self.temp.name) / "matplotlib"))
        env.insert("PYTHONUNBUFFERED", "1")
        self.process.setProcessEnvironment(env)
        args = ["--worker"] if getattr(sys, "frozen", False) else ["-m", "poreckit_gui.worker"]
        self.process.start(sys.executable, [*args, str(job_path), str(self.result_path)])

    def read_log(self):
        value = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if value:
            self.log.appendPlainText(value.rstrip())

    def process_error(self, error):
        if error == QProcess.ProcessError.FailedToStart and self.job:
            self.log.appendPlainText(self.process.errorString())
            self.finished(-1, QProcess.ExitStatus.CrashExit)

    def reap_runs(self):
        """Keep only the current and previous run directories; older ones are
        removed during the session instead of at shutdown."""
        import re
        keep = {self.serial, max(1, self.serial - 1)}
        for path in Path(self.temp.name).glob("run-*"):
            match = re.fullmatch(r"run-(\d+)", path.name)
            if match and int(match.group(1)) not in keep:
                shutil.rmtree(path, ignore_errors=True)

    def discard_run(self, job):
        directory = Path(job["directory"]) if job.get("directory") else None
        if directory and directory.is_dir():
            shutil.rmtree(directory, ignore_errors=True)

    def finished(self, code, status):
        self.read_log()
        job, self.job = self.job, None
        self.settings_panel.setEnabled(True)
        self.tree.setEnabled(True)
        self.use_browser.setEnabled(True)
        self.progress.hide()
        self.stop.hide()
        if job is None:
            return
        if self.cancelled:
            # A cancelled run leaves no usable figures; drop its directory.
            self.discard_run(job)
            self.status.setText("Task cancelled · Previous results retained")
            return
        if code != 0 or not self.result_path.exists():
            self.discard_run(job)
            self.status.setText("Rendering failed · See the run log")
            self.error("Rendering failed. See the run log for the original poreckit error.")
            return
        try:
            result = json.loads(self.result_path.read_text(encoding="utf-8"))
            if not result.get("figures"):
                raise ValueError("No figures in result")
            self.result = result
            self.reap_runs()
            self.figure_list.blockSignals(True)
            self.figure_list.clear()
            for index, figure in enumerate(result["figures"]):
                self.figure_list.addItem(f"{index + 1}. {figure['title']}")
            self.figure_list.blockSignals(False)
            self.show_figure(0)
            count = len(result["figures"])
            self.status.setText(f"Ready · {count} {'figure' if count == 1 else 'figures'} · {len(result['sidecars'])} accompanying files")
        except (ValueError, KeyError, OSError) as exc:
            self.error(exc)

    def show_figure(self, index):
        if self.result and 0 <= index < len(self.result["figures"]):
            try:
                files = self.result["figures"][index]["files"]
                self.image.show_image(files.get("svg") or files["png"])
            except (ValueError, OSError) as exc:
                self.error(exc)

    def export(self):
        if not self.result:
            return
        path, selected = QFileDialog.getSaveFileName(self, "Export full figure", "poreckit.svg",
            "SVG (*.svg);;PDF (*.pdf);;PNG (*.png)", options=DIALOG)
        if path:
            try:
                suffix = Path(path).suffix.lower().lstrip(".")
                if not suffix:
                    suffix = {"SVG (*.svg)": "svg", "PDF (*.pdf)": "pdf", "PNG (*.png)": "png"}[selected]
                    path += "." + suffix
                files = self.result["figures"][self.figure_list.currentIndex()]["files"]
                if suffix not in files:
                    raise ValueError("Select SVG, PDF or PNG")
                shutil.copyfile(files[suffix], path)
                self.status.setText(f"Exported: {path}")
            except (ValueError, OSError) as exc:
                self.error(exc)

    def export_viewport(self):
        if not self.result:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export visible viewport", "viewport.png", "PNG (*.png)", options=DIALOG)
        if path:
            path = path if path.lower().endswith(".png") else path + ".png"
            if not self.image.viewport().grab().save(path):
                self.error("Unable to save the viewport")

    def export_all(self):
        if not self.result:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export all figures and accompanying files", "poreckit-results.zip", "ZIP (*.zip)", options=DIALOG)
        if path:
            try:
                path = path if path.lower().endswith(".zip") else path + ".zip"
                files = {Path(p) for fig in self.result["figures"] for p in fig["files"].values()}
                files.update(Path(p) for p in self.result["sidecars"])
                root = Path(self.result["figures"][0]["files"]["png"]).parent
                with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
                    for file in sorted(files):
                        archive.write(file, file.relative_to(root))
                self.status.setText(f"Exported result set: {path}")
            except OSError as exc:
                self.error(exc)

    def reset_form(self, key, values):
        if key in self.form_cache:
            old = self.form_cache.pop(key)
            self.forms.removeWidget(old)
            old.deleteLater()
        self.form = ParameterForm(self.entries[key])
        self.form_cache[key] = self.form
        self.forms.addWidget(self.form)
        self.forms.setCurrentWidget(self.form)
        self.form.set_values(values)

    def save_preset(self):
        if self.busy() or not self.selected:
            return
        try:
            data = {"version": 2, "entry": self.selected.key, "values": self.form.values(),
                    "tracks": self.tracks(), "cwd": self.cwd.text(), "locus": self.locus.text()}
            path, _ = QFileDialog.getSaveFileName(self, "Save plot preset", "plot-preset.json", "JSON (*.json)", options=DIALOG)
            if path:
                from .sessionio import save_json
                save_json(path if Path(path).suffix else path + ".json", data)
        except (ValueError, OSError) as exc:
            self.error(exc)

    def load_preset(self):
        if self.busy():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open plot preset", "", "JSON (*.json)", options=DIALOG)
        if not path:
            return
        try:
            from .sessionio import load_json
            data = load_json(path)
            if data["version"] != 2 or data["entry"] not in self.entries or not isinstance(data["values"], dict):
                raise ValueError("Unsupported plot preset")
            Region.parse(data["locus"])
            for track in data["tracks"]:
                if self.entries[track["entry"]].kind != "track" or not isinstance(track["values"], dict):
                    raise ValueError("Invalid track stack")
            self.select_entry(data["entry"])
            self.reset_form(data["entry"], data["values"])
            self.locus.setText(data["locus"])
            self.cwd.setText(data["cwd"])
            self.track_stack.clear()
            for track in data["tracks"]:
                item = QListWidgetItem(track["values"].get("name") or self.entries[track["entry"]].title)
                item.setData(Qt.ItemDataRole.UserRole, track)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(Qt.CheckState.Checked if track.get("visible", True) else Qt.CheckState.Unchecked)
                self.track_stack.addItem(item)
            self.status.setText("Preset loaded · Review parameters and click Render")
        except (ValueError, KeyError, TypeError, OSError) as exc:
            self.error(f"Unable to open preset: {exc}")

    def cancel(self):
        if self.busy():
            self.cancelled = True
            self.process.kill()

    def shutdown(self):
        self.process.blockSignals(True)
        if self.busy():
            self.process.kill()
            self.process.waitForFinished(3000)
            self.job = None
        self.process.shutdown()
        from .sessionio import cleanup_temp
        cleanup_temp(self.temp.name)
