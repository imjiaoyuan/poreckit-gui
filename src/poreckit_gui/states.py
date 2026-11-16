"""Named browser states and timed navigation."""
import copy
import json
import os
import uuid
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QTimer
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDockWidget,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


def file_references(state, transform):
    state = copy.deepcopy(state)
    for key in ("input", "control"):
        if state.get(key):
            state[key] = transform(state[key])
    for track in state.get("tracks", []):
        track["path"] = transform(track["path"])
    if "gene_annotations" in state:
        state["gene_annotations"] = [transform(path) for path in state["gene_annotations"]]
    for saved in state.get("bookmarks", []):
        saved["view"] = file_references(saved["view"], transform)
    return state


def relative_references(state, folder):
    def relative(path):
        try:
            return os.path.relpath(path, folder)
        except ValueError:
            return path
    return file_references(state, relative)


def absolute_references(state, folder):
    return file_references(state, lambda path: str((Path(folder) / path).resolve()))


class StateController(QObject):
    def __init__(self, browser):
        super().__init__(browser)
        self.browser = browser
        self.entries = []
        self.saved = "[]"
        self.timer = QTimer(self)
        self.timer.setInterval(5000)
        self.timer.timeout.connect(self.advance)
        menu = browser.menuBar().addMenu("Bookmarks")
        browser.action(menu, "Bookmark current view…", self.add, "Ctrl+B")
        self.dock = QDockWidget("Bookmarks", browser)
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.list = QListWidget()
        self.list.itemDoubleClicked.connect(lambda item: self.open_current())
        layout.addWidget(self.list)
        row = QHBoxLayout()
        for text, callback in [("Open", self.open_current), ("Add…", self.add), ("Update", self.update),
                               ("Rename…", self.rename), ("Remove", self.remove)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            row.addWidget(button)
        layout.addLayout(row)
        row = QHBoxLayout()
        for text, callback in [("Up", lambda: self.move(-1)), ("Down", lambda: self.move(1)),
                               ("Import…", self.load), ("Export…", self.export)]:
            button = QPushButton(text)
            button.clicked.connect(callback)
            row.addWidget(button)
        layout.addLayout(row)
        row = QHBoxLayout()
        self.play = QPushButton("Start slideshow")
        self.play.setCheckable(True)
        self.play.toggled.connect(self.play_changed)
        row.addWidget(self.play)
        self.next = QPushButton("Next")
        self.next.clicked.connect(self.advance)
        row.addWidget(self.next)
        row.addWidget(QLabel("Seconds"))
        self.interval = QSpinBox()
        self.interval.setRange(1, 3600)
        self.interval.setValue(5)
        self.interval.valueChanged.connect(lambda seconds: self.timer.setInterval(seconds*1000))
        row.addWidget(self.interval)
        self.indicator = QLabel("")
        self.indicator.hide()
        layout.addWidget(self.indicator)
        self.escape = QShortcut(QKeySequence("Esc"), panel)
        self.escape.activated.connect(lambda: self.play.setChecked(False))
        layout.addLayout(row)
        self.dock.setWidget(panel)
        browser.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.dock)
        self.dock.hide()
        menu.addAction(self.dock.toggleViewAction())
        browser.action(menu, "Import bookmarks…", self.load)
        browser.action(menu, "Export bookmarks…", self.export)
        self.restore_settings()

    @property
    def settings_key(self):
        return "bookmarks"

    def restore_settings(self):
        raw = self.browser.settings.value(self.settings_key)
        if not raw:
            return
        try:
            entries = json.loads(raw)
            valid = []
            if isinstance(entries, list):
                for entry in entries:
                    if (isinstance(entry, dict) and isinstance(entry.get("name"), str)
                            and isinstance(entry.get("view"), dict)):
                        valid.append(entry)
            self.entries = valid
            self.refresh()
            self.mark_saved()
        except (ValueError, TypeError):
            pass

    def persist(self):
        self.browser.settings.setValue(self.settings_key, json.dumps(self.entries))

    @property
    def dirty(self):
        return json.dumps(self.entries, sort_keys=True) != self.saved

    def mark_saved(self):
        self.saved = json.dumps(self.entries, sort_keys=True)

    def restore(self, entries):
        for entry in entries:
            if not isinstance(entry.get("name"), str) or not isinstance(entry.get("view"), dict):
                raise ValueError("Invalid bookmark")
            self.browser.validate_session(entry["view"])
        self.entries = copy.deepcopy(entries)
        self.refresh()
        self.mark_saved()
        self.persist()

    def refresh(self, index=None):
        index = self.list.currentRow() if index is None else index
        self.list.clear()
        for entry in self.entries:
            item = QListWidgetItem(entry["name"])
            view = entry["view"]
            item.setToolTip(f'{Path(view["input"]).name}\nX: {view["region"]}\nY: {view.get("region_y") or view["region"]}')
            self.list.addItem(item)
        if self.entries:
            self.list.setCurrentRow(max(0, min(index, len(self.entries)-1)))

    def add(self):
        if self.browser.busy() or not self.browser.input_path:
            return
        name, ok = QInputDialog.getText(self.browser, "Bookmark current view", "Name")
        if ok:
            try:
                self.entries.append(dict(id=uuid.uuid4().hex, name=name or self.browser.locus.text(),
                                         view=self.browser.session_data()))
                self.refresh(len(self.entries)-1)
                self.dock.show()
                self.persist()
            except (ValueError, OSError) as exc:
                self.browser.error(exc)

    def update(self):
        row = self.list.currentRow()
        if row >= 0 and not self.browser.busy():
            try:
                self.entries[row]["view"] = self.browser.session_data()
                self.refresh()
                self.persist()
            except (ValueError, OSError) as exc:
                self.browser.error(exc)

    def rename(self):
        row = self.list.currentRow()
        if row >= 0:
            name, ok = QInputDialog.getText(self.browser, "Rename bookmark", "Name", text=self.entries[row]["name"])
            if ok and name.strip():
                self.entries[row]["name"] = name
                self.refresh()
                self.persist()

    def remove(self):
        row = self.list.currentRow()
        if row >= 0:
            del self.entries[row]
            self.refresh(row)
            self.persist()

    def move(self, delta):
        row = self.list.currentRow()
        target = row+delta
        if row >= 0 and 0 <= target < len(self.entries):
            self.entries.insert(target, self.entries.pop(row))
            self.refresh(target)
            self.persist()

    def open_current(self):
        row = self.list.currentRow()
        if row < 0 or self.browser.busy():
            return
        try:
            self.browser.apply_session(copy.deepcopy(self.entries[row]["view"]))
        except (ValueError, OSError, KeyError) as exc:
            self.play.setChecked(False)
            self.browser.error(exc)

    def play_changed(self, enabled):
        if enabled and not self.entries:
            self.play.setChecked(False)
            return
        self.play.setText("Stop slideshow" if enabled else "Start slideshow")
        if enabled:
            self.open_current()
            self.update_indicator()
            self.timer.start()
        else:
            self.timer.stop()
            self.indicator.hide()

    def update_indicator(self):
        row = self.list.currentRow()
        if self.play.isChecked() and 0 <= row < len(self.entries):
            self.indicator.setText(f"{row+1}/{len(self.entries)} — {self.entries[row]['name']}")
            self.indicator.show()

    def advance(self):
        if self.browser.busy() or not self.entries:
            return
        row = self.list.currentRow()+1
        if row >= len(self.entries):
            # No wraparound: the slideshow stops at the last bookmark.
            self.play.setChecked(False)
            self.browser.statusBar().showMessage("Slideshow finished")
            return
        self.list.setCurrentRow(row)
        self.open_current()
        self.update_indicator()

    def export(self):
        path, _ = QFileDialog.getSaveFileName(self.browser, "Export bookmarks", "bookmarks.json", "JSON (*.json)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            try:
                data = relative_references({"version": 1, "bookmarks": self.entries}, Path(path).parent)
                Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")
                self.mark_saved()
            except (ValueError, OSError) as exc:
                self.browser.error(exc)

    def load(self):
        path, _ = QFileDialog.getOpenFileName(self.browser, "Import bookmarks", "", "JSON (*.json)",
                                             options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            try:
                data = absolute_references(json.loads(Path(path).read_text(encoding="utf-8")), Path(path).parent)
                if data.get("version") != 1:
                    raise ValueError("Unsupported bookmark version")
                entries = data["bookmarks"]
                valid, skipped = [], []
                for entry in entries:
                    try:
                        if not isinstance(entry.get("name"), str) or not isinstance(entry.get("view"), dict):
                            raise ValueError("Invalid bookmark")
                        self.browser.validate_session(entry["view"])
                        valid.append(entry)
                    except (ValueError, OSError, KeyError, TypeError, AttributeError) as exc:
                        # Skip bookmarks whose data files are missing or whose
                        # view state is unusable; never abort the whole import.
                        skipped.append(f'{entry.get("name", "Unnamed")}: {exc}')
                self.entries.extend(valid)
                self.refresh()
                self.dock.show()
                self.persist()
                if skipped:
                    self.browser.error("Some bookmarks were skipped:\n" + "\n".join(skipped))
            except (ValueError, OSError, KeyError, TypeError, AttributeError) as exc:
                self.browser.error(exc)
