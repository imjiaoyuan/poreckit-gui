"""Independent browser windows with optional genomic navigation linking."""
import shutil
import uuid
from pathlib import Path

from .model import Region
from .viewport import Viewport


def canonical_assembly(value):
    value = str(value or "").casefold().replace("-", "").replace("_", "")
    return {"hg19": "grch37", "hs37d5": "grch37", "hg38": "grch38",
            "mm10": "grcm38", "mm39": "grcm39"}.get(value, value)


def translated_viewport(viewport, source, target):
    a, b = canonical_assembly(source.get("assembly")), canonical_assembly(target.get("assembly"))
    if a and b and a != b:
        raise ValueError("Genome assemblies differ")
    def names(meta):
        return {name.removeprefix("chr"): (name, length) for name, length in meta.get("chromsizes", [])}
    src, dst = names(source), names(target)
    def axis(region):
        key = region.chrom.removeprefix("chr")
        if key not in src or key not in dst:
            raise ValueError("A linked chromosome is missing")
        name, length = dst[key]
        if src[key][1] != length:
            raise ValueError("Chromosome lengths differ")
        if region.end > length:
            raise ValueError("Linked region exceeds chromosome length")
        return Region(name, region.start, region.end)
    return Viewport(axis(viewport.x), axis(viewport.y))


class WindowController:
    # Strong ownership keeps windows created by menu actions alive. Closing a
    # window unregisters it after the close event has actually been accepted.
    instances: list["WindowController"] = []

    def __init__(self, browser):
        self.browser = browser
        self.previous = None
        self.instances.append(self)
        menu = browser.menuBar().addMenu("Window")
        browser.action(menu, "New window", lambda: self.create(), "Ctrl+N")
        browser.action(menu, "Duplicate view in new window", lambda: self.create(duplicate=True), "Ctrl+Shift+N")
        menu.addSeparator()
        self.link = browser.action(menu, "Link genomic navigation", lambda: None)
        self.link.setCheckable(True)
        self.link.setToolTip("Synchronize genomic regions between linked windows with matching genomes. Each window keeps its own data and display settings.")
        browser.action(menu, "Link all open windows", lambda: self.link_all(True))
        browser.action(menu, "Unlink all windows", lambda: self.link_all(False))
        self.window_menu = menu.addMenu("Activate window")
        self.window_menu.aboutToShow.connect(self.populate)

    def populate(self):
        self.window_menu.clear()
        for index, controller in enumerate(list(self.instances), 1):
            b = controller.browser
            label = f"{index} · {Path(b.input_path).name if b.input_path else 'Empty browser'}"
            action = self.browser.action(self.window_menu, label, lambda checked=False, window=b: self.activate(window))
            action.setCheckable(True)
            action.setChecked(b is self.browser)

    @staticmethod
    def activate(window):
        window.showNormal()
        window.raise_()
        window.activateWindow()

    def create(self, duplicate=False):
        from .states import file_references
        if duplicate and self.browser.tabs.currentWidget() is self.browser.sequence_browser:
            return self.browser.sequence_browser.new_window()
        if duplicate and (self.browser.busy() or not self.browser.input_path):
            return None
        b = type(self.browser)()
        try:
            if duplicate:
                view = self.browser.session_data()
                copies = {}
                def copy_temporary(path):
                    source = Path(path).resolve()
                    root = Path(self.browser.temp.name).resolve()
                    if not source.is_relative_to(root):
                        return str(source)
                    if source not in copies:
                        destination = Path(b.temp.name) / (uuid.uuid4().hex+"-"+source.name)
                        shutil.copyfile(source, destination)
                        copies[source] = str(destination)
                    return copies[source]
                b.apply_session(file_references(view, copy_temporary))
            b.show()
            b.move(self.browser.pos().x()+32, self.browser.pos().y()+32)
            return b
        except (ValueError, OSError, KeyError) as exc:
            b.close()
            self.browser.error(exc)
            return None

    @classmethod
    def link_all(cls, checked):
        for controller in list(cls.instances):
            controller.link.setChecked(checked)

    def completed(self, viewport, linked=False):
        changed = self.previous != viewport
        self.previous = viewport
        if not changed or linked or not self.link.isChecked():
            return
        for controller in list(self.instances):
            if controller is not self and controller.link.isChecked() and controller.browser.input_path:
                controller.receive(self.browser, viewport)

    def receive(self, source, viewport):
        b = self.browser
        try:
            target = translated_viewport(viewport, source.metadata, b.metadata)
            mode = b.map_display.currentData()
            if not target.diagonal and mode in {"split", "split_oe", "oe_local", "pearson_local", "pearson", "compartment"}:
                raise ValueError("This display requires identical X and Y regions")
            from .contactview import CONTROL_EXPECTED_MODES, EXPECTED_MODES
            if target.x.chrom != target.y.chrom and mode in EXPECTED_MODES | CONTROL_EXPECTED_MODES:
                raise ValueError("This display requires a single chromosome")
            text = b.resolution.currentText()
            if text != "auto" and max(r.end-r.start for r in (target.x, target.y)) > int(text)*2500:
                raise ValueError("Choose auto or a coarser resolution in this window")
            b.navigate_viewport(target, linked=True)
        except ValueError as exc:
            b.statusBar().showMessage(f"Linked navigation skipped · {exc}")

    def close(self):
        if self in self.instances:
            self.instances.remove(self)
