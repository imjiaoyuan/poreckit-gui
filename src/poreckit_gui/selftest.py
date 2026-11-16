"""Exercise the actual packaged GUI and its persistent subprocesses."""
import json
import platform
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from .app import Browser
from .catalog import catalog
from .demo import create_demo


def run(report_path):
    report: dict[str, Any] = dict(ok=False, frozen=bool(getattr(sys, "frozen", False)), platform=platform.platform(), checks=[])
    browser: Browser | None = None
    config = tempfile.TemporaryDirectory(prefix="poreckit-selftest-config-")
    try:
        app = QApplication.instance() or QApplication([])
        browser = Browser(settings=QSettings(str(Path(config.name)/"settings.ini"), QSettings.Format.IniFormat))
        browser.auto_update.setChecked(False)
        browser.show()
        errors = []
        # Intentional test-hook patching of an existing method.
        browser.error = lambda value: errors.append(str(value))  # type: ignore[method-assign]
        browser.plot_workbench.error = lambda value: errors.append(str(value))  # type: ignore[method-assign]
        def wait():
            deadline = time.monotonic()+90
            def busy():
                return any((browser.busy(), browser.plot_workbench.busy(), browser.sequence_browser.busy(), browser.assembly_editor.busy()))
            while busy() and time.monotonic() < deadline:
                app.processEvents()
                time.sleep(.01)
            if busy() or errors:
                raise RuntimeError("; ".join(errors)+"\n"+browser.log_view.toPlainText()+"\n"+browser.plot_workbench.log.toPlainText())
            app.processEvents()
        with tempfile.TemporaryDirectory(prefix="poreckit-selftest-") as folder:
            browser.load_file(create_demo(folder, [("chr1", 1200000), ("chr2", 800000)]))
            wait()
            assert browser.image.map_data is not None and browser.image.map_data.observed.shape == (120, 120)
            report["checks"].append("native reader and GUI rendering")
            pid = browser.process.processId()
            signal = Path(folder)/"signal.bg"
            signal.write_text("chr1\t0\t1200000\t3\nchr2\t0\t800000\t7\n")
            browser.append_track(dict(kind="signal", path=str(signal), placement="auto"))
            browser.cmap.setCurrentText("viridis")
            browser.render()
            wait()
            assert browser.process.processId() == pid
            report["checks"].append("persistent worker reuse and signal track")
            browser.genome_dialog.open()
            wait()
            assert browser.genome_dialog.image.map_data is not None and browser.genome_dialog.image.map_data.observed.shape == (200, 200)
            report["checks"].append("whole-genome sparse aggregation")
            from .demo import create_track_demo
            sequence = browser.sequence_browser
            sequence.add_files(create_track_demo(Path(folder)/"tracks"))
            wait()
            from .model import Region
            sequence.navigate(Region("chr1", 90, 220))
            wait()
            assert not any("error" in t for t in sequence.result["tracks"]), sequence.result
            assert len(sequence.result["sequence"]) == 130
            assert next(t for t in sequence.result["tracks"] if t["track"]["kind"] == "alignment")["data"]["accepted"] == 4
            report["checks"].append("standalone signal, reference, reads and variant queries")
            editor = browser.assembly_editor
            editor.from_map()
            wait()
            assert editor.preview_result, editor.status.text()
            editor.document.invert(1, 0, 0)
            editor.refresh()
            wait()
            assert editor.document.lift(1200000)["position"] == 799999
            editor.document.mark_saved()
            report["checks"].append("assembly inversion and lifted contact rendering")
            entries = catalog()
            entry = next(e for e in entries.values() if e.title == "plot_aggregate_methods" and e.kind == "api")
            workbench = browser.plot_workbench
            workbench.select_entry(entry.key)
            workbench.run()
            wait()
            assert workbench.result["figures"]
            report["checks"].append("dynamic plot catalog and workbench subprocess")
            report["catalog_entries"] = len(entries)
            report["ok"] = True
    except Exception:
        report["error"] = traceback.format_exc()
    finally:
        if browser is not None:
            if browser.assembly_editor.document:
                browser.assembly_editor.document.mark_saved()
            browser.annotations.model.mark_saved()
            browser.bookmarks.mark_saved()
            browser.job = None
            browser.process.shutdown()
            browser.plot_workbench.shutdown()
            browser.close()
            browser.settings.sync()
        config.cleanup()
        Path(report_path).write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0 if report["ok"] else 1
