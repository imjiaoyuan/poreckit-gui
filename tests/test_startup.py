import os
import subprocess
import sys


def test_application_entry_point_starts_and_closes():
    script = '''
import runpy
import sys
from PySide6.QtCore import QTimer, QLocale
from poreckit_gui.app import Browser

original_show = Browser.show
def show_and_close(self):
    original_show(self)
    assert QLocale().name() == "en_US"
    QTimer.singleShot(50, self.close)

Browser.show = show_and_close
sys.argv = ["poreckit-gui"]
runpy.run_module("poreckit_gui", run_name="__main__")
'''
    result = subprocess.run([sys.executable, "-c", script], env={**os.environ, "QT_QPA_PLATFORM": "offscreen"},
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
