import sys
import os
from pathlib import Path
from PyInstaller.utils.hooks import collect_all, copy_metadata

datas, binaries, hiddenimports = [], [], []
for package in ("poreckit", "poreckit_gui"):
    data, binary, hidden = collect_all(package, include_py_files=False)
    datas += data
    binaries += binary
    hiddenimports += hidden
for package in ("poreckit-gui", "poreckit", "numpy", "matplotlib", "h5py", "PySide6-Essentials"):
    datas += copy_metadata(package)

a = Analysis([str(SPECPATH) + "/launch.py"], pathex=[], binaries=binaries,
             datas=datas, hiddenimports=hiddenimports,
             hooksconfig={"matplotlib": {"backends": ["Agg", "svg", "pdf"]}},
             excludes=["PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtQml", "tkinter", "pytest", "IPython", "PyQt5", "PyQt6", "PySide2"])
pyz = PYZ(a.pure)
icon = str(Path(SPECPATH) / "poreckit.ico") if sys.platform == "win32" else None
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="poreckit-gui", console=False,
          icon=icon, upx=False, codesign_identity=os.environ.get("MACOS_CODESIGN_IDENTITY"))
bundle = COLLECT(exe, a.binaries, a.datas, name="poreckit-gui", upx=False)
if sys.platform == "darwin":
    app = BUNDLE(bundle, name="poreckit-gui.app", bundle_identifier="org.poreckit.browser",
                 codesign_identity=os.environ.get("MACOS_CODESIGN_IDENTITY"),
                 icon=str(Path(SPECPATH)/"poreckit.icns"),
                 info_plist={"NSHighResolutionCapable": True, "CFBundleShortVersionString": os.environ.get("PORECKIT_VERSION", "0.1.0")})
