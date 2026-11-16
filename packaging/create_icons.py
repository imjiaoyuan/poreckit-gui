"""Rasterize the source SVG into native installer and application icons."""
import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtGui import QGuiApplication, QImage, QPainter
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtCore import Qt
from PIL import Image

root = Path(__file__).resolve().parent
app = QGuiApplication.instance() or QGuiApplication([])
image = QImage(1024, 1024, QImage.Format.Format_ARGB32)
image.fill(Qt.GlobalColor.transparent)
painter = QPainter(image)
QSvgRenderer(str(root/"poreckit.svg")).render(painter)
painter.end()
image.save(str(root/"poreckit.png"))
with Image.open(root/"poreckit.png") as source:
    source.save(root/"poreckit.ico", sizes=[(n, n) for n in (16, 24, 32, 48, 64, 128, 256)])
    source.save(root/"poreckit.icns")
