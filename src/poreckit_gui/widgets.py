from typing import TYPE_CHECKING

from PySide6.QtCore import QRect, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen, QPixmap
from PySide6.QtSvgWidgets import QGraphicsSvgItem
from PySide6.QtWidgets import QGraphicsScene, QGraphicsView, QRubberBand, QWidget

if TYPE_CHECKING:
    from .mapdata import MapData



class ImageView(QGraphicsView):
    zoomChanged = Signal(float)
    dataZoom = Signal(float, float, float)
    dataPan = Signal(float, float)
    dataSelection = Signal(float, float, float, float)
    binHovered = Signal(object)
    binActivated = Signal(object)
    annotationSelected = Signal(float, float, float, float)
    featureSelected = Signal(object)
    featureEdited = Signal(object, float, float, bool)

    def __init__(self):
        super().__init__()
        self.setScene(QGraphicsScene(self))
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        self.setMinimumSize(400, 280)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.fit_mode = True
        self.fit_scale = 1.0
        self.map_data: MapData | None = None
        self.data_rect = QRectF()
        self.data_navigation = True
        self.crosshair = False
        self.annotation_mode = False
        self.cross_lines = []
        self.press = None
        self.feature_press = None
        self.rubber = QRubberBand(QRubberBand.Shape.Rectangle, self.viewport())
        self.setMouseTracking(True)
        self.setToolTip("Scroll to zoom around the pointer; drag to pan. Use Fit to reset the view.")

    def show_image(self, path):
        self.map_data = None
        self.data_rect = QRectF()
        self.cross_lines.clear()
        self.press = None
        self.feature_press = None
        self.rubber.hide()
        if str(path).lower().endswith(".svg"):
            item = QGraphicsSvgItem(str(path))
            if not item.renderer().isValid():
                raise ValueError("Unable to read the SVG preview")
            self.scene().clear()
            self.scene().addItem(item)
            self.scene().setSceneRect(item.boundingRect())
            self.fit()
            return
        pixmap = QPixmap(str(path))
        if pixmap.isNull():
            raise ValueError("Unable to read the preview image")
        self.scene().clear()
        self.scene().addPixmap(pixmap)
        self.scene().setSceneRect(QRectF(pixmap.rect()))
        self.fit()

    def clear_image(self):
        self.map_data = None
        self.data_rect = QRectF()
        self.cross_lines.clear()
        self.press = None
        self.feature_press = None
        self.rubber.hide()
        self.scene().clear()

    def set_map_data(self, spec):
        from .mapdata import MapData
        if not spec or not spec.get("bounds"):
            self.map_data = None
            return
        self.map_data = MapData(spec)
        width, height = spec["svg_size"]
        scene = self.scene().sceneRect()
        x, y, w, h = spec["bounds"]
        self.data_rect = QRectF(x*scene.width()/width, y*scene.height()/height, w*scene.width()/width, h*scene.height()/height)

    def fractions(self, point):
        if self.map_data is None or self.data_rect.isEmpty():
            return None
        scene = self.mapToScene(point)
        if not self.data_rect.contains(scene):
            return None
        return ((scene.x()-self.data_rect.x()) / self.data_rect.width(), (scene.y()-self.data_rect.y()) / self.data_rect.height())

    def mousePressEvent(self, event):
        fraction = self.fractions(event.position().toPoint())
        if fraction and event.button() == Qt.MouseButton.LeftButton and event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            hit = self.feature_at(*fraction)
            if hit:
                self.feature_press = (hit, fraction, bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier))
                event.accept()
                return
        if (self.data_navigation or self.annotation_mode) and fraction and event.button() == Qt.MouseButton.LeftButton:
            point = event.position().toPoint()
            selecting = self.annotation_mode or bool(event.modifiers() & Qt.KeyboardModifier.AltModifier)
            self.press = (point, fraction, selecting, self.annotation_mode)
            if selecting:
                self.rubber.setGeometry(QRect(point, point))
                self.rubber.show()
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        point = event.position().toPoint()
        fraction = self.fractions(point)
        if fraction and self.map_data is not None:
            self.binHovered.emit(self.map_data.at(*fraction))
        if self.press and self.press[2]:
            self.rubber.setGeometry(QRect(self.press[0], point).normalized())
        if self.crosshair and fraction:
            if not self.cross_lines:
                pen = QPen(QColor("#555555"), 0, Qt.PenStyle.DashLine)
                self.cross_lines = [self.scene().addLine(0, 0, 0, 0, pen) for _ in range(2)]
                for line in self.cross_lines:
                    line.setZValue(100)
            position = self.mapToScene(point)
            self.cross_lines[0].setLine(self.data_rect.left(), position.y(), self.data_rect.right(), position.y())
            self.cross_lines[1].setLine(position.x(), self.data_rect.top(), position.x(), self.data_rect.bottom())
        for line in self.cross_lines:
            line.setVisible(self.crosshair and fraction is not None)
        if self.press or self.feature_press:
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.feature_press and event.button() == Qt.MouseButton.LeftButton:
            hit, start, resize = self.feature_press
            self.feature_press = None
            end = self.fractions(event.position().toPoint())
            if end and max(abs(start[0]-end[0]), abs(start[1]-end[1])) > .001:
                self.featureEdited.emit(hit, end[0]-start[0], end[1]-start[1], resize)
            self.featureSelected.emit(hit)
            event.accept()
            return
        if self.press and event.button() == Qt.MouseButton.LeftButton:
            start, fraction, selecting, annotating = self.press
            self.press = None
            self.rubber.hide()
            scene = self.mapToScene(event.position().toPoint())
            end = (max(0, min(1, (scene.x()-self.data_rect.x())/self.data_rect.width())),
                   max(0, min(1, (scene.y()-self.data_rect.y())/self.data_rect.height())))
            if (event.position().toPoint()-start).manhattanLength() >= 6:
                if selecting:
                    (self.annotationSelected if annotating else self.dataSelection).emit(fraction[0], end[0], fraction[1], end[1])
                else:
                    self.dataPan.emit(fraction[0]-end[0], fraction[1]-end[1])
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event):
        fraction = self.fractions(event.position().toPoint())
        if fraction and self.map_data is not None and self.map_data.spec.get("chromosomes"):
            self.binActivated.emit(self.map_data.at(*fraction))
            event.accept()
            return
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            event.accept()
            return
        if self.data_navigation and fraction:
            self.press = None
            self.dataZoom.emit(.5, *fraction)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def feature_at(self, fx, fy):
        if self.map_data is None:
            return None
        spec = self.map_data.spec
        x = spec["x"][1] + fx*(spec["x"][2]-spec["x"][1])
        y = spec["y"][1] + fy*(spec["y"][2]-spec["y"][1])
        tx = (spec["x"][2]-spec["x"][1])*4/max(1, self.data_rect.width()*self.transform().m11())
        ty = (spec["y"][2]-spec["y"][1])*4/max(1, self.data_rect.height()*self.transform().m22())
        for item in reversed(spec.get("features", [])):
            x0, x1, y0, y1 = item["coords"]
            if x0-tx <= x <= x1+tx and y0-ty <= y <= y1+ty:
                return item
        return None

    def fit(self):
        if not self.scene().sceneRect().isEmpty():
            self.fit_mode = True
            self.fitInView(self.scene().sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
            self.fit_scale = self.transform().m11()
            self.zoomChanged.emit(1.0)

    def zoom_by(self, factor):
        if self.scene().sceneRect().isEmpty():
            return
        current = self.transform().m11() / self.fit_scale
        target = min(64, max(0.1, current * factor))
        self.fit_mode = False
        self.scale(target / current, target / current)
        self.zoomChanged.emit(target)

    def actual_size(self):
        if not self.scene().sceneRect().isEmpty():
            self.fit_mode = False
            self.resetTransform()
            self.zoomChanged.emit(1 / self.fit_scale)

    def wheelEvent(self, event):
        delta = event.angleDelta().y() or event.pixelDelta().y()
        fraction = self.fractions(event.position().toPoint())
        if delta and fraction and self.data_navigation and not event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.dataZoom.emit(.8 if delta > 0 else 1.25, *fraction)
            event.accept()
            return
        if delta and not self.scene().sceneRect().isEmpty():
            self.zoom_by(1.25 if delta > 0 else 0.8)
            event.accept()
        else:
            super().wheelEvent(event)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.fit_mode:
            self.fit()


class ChromosomeOverview(QWidget):
    navigate = Signal(float)

    def __init__(self):
        super().__init__()
        self.region = None
        self.length = None
        self.setMinimumHeight(60)
        self.setToolTip("Chromosome overview; click to navigate. Coordinates are zero-based, half-open.")

    def update_region(self, region, length):
        self.region, self.length = region, length
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(self.palette().text().color())
        if self.region is None:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Load contact data to browse genomic regions")
            return
        painter.drawText(16, 18, self.region.chrom)
        painter.drawText(100, 18, f"{self.region.end - self.region.start:,} bp")
        width = self.width() - 32
        length = self.length or self.region.end
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(self.palette().mid())
        painter.drawRoundedRect(QRectF(16, 30, width, 13), 5, 5)
        painter.setBrush(self.palette().highlight())
        x = 16 + width * self.region.start / length
        span = max(3, width * (self.region.end - self.region.start) / length)
        painter.drawRoundedRect(QRectF(x, 28, span, 17), 3, 3)

    def mousePressEvent(self, event):
        if self.region and self.length and event.button() == Qt.MouseButton.LeftButton:
            self.navigate.emit(max(0, min(1, (event.position().x() - 16) / max(1, self.width() - 32))))
