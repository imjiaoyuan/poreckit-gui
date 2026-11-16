"""Editable BEDPE layers with reversible operations and portable state."""
import copy
import csv
import functools
import json
import uuid
from pathlib import Path


def feature(value):
    value = dict(value)
    for key in ("x0", "x1", "y0", "y1"):
        value[key] = int(value[key])
    for axis in ("x", "y"):
        if not str(value[axis + "chrom"]).strip() or not 0 <= value[axis+"0"] < value[axis+"1"]:
            raise ValueError("Annotation intervals must satisfy 0 ≤ start < end")
    value.setdefault("id", uuid.uuid4().hex)
    value.setdefault("name", "")
    return value


def read_bedpe(path):
    rows = []
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            cells = line.rstrip().split()
            if len(cells) < 6:
                raise ValueError("BEDPE annotations require at least six columns")
            rows.append(feature(dict(xchrom=cells[0], x0=cells[1], x1=cells[2],
                                     ychrom=cells[3], y0=cells[4], y1=cells[5], name=cells[6] if len(cells)>6 else "", extra=cells[7:])))
    return rows


def write_bedpe(path, layer):
    with open(path, "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["#chrom1", "start1", "end1", "chrom2", "start2", "end2", "name"])
        for item in layer["features"]:
            writer.writerow([item[k] for k in ("xchrom", "x0", "x1", "ychrom", "y0", "y1", "name")] + item.get("extra", []))


class AnnotationModel:
    def __init__(self, on_change=lambda: None):
        self.layers = []
        self.undo_stack = []
        self.redo_stack = []
        self.saved = {}
        self.revision = 0
        self.on_change = on_change

    def layer(self, key):
        return next(value for value in self.layers if value["id"] == key)

    def _signature(self, value):
        return json.dumps({k: v for k, v in value.items() if k != "visible"}, sort_keys=True)

    @property
    def dirty(self):
        return {v["id"]: self._signature(v) for v in self.layers} != self.saved

    def mark_saved(self, key=None):
        if key is None:
            self.saved = {v["id"]: self._signature(v) for v in self.layers}
        else:
            self.saved[key] = self._signature(self.layer(key))
        self.on_change()

    def _notify(self):
        self.revision += 1
        self.on_change()

    def _execute(self, do, undo):
        do()
        self.undo_stack.append((do, undo))
        self.undo_stack = self.undo_stack[-100:]
        self.redo_stack.clear()
        self._notify()

    def add_layer(self, name="Annotations", values=(), color="#299fca"):
        layer = dict(id=uuid.uuid4().hex, name=name, color=color, visible=True, mirror=True,
                     features=[feature(v) for v in values])
        self._execute(lambda: self.layers.append(layer), lambda: self.layers.remove(layer))
        return layer["id"]

    def remove_layer(self, key):
        layer = self.layer(key)
        index = self.layers.index(layer)
        self._execute(lambda: self.layers.remove(layer), lambda: self.layers.insert(index, layer))

    def set_layer(self, key, **changes):
        layer = self.layer(key)
        before = {k: layer[k] for k in changes}
        self._execute(lambda: layer.update(changes), lambda: layer.update(before))

    def add(self, key, value):
        layer, value = self.layer(key), feature(value)
        self._execute(lambda: layer["features"].append(value), lambda: layer["features"].remove(value))
        return value["id"]

    def update(self, key, feature_id, **changes):
        target = next(v for v in self.layer(key)["features"] if v["id"] == feature_id)
        before = dict(target)
        after = feature({**target, **changes})
        self._execute(lambda: target.update(after), lambda: target.update(before))

    def remove(self, key, feature_id):
        values = self.layer(key)["features"]
        target = next(v for v in values if v["id"] == feature_id)
        index = values.index(target)
        self._execute(lambda: values.remove(target), lambda: values.insert(index, target))

    def undo(self):
        if self.undo_stack:
            item = self.undo_stack.pop()
            item[1]()
            self.redo_stack.append(item)
            self._notify()

    def redo(self):
        if self.redo_stack:
            item = self.redo_stack.pop()
            item[0]()
            self.undo_stack.append(item)
            self._notify()

    def restore(self, layers):
        new = copy.deepcopy(layers)
        for layer in new:
            layer.setdefault("id", uuid.uuid4().hex)
            layer.setdefault("visible", True)
            layer.setdefault("mirror", True)
            layer.setdefault("color", "#299fca")
            layer["features"] = [feature(v) for v in layer.get("features", [])]
        self.layers = new
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.mark_saved()
        self._notify()


@functools.lru_cache(maxsize=4)
def _layers_from_file(path, mtime, size):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def layers_for_view(view):
    if view.get("layers_file"):
        path = Path(view["layers_file"])
        stat = path.stat()
        return _layers_from_file(str(path), stat.st_mtime_ns, stat.st_size)
    return view.get("layers", [])


def visible_features(layer, x, y):
    def same(a, b):
        return a.removeprefix("chr") == b.removeprefix("chr")
    for item in layer["features"]:
        variants = [(item["xchrom"], item["x0"], item["x1"], item["ychrom"], item["y0"], item["y1"])]
        if layer.get("mirror", True):
            variants.append((item["ychrom"], item["y0"], item["y1"], item["xchrom"], item["x0"], item["x1"]))
        seen = set()
        for cx, lo1, hi1, cy, lo2, hi2 in variants:
            coords = lo1, hi1, lo2, hi2
            if coords in seen or not same(cx, x.chrom) or not same(cy, y.chrom):
                continue
            if hi1 <= x.start or lo1 >= x.end or hi2 <= y.start or lo2 >= y.end:
                continue
            seen.add(coords)
            yield item, coords


def draw_layers(ax, layers, x, y, triangle=False):
    from matplotlib.patches import Rectangle, Polygon
    for layer in layers:
        if not layer.get("visible", True):
            continue
        for item, (x0, x1, y0, y1) in visible_features(layer, x, y):
            if triangle:
                points = [((a+b)/2, abs(a-b)) for a,b in ((x0,y0),(x1,y0),(x1,y1),(x0,y1))]
                artist = Polygon(points, closed=True, fill=False)
            else:
                artist = Rectangle((x0, y0), x1-x0, y1-y0, fill=False)
            artist.set_edgecolor(layer["color"])
            artist.set_linewidth(float(layer.get("width", 1)))
            artist.set_zorder(7)
            ax.add_patch(artist)


def feature_payload(layers, x, y):
    result = []
    for layer in layers:
        if not layer.get("visible", True):
            continue
        for item, coords in visible_features(layer, x, y):
            original = (item["x0"], item["x1"], item["y0"], item["y1"])
            reflected = coords != original or item["xchrom"].removeprefix("chr") != x.chrom.removeprefix("chr")
            result.append({"layer": layer["id"], "id": item["id"], "coords": coords, "reflected": reflected})
    return result
