"""Typed, inspectable forms for Click commands and Python plot signatures."""
from __future__ import annotations

import json

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from .catalog import fields


class ValueField(QWidget):
    def __init__(self, spec):
        super().__init__()
        self.spec = spec
        self.array = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        kind, default = spec["kind"], spec.get("default")
        if kind in {"bool", "flag", "choice"}:
            self.editor: QComboBox | QPlainTextEdit | QLineEdit = QComboBox()
            self.editor.addItem(f"Default ({default})" if default is not None else "Default", None)
            if kind == "bool":
                self.editor.addItem("True", True)
                self.editor.addItem("False", False)
            elif kind == "flag":
                self.editor.addItem("Apply flag", True)
            else:
                for choice in spec["choices"]:
                    self.editor.addItem(str(choice), choice)
            layout.addWidget(self.editor)
        else:
            structured = kind in {"json", "data"}
            self.editor = QPlainTextEdit() if structured else QLineEdit()
            if structured:
                self.editor.setFixedHeight(68)
            placeholder = json.dumps(default, ensure_ascii=False) if isinstance(default, (tuple, list, dict)) else str(default) if default is not None else ""
            if not placeholder:
                placeholder = "JSON value or data file" if structured else "Required" if spec["required"] else "Default"
            self.editor.setPlaceholderText(placeholder)
            layout.addWidget(self.editor)
            if kind in {"path", "reader", "contacts", "data", "json"}:
                row = QHBoxLayout()
                browse = QPushButton("Browse…")
                browse.clicked.connect(self.browse)
                row.addWidget(browse)
                if kind == "data":
                    self.array = QCheckBox("Array")
                    self.array.setToolTip("Convert inline JSON numeric values to a NumPy array")
                    row.addWidget(self.array)
                    matrix = QPushButton("Current matrix")
                    matrix.clicked.connect(self.current_matrix)
                    row.addWidget(matrix)
                row.addStretch()
                layout.addLayout(row)
        tooltip = spec.get("help", "")
        if kind == "json":
            tooltip += "\nUse a JSON array for repeated/tuple arguments; nested arrays for repeated tuples."
        elif kind == "data":
            tooltip += '\nJSON, NPY/NPZ or numeric CSV/TSV. Nested arrays: {"$array": [...]}. Files: {"$file": "path"}.'
        self.setToolTip(tooltip.strip())

    def current_matrix(self):
        entry = {"$matrix": True}
        value: object = [entry] if self.spec["name"] in {"matrices", "stacks", "aggregates"} else entry
        self.set_value(value)

    def browse(self):
        kind = self.spec["kind"]
        structured = kind in {"data", "json"}
        path, _ = QFileDialog.getOpenFileName(self, "Select data", "",
            "Data (*.json *.npy *.npz *.csv *.tsv *.txt);;All files (*)" if structured else "All files (*)",
            options=QFileDialog.Option.DontUseNativeDialog)
        if path:
            self.set_value({"$file": path} if structured else path)

    def value(self):
        kind = self.spec["kind"]
        if isinstance(self.editor, QComboBox):
            return self.editor.currentData()
        raw = (self.editor.toPlainText() if isinstance(self.editor, QPlainTextEdit) else self.editor.text()).strip()
        if not raw:
            return None
        try:
            if kind in {"json", "data", "number"}:
                value = json.loads(raw)
                if kind == "number" and (isinstance(value, bool) or not isinstance(value, (int, float))):
                    raise ValueError("Expected a number")
                if kind == "data" and self.array is not None and self.array.isChecked():
                    value = {"$array": value}
                return value
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{self.spec['label']}: {exc}") from exc
        if kind == "reader":
            return {"$reader": raw}
        if kind == "contacts":
            return {"$contacts": raw}
        return raw

    def set_value(self, value):
        if isinstance(self.editor, QComboBox):
            index = self.editor.findData(value)
            self.editor.setCurrentIndex(max(0, index))
            return
        if self.spec["kind"] in {"reader", "contacts"} and isinstance(value, dict):
            value = next(iter(value.values()))
        if self.spec["kind"] in {"json", "data"}:
            text = json.dumps(value, ensure_ascii=False)
        else:
            text = str(value)
        if isinstance(self.editor, QPlainTextEdit):
            self.editor.setPlainText(text)
        else:
            self.editor.setText(text)


class ParameterForm(QWidget):
    def __init__(self, entry):
        super().__init__()
        self.entry = entry
        self.inputs = {}
        self.optional = []
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        self.show_optional = QCheckBox("Show all optional parameters")
        self.show_optional.toggled.connect(self.toggle_optional)
        outer.addWidget(self.show_optional)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        container = QWidget()
        form = QFormLayout(container)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        basic = {"chrom", "start", "end", "binsize", "bin_size", "chrom2", "start2", "end2",
                 "norm", "cmap", "title", "name", "reader", "reader1", "reader2", "matrix_override",
                 "anchor", "viewpoint", "analysis", "style", "mode", "colour", "color", "overlay"}
        for spec in fields(entry):
            field = ValueField(spec)
            label = QLabel(spec["label"] + (" *" if spec["required"] else ""))
            label.setTextFormat(Qt.TextFormat.PlainText)
            label.setWordWrap(True)
            label.setToolTip(field.toolTip())
            form.addRow(label, field)
            self.inputs[spec["name"]] = field
            if not spec["required"] and spec["name"] not in basic:
                self.optional.append((label, field))
        outer.addWidget(scroll, 1)
        scroll.setWidget(container)
        self.toggle_optional(False)

    def toggle_optional(self, checked):
        for label, field in self.optional:
            label.setVisible(checked)
            field.setVisible(checked)

    def values(self):
        values = {}
        for name, field in self.inputs.items():
            value = field.value()
            if value is not None:
                values[name] = value
            elif field.spec["required"]:
                raise ValueError(f"Enter {field.spec['label']}")
        return values

    def set_values(self, values):
        for name, value in values.items():
            if name in self.inputs:
                self.inputs[name].set_value(value)
        if any(field.spec["name"] in values for _, field in self.optional):
            self.show_optional.setChecked(True)
