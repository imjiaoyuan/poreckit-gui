"""Regenerate the plot catalog from the installed poreckit package."""
from collections import defaultdict
from pathlib import Path

from poreckit_gui.catalog import catalog

groups = defaultdict(list)
for entry in catalog().values():
    groups[entry.group].append(entry)
lines = ["# Plot catalog", "", "Generated from the installed poreckit definitions by `python scripts/write_catalog.py`.",
         "All entries use the same zoom, pan, full-figure export and viewport export controls.", ""]
for group, entries in groups.items():
    lines.extend([f"## {group}", "", "| Plot / track | Entry key |", "| --- | --- |"])
    for entry in entries:
        lines.append(f"| {entry.title} | `{entry.key}` |")
    lines.append("")
path = Path(__file__).resolve().parents[1] / "docs" / "plot-catalog.md"
path.parent.mkdir(exist_ok=True)
path.write_text("\n".join(lines), encoding="utf-8")
print(path)
