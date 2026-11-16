"""Run every plot through its original API and collect portable figure assets."""
from __future__ import annotations

import contextlib
import os
import shutil
from pathlib import Path
from unittest.mock import patch


class FigureCollector:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.figures = []
        self.assets = []

    def capture(self, fig, source=None):
        from matplotlib.figure import Figure
        if not isinstance(fig, Figure) or not fig.axes:
            return
        if fig in self.figures:
            return
        self.figures.append(fig)
        number = len(self.figures)
        title = fig._suptitle.get_text() if fig._suptitle is not None else ""
        if not title:
            title = next((ax.get_title() or ax.get_title(loc="left") for ax in fig.axes
                          if ax.get_title() or ax.get_title(loc="left")), f"Figure {number}")
        files = {}
        for extension in ("png", "svg", "pdf"):
            path = self.directory / f"plot-{number:03d}.{extension}"
            if isinstance(source, (str, os.PathLike)) and Path(source).suffix.lower() == "." + extension:
                if Path(source).resolve() != path.resolve():
                    shutil.copyfile(source, path)
            else:
                self.savefig(fig, path, format=extension, bbox_inches="tight", dpi=160)
            files[extension] = str(path)
        self.assets.append({"title": str(title), "files": files})

    @contextlib.contextmanager
    def intercept(self):
        import matplotlib.pyplot as plt
        from matplotlib.figure import Figure
        self.savefig = Figure.savefig
        original_close = plt.close

        def save(fig, *args, **kwargs):
            result = self.savefig(fig, *args, **kwargs)
            self.capture(fig, args[0] if args else kwargs.get("fname"))
            return result

        def close(fig=None):
            if fig == "all":
                for num in plt.get_fignums():
                    self.capture(plt.figure(num))
            elif isinstance(fig, Figure):
                self.capture(fig)
            elif plt.get_fignums():
                self.capture(plt.gcf() if fig is None else plt.figure(fig))
            return original_close(fig)

        with patch.object(Figure, "savefig", save), patch.object(plt, "close", close):
            yield
            for num in plt.get_fignums():
                self.capture(plt.figure(num))


def _kwargs(values, context, resources):
    from .plotdata import load_value
    result = {k: load_value(v, context, resources) for k, v in values.items() if k != "__extra__"}
    extra = load_value(values.get("__extra__", {}), context, resources)
    if not isinstance(extra, dict):
        raise ValueError("Extra options must be a JSON object")
    if "ax" in extra or any(k.startswith("_") for k in extra):
        raise ValueError("Axes and private parameters are managed by the viewer")
    result.update(extra)
    return result


def _guard_region(values):
    for a, b in (("start", "end"), ("start2", "end2")):
        if values.get(b) is None:
            continue
        start, end = int(values.get(a) or 0), int(values[b])
        if start < 0 or end <= start:
            raise ValueError("Genomic regions must satisfy 0 ≤ start < end")
        size = values.get("binsize", values.get("bin_size"))
        if size is not None and str(size) != "auto":
            size = int(size)
            if size <= 0 or (end - start + size - 1) // size > 2500:
                raise ValueError("Region exceeds 2500 bins; select a coarser resolution or smaller region")


def execute_plot(job):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from poreckit import Canvas
    from .catalog import catalog, cli_command, cli_arguments, output_options, resolve

    directory = Path(job["directory"]).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    collector = FigureCollector(directory)
    context = job.get("context", {})
    previous = Path.cwd()
    try:
        os.chdir(job.get("cwd", directory))
        with matplotlib.rc_context(), contextlib.ExitStack() as resources, collector.intercept():
            if job["kind"] == "track":
                region = context.get("region")
                from .model import Region
                if not region:
                    raise ValueError("Set a genomic locus for the track stack")
                r = Region.parse(region)
                canvas = Canvas(figsize=tuple(job.get("figsize", [12, 8])), dpi=150)
                canvas._chrom, canvas._start, canvas._end = r.chrom, r.start, r.end
                for track in job["tracks"]:
                    if not track.get("visible", True):
                        continue
                    entry = catalog()[track["entry"]]
                    if entry.kind != "track":
                        raise ValueError("Track stacks only accept track entries")
                    values = _kwargs(track["values"], context, resources)
                    _guard_region(values)
                    canvas.add(resolve(entry.target)(**values))
                if not canvas.tracks:
                    raise ValueError("Add at least one visible track to the stack")
                collector.capture(canvas.render())
            else:
                entry = catalog()[job["entry"]]
                values = dict(job.get("values", {}))
                if entry.kind == "cli":
                    from .plotdata import load_value
                    values = {k: load_value(v, context, resources) for k, v in values.items()}
                    # Outputs are generated automatically; export copies them later.
                    values.update({k: str(directory / v) for k, v in output_options(entry.target).items()})
                    values.pop("html", None)
                    _guard_region(values)
                    cli_command(entry.target).main(args=cli_arguments(entry.target, values), standalone_mode=False)
                elif entry.kind == "api":
                    kwargs = _kwargs(values, context, resources)
                    if entry.title == "plot_runtime_benchmark":
                        kwargs.setdefault("labels", ("Method 1", "Method 2"))
                    _guard_region(kwargs)
                    result = resolve(entry.target)(**kwargs)
                    _collect_return(result, collector)
                else:
                    raise ValueError("Use a track stack to render this entry")
        if not collector.assets:
            raise ValueError("No figures were produced. Check the selected parameters and input data.")
        figures = collector.assets
        captured = {Path(path) for item in figures for path in item["files"].values()}
        sidecars = [str(p) for p in directory.rglob("*") if p.is_file() and p not in captured]
        return {"figures": figures, "sidecars": sidecars}
    finally:
        os.chdir(previous)
        plt.close("all")


def _collect_return(result, collector):
    from matplotlib.figure import Figure
    if isinstance(result, Figure):
        collector.capture(result)
    elif hasattr(result, "figure") and isinstance(result.figure, Figure):
        collector.capture(result.figure)
    elif isinstance(result, (tuple, list)):
        for item in result:
            _collect_return(item, collector)
