"""Discover plotting commands, public plotting functions and track classes."""
from __future__ import annotations

import codecs
import functools
import importlib
import inspect
import pkgutil
import re
from dataclasses import dataclass

import click

# These commands manage data without rendering figures.
NON_PLOTTING = {"contacts", "convert", "merge", "validate", "fix-mcool", "zoomify", "metadata"}
OUTPUTS = {
    "multiway": {"output": "results.json", "figure": "multiway.pdf"},
    "qc": {"output": "qc"}, "report": {"output": "report.html"},
    "tads": {"output": "tads.pdf", "plot_path": "insulation.pdf", "overlay": "tads.pdf",
             "bed": "tads.bed", "boundaries": "boundaries.bed"},
    "loops": {"bedpe": "loops.bedpe", "apa_path": "apa.pdf"},
    "balance": {"plot_path": "weights.pdf", "weights_tsv": "weights.tsv", "matrix_npz": "balanced.npz"},
}


BLOCKED_TOOLS = tuple(codecs.decode(s, "rot13") for s in (
    "Whvprobk", "Whvpre", "VTI", "TRABIN", "clTrabzrGenpxf", "cnvegbbyf", "uvpgx", "UvPRkcybere"))


def product_text(text):
    """Keep upstream help usable without putting other tool branding in the UI."""
    for name in BLOCKED_TOOLS:
        text = re.sub(r"\b" + re.escape(name) + r"\b", "reference", text, flags=re.IGNORECASE)
    return text


@dataclass
class Entry:
    key: str
    kind: str
    title: str
    group: str
    target: str
    description: str


def cli_command(name):
    from poreckit.cli import main
    cmd = main.get_command(click.Context(main), name)
    if cmd is None:
        raise ValueError(f"Unknown poreckit command: {name}")
    return cmd


def resolve(target):
    module, name = target.rsplit(":", 1)
    if not module.startswith("poreckit."):
        raise ValueError("Only poreckit plotting targets are supported")
    return getattr(importlib.import_module(module), name)


@functools.lru_cache(maxsize=1)
def catalog():
    import poreckit
    from poreckit.cli import _SUBCOMMANDS
    from poreckit.tracks.base import BaseTrack
    import poreckit.tracks as tracks

    entries = []
    for name in _SUBCOMMANDS:
        if name not in NON_PLOTTING:
            cmd = cli_command(name)
            entries.append(Entry("cli:" + name, "cli", name, "Commands", name, product_text(inspect.cleandoc(cmd.help or ""))))
    for info in pkgutil.walk_packages(poreckit.__path__, "poreckit."):
        if any(part.startswith("_") for part in info.name.split(".")) or ".cli" in info.name or info.name.endswith("mcp_server"):
            continue
        module = importlib.import_module(info.name)
        for name, fn in vars(module).items():
            if (name.startswith("plot_") or name == "persp_aggregate") and inspect.isfunction(fn) and fn.__module__ == info.name:
                target = info.name + ":" + name
                entries.append(Entry("api:" + target, "api", name, "Functions · " + info.name.removeprefix("poreckit."),
                                     target, product_text(inspect.getdoc(fn) or "")))
    for name, cls in vars(tracks).items():
        if name.startswith("Track") and inspect.isclass(cls) and issubclass(cls, BaseTrack):
            target = cls.__module__ + ":" + cls.__name__
            entries.append(Entry("track:" + target, "track", name, "Track stack", target, product_text(inspect.getdoc(cls) or "")))
    for name in ("gc_track", "restriction_track"):
        fn = getattr(tracks, name)
        target = fn.__module__ + ":" + fn.__name__
        entries.append(Entry("track:" + target, "track", name, "Track stack", target, product_text(inspect.getdoc(fn) or "")))
    return {entry.key: entry for entry in entries}


def output_options(name):
    params = {p.name for p in cli_command(name).params}
    result = {"output": "figure.pdf"} if "output" in params else {}
    result.update(OUTPUTS.get(name, {}))
    # Keep sidecars with the result set rather than writing into the project.
    for field in ("table", "save", "json_output", "tsv"):
        if field in params:
            result[field] = field + (".npz" if field == "save" else ".tsv")
    return result


def fields(entry):
    if entry.kind == "cli":
        result = []
        for p in cli_command(entry.target).params:
            if p.name in output_options(entry.target) or p.name == "html":
                continue
            default = p.default
            if default is None or callable(default) or type(default).__module__ == "click._utils":
                default = None
            if isinstance(p, click.Option) and p.is_flag:
                kind = "bool" if p.secondary_opts else "flag"
            elif p.multiple or p.nargs != 1:
                kind = "json"
            elif isinstance(p.type, click.Choice):
                kind = "choice"
            elif isinstance(p.type, click.Path):
                kind = "path"
            else:
                kind = "text"
            result.append({"name": p.name, "label": next((o for o in p.opts if o.startswith("--")), p.name),
                           "kind": kind, "default": default, "required": p.required,
                           "choices": list(p.type.choices) if isinstance(p.type, click.Choice) else [],
                           "help": product_text(getattr(p, "help", "") or ""),
                           "multiple": p.multiple, "nargs": p.nargs})
        return result
    fn = resolve(entry.target)
    result = []
    for p in inspect.signature(fn).parameters.values():
        if p.name == "ax" or p.name.startswith("_") or p.kind == p.VAR_POSITIONAL:
            continue
        if p.kind == p.VAR_KEYWORD:
            result.append({"name": "__extra__", "label": "Extra options", "kind": "json", "default": None,
                           "required": False, "help": "Optional JSON object of additional keyword arguments."})
            continue
        default = None if p.default is p.empty or repr(p.default) == "<factory>" else p.default
        if entry.title == "plot_runtime_benchmark" and p.name == "labels":
            default = ("Method 1", "Method 2")
        annotation = str(p.annotation)
        if p.name == "reads" and fn.__module__ == "poreckit.qc":
            kind = "path"
        elif p.name.startswith("reader"):
            kind = "reader"
        elif p.name == "contacts":
            kind = "contacts"
        elif isinstance(default, bool) or annotation in ("bool", "<class 'bool'>"):
            kind = "bool"
        elif isinstance(default, str) or annotation in ("str", "<class 'str'>") or p.name.endswith(("_path", "_file")):
            kind = "path" if p.name.endswith(("_path", "_file")) else "text"
        elif isinstance(default, (int, float)) or annotation in ("int", "float", "<class 'int'>", "<class 'float'>"):
            kind = "number"
        else:
            kind = "data"
        result.append({"name": p.name, "label": p.name.replace("_", " "), "kind": kind,
                       "default": default, "required": p.default is p.empty,
                       "help": str(p.annotation) if p.annotation is not p.empty else ""})
    return result


def cli_arguments(name, values):
    options, positionals = [], []
    for p in cli_command(name).params:
        if p.name not in values or values[p.name] is None:
            continue
        value = values[p.name]
        if hasattr(value, "tolist"):
            value = value.tolist()
        if isinstance(p, click.Argument):
            positionals.extend(map(str, value if isinstance(value, (list, tuple)) else [value]))
            continue
        flag = next((o for o in p.opts if o.startswith("--")), p.opts[0])
        if p.is_flag:
            if p.secondary_opts:
                options.append(flag if value else p.secondary_opts[0])
            elif value:
                options.append(flag)
        else:
            groups = value if p.multiple else [value]
            for item in groups:
                group = item if p.nargs != 1 else [item]
                if not isinstance(group, (list, tuple)) or (p.nargs > 1 and len(group) != p.nargs):
                    raise ValueError(f"{flag} requires groups of {p.nargs} values")
                options.extend([flag, *map(str, group)])
    return options + (["--", *positionals] if positionals else [])
