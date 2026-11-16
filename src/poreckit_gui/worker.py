"""Worker entry point shared by source and frozen desktop applications."""
from __future__ import annotations

import gzip
import json
import sys
import traceback
from pathlib import Path


def inspect_file(path):
    from poreckit import make_reader
    reader = make_reader(path)
    try:
        resolutions = [int(n) for n in reader.resolutions()]
        normalizations = {}
        assembly = None
        if path.lower().endswith(".hic"):
            meta = reader.meta()
            normalizations = meta.get("norms", {})
            assembly = meta.get("assembly")
        elif path.lower().endswith((".cool", ".mcool")):
            import h5py
            with h5py.File(path, "r") as source:
                assembly = source.attrs.get("genome-assembly", source.attrs.get("assembly"))
                groups = [(int(n), g) for n, g in source["resolutions"].items()] if "resolutions" in source else [(resolutions[0], source)]
                for resolution, group in groups:
                    if assembly is None:
                        assembly = group.attrs.get("genome-assembly", group.attrs.get("assembly"))
                    for name in group["bins"]:
                        if name not in {"chrom", "start", "end"} and group["bins"][name].ndim == 1:
                            normalizations.setdefault("balanced" if name == "weight" else name, []).append(resolution)
        elif path.lower().endswith(".h5"):
            import h5py
            with h5py.File(path, "r") as source:
                if "correction_factors" in source:
                    normalizations = {"balanced": resolutions}
        chroms = []
        if hasattr(reader, "chromsizes"):
            chroms = [(str(n), int(size)) for n, size in reader.chromsizes()]
        elif path.lower().endswith(".hic"):
            sizes = reader.meta()["chromsizes"]
            chroms = list(sizes.items()) if isinstance(sizes, dict) else list(sizes)
        elif path.lower().endswith((".pairs", ".pairs.gz")):
            opener = gzip.open if path.endswith(".gz") else open
            with opener(path, "rt") as source:
                for line in source:
                    if not line.startswith("#"):
                        break
                    if line.startswith("#chromsize:"):
                        _, chrom, length = line.split()
                        chroms.append((chrom, int(length)))
        if isinstance(assembly, bytes):
            assembly = assembly.decode()
        return {"input": path, "resolutions": resolutions, "normalizations": normalizations, "assembly": assembly,
                "chromsizes": [(str(n), int(s)) for n, s in chroms if int(s) > 0 and str(n).lower() != "all"]}
    finally:
        if hasattr(reader, "close"):
            reader.close()


def render(view, output, capture_formats=False, state=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import FuncFormatter
    from poreckit import (Canvas, TrackHeatmap, TrackGene,
                          TrackTAD, TrackLoop, TrackXAxis, make_reader)
    from .model import validate_view
    from .tracks import SignalPanel, SquareLoopOverlay
    from .axistracks import AxisTracks, placement
    from .viewport import Viewport
    from .contactview import ContactPanel, transform_matrix, EXPECTED_MODES, CONTROL_EXPECTED_MODES, SIGNED_MODES
    from .cache import source_chromosome

    region = validate_view(view)
    viewport = Viewport.from_view(view)
    region_y = viewport.y
    reader = state.reader(view["input"]) if state else make_reader(view["input"])
    secondary_readers = []
    try:
        resolution = view.get("binsize")
        if resolution is None:
            available = sorted(reader.resolutions())
            if view.get("control"):
                control_reader = state.reader(view["control"]) if state else make_reader(view["control"])
                secondary_readers.append(control_reader)
                common_levels = set(control_reader.resolutions())
                available = [n for n in available if n in common_levels] if common_levels else available
                if reader.resolutions() and common_levels and not available:
                    raise ValueError("Observed and control maps have no common resolution")
            target = max(1, (max(region.end-region.start, region_y.end-region_y.start) + 249) // 250)
            resolution = next((n for n in available if n >= target),
                              available[-1] if available else target)
        if max((r.end-r.start+resolution-1)//resolution for r in (region, region_y)) > 2500:
            raise ValueError("The file resolution produces more than 2500 bins. Select a smaller region.")
        def align(r):
            from .model import Region
            length = dict(reader.chromsizes()).get(source_chromosome(reader, r.chrom)) if hasattr(reader, "chromsizes") else None
            end = ((r.end+resolution-1)//resolution)*resolution
            return Region(r.chrom, (r.start//resolution)*resolution, min(end, length) if length else end)
        region, region_y = align(region), align(region_y)
        viewport = Viewport(region, region_y)
        if max((r.end-r.start+resolution-1)//resolution for r in (region, region_y)) > 2500:
            raise ValueError("Aligned region exceeds 2500 bins; select a coarser resolution")
        if view.get("control"):
            observed_meta = state.metadata(view["input"]) if state else inspect_file(view["input"])
            control_meta = state.metadata(view["control"]) if state else inspect_file(view["control"])
            def canonical(value):
                value = str(value or "").lower().replace("-", "").replace("_", "")
                return {"hg19": "grch37", "hs37d5": "grch37", "hg38": "grch38", "mm10": "grcm38", "mm39": "grcm39"}.get(value, value)
            a, b = canonical(observed_meta.get("assembly")), canonical(control_meta.get("assembly"))
            if a and b and a != b:
                raise ValueError("Observed and control maps use different genome assemblies")
        def extract(path, rd, x, y, norm, binsize=resolution):
            if state:
                if view.get("tile_cache") and str(path).lower().endswith((".hic", ".cool", ".mcool")):
                    return state.tiled_rectangle(path, x, y, binsize, norm)
                return state.rectangle(path, x, y, binsize, norm)
            if x == y:
                return rd.matrix(source_chromosome(rd, x.chrom), x.start, x.end, binsize, norm)
            return rd.matrix(source_chromosome(rd, y.chrom), y.start, y.end, binsize, norm,
                             chrom2=source_chromosome(rd, x.chrom), start2=x.start, end2=x.end)
        matrix = extract(view["input"], reader, region, region_y, view.get("norm"))
        control_matrix = None
        if view.get("control"):
            control_reader = state.reader(view["control"]) if state else make_reader(view["control"])
            secondary_readers.append(control_reader)
            control_matrix = extract(view["control"], control_reader, region, region_y, view.get("control_norm", view.get("norm")))
        display_mode = view.get("display_mode", "observed")
        if display_mode in {"split", "split_oe", "oe_local", "pearson_local"} and not viewport.diagonal:
            raise ValueError("The selected display needs identical X and Y regions")
        from .expected import expected_matrix
        def file_expected(path, norm):
            if not path or not path.lower().endswith(".hic"):
                raise ValueError("File expected displays require a HIC map containing expected curves; use visible-region O/E for other formats")
            return expected_matrix(path, region, region_y, resolution, norm)
        expected = file_expected(view["input"], view.get("norm")) if display_mode in EXPECTED_MODES else None
        control_expected = file_expected(view.get("control"), view.get("control_norm", view.get("norm"))) if display_mode in CONTROL_EXPECTED_MODES else None
        shown = transform_matrix(matrix, control_matrix, display_mode, float(view.get("pseudocount", 1)), expected, control_expected)
        tracks = [t for t in view.get("tracks", []) if t.get("visible", True)]
        square = view["mode"] == "square"
        axis_tracks = AxisTracks(tracks, square)
        below_count = sum(placement(t, square) == "panel" for t in tracks)
        map_count = sum(t["kind"] == "contacts" for t in tracks)
        canvas = Canvas(figsize=(10 + axis_tracks.left_inches, 9 + axis_tracks.top_inches + below_count * 1.5 + map_count * 3)
                        if square else (12, 5.4 + below_count + map_count * 3), dpi=140)
        if square:
            primary = ContactPanel(x=region, y=region_y, matrix=shown, name=Path(view["input"]).name,
                                   cmap=view.get("cmap", "hic_whitered"), log=view.get("log", False),
                                   vmax=view.get("vmax"), percentile=view.get("vmax_percentile"),
                                   signed=display_mode in SIGNED_MODES,
                                   display=display_mode, grid=view.get("grid", False), resolution=resolution)
        else:
            # Triangle path: TrackHeatmap only handles nonnegative matrices
            # (plain log1p, vmin=0), so signed displays pre-transform here and
            # draw linear limits over the signed-log values instead.
            signed = display_mode in SIGNED_MODES
            values = np.asarray(shown, dtype=np.float32)
            vmin, vmax, log = 0, view.get("vmax"), view.get("log", True)
            percentile = view.get("vmax_percentile")
            if signed:
                finite = np.abs(values[np.isfinite(values)])
                upper = vmax
                if upper is None:
                    upper = (float(np.percentile(finite, percentile)) if finite.size and percentile is not None
                             else float(finite.max()) if finite.size else 1.0)
                upper = max(upper, 1e-12)
                if log:
                    values = np.sign(values) * np.log1p(np.abs(values))
                    vmin, vmax = -np.log1p(upper), np.log1p(upper)
                else:
                    vmin, vmax = -upper, upper
                log = False
                percentile = None
            primary = TrackHeatmap(reader=reader, chrom=region.chrom,
                               start=region.start, end=region.end, bin_size=resolution,
                               norm=view.get("norm") or None, cmap=view.get("cmap", "Reds"),
                               log=log, square=view["mode"] == "square",
                               colorbar=True, max_pixels=800, aspect="equal" if square else "auto",
                               vmin=vmin, vmax=vmax, vmax_percentile=percentile,
                               matrix_override=values,
                               name=Path(view["input"]).name)
        canvas.add(primary)
        # TAD/loop overlays belong to the contact panel, before independent rows.
        def is_overlay(track):
            return track["kind"] in {"tads", "loops"} and placement(track, square) == "overlay"
        tracks = [t for t in tracks if is_overlay(t)] + [t for t in tracks if not is_overlay(t)]
        labels = [Path(view["input"]).name]
        panel_tracks = []
        for track in tracks:
            kind, path = track["kind"], track["path"]
            common = {"chrom": track.get("chrom") or region.chrom, "name": track.get("name") or Path(path).name}
            color = track.get("color", "#397ca8")
            if kind == "genes":
                gene_data = {"gtf_path": "", "genes": state.records("genes", path, common["chrom"])} if state else {"gtf_path": path}
                obj = TrackGene(**gene_data, style=track.get("style", "arrow"), maxgap=0,
                                colour=color, cds_colour=color, **common)
            elif kind == "signal":
                values = state.signal(path, common["chrom"], region.start, region.end, index=track.get("index")) if state else None
                obj = SignalPanel(path=path, mode=track.get("style", "fill"), color=color, values_override=values, **common)
            elif kind == "tads":
                if square and is_overlay(track) and region.chrom != region_y.chrom:
                    continue
                tad_data = {"tad_path": "", "tads": state.records("tads", path, common["chrom"])} if state else {"tad_path": path}
                obj = TrackTAD(**tad_data, bin_size=resolution,
                               colour=color, style="box" if square and is_overlay(track) else track.get("style", "triangle"), **common)
            elif kind == "loops":
                cls = SquareLoopOverlay if square and is_overlay(track) else TrackLoop
                style = track.get("style", "box" if square and is_overlay(track) else "arc")
                if cls is TrackLoop and style == "box":
                    style = "arc"
                loop_data = {"loop_path": "", "loops": state.records("loops", path, None if cls is SquareLoopOverlay else common["chrom"])} if state else {"loop_path": path}
                obj = cls(**loop_data, colour=color, style=style,
                          **({"y_chrom": track.get("chrom_y") or region_y.chrom} if cls is SquareLoopOverlay else {}), **common)
            elif kind == "contacts":
                panel_reader = state.reader(path) if state else make_reader(path)
                from .cache import source_chromosome
                common["chrom"] = source_chromosome(panel_reader, common["chrom"])
                secondary_readers.append(panel_reader)
                panel_res = track.get("binsize", "auto")
                if panel_res == "auto":
                    available = sorted(panel_reader.resolutions())
                    target = max(1, (region.end - region.start + 249) // 250)
                    panel_res = next((r for r in available if r >= target), available[-1] if available else target)
                panel_res = int(panel_res)
                if panel_res <= 0 or (region.end - region.start + panel_res - 1) // panel_res > 2500:
                    raise ValueError("Contact panel exceeds 2500 bins; choose a coarser resolution")
                panel_norm = track.get("norm", view.get("norm"))
                panel_norm = None if panel_norm == "raw" else panel_norm
                if track.get("style", "triangle") == "square":
                    from .model import Region
                    px = Region(common["chrom"], region.start, region.end)
                    py = Region(track.get("chrom_y") or region_y.chrom, region_y.start, region_y.end)
                    if max((r.end-r.start+panel_res-1)//panel_res for r in (px, py)) > 2500:
                        raise ValueError("Contact panel exceeds 2500 bins on Y; choose a coarser resolution")
                    qx = Region(px.chrom, px.start//panel_res*panel_res, px.end)
                    qy = Region(py.chrom, py.start//panel_res*panel_res, py.end)
                    panel_matrix = extract(path, panel_reader, qx, qy, panel_norm, panel_res)
                    obj = ContactPanel(x=px, y=py, matrix=panel_matrix, resolution=panel_res,
                        matrix_x_start=qx.start, matrix_y_start=qy.start,
                        name=common["name"], cmap=track.get("cmap") or view.get("cmap", "Reds"),
                        log=view.get("log", False), vmax=view.get("vmax"), percentile=view.get("vmax_percentile"))
                else:
                    panel_matrix = state.matrix(path, common["chrom"], region.start, region.end, panel_res, panel_norm) if state else None
                    obj = TrackHeatmap(reader=panel_reader, start=region.start, end=region.end, bin_size=panel_res,
                        norm=panel_norm, matrix_override=panel_matrix, square=False, aspect="auto",
                        cmap=track.get("cmap") or view.get("cmap", "Reds"), log=view.get("log", True),
                        vmin=0, vmax=view.get("vmax"), vmax_percentile=view.get("vmax_percentile"),
                        colorbar=True, max_pixels=800, **common)
                common["name"] += f" · {panel_res:,} bp / bin"
            else:
                raise ValueError(f"Unknown track type: {kind}")
            if placement(track, square) in {"top", "left", "both"}:
                other = obj
                if region_y != region:
                    from dataclasses import replace
                    y_chrom = track.get("chrom_y") or region_y.chrom
                    if kind == "signal":
                        values = state.signal(path, y_chrom, region_y.start, region_y.end, index=track.get("index")) if state else None
                        other = replace(obj, chrom=y_chrom, values_override=values)
                    elif kind == "genes" and y_chrom != common["chrom"]:
                        gene_data = {"gtf_path": "", "genes": state.records("genes", path, y_chrom)} if state else {"gtf_path": path}
                        other = TrackGene(**gene_data, chrom=y_chrom, name=common["name"], style=track.get("style", "arrow"),
                                          maxgap=0, colour=color, cds_colour=color)
                axis_tracks.register(track, obj, other)
                continue
            if not is_overlay(track):
                obj.overlay = False
                obj.height = float(track.get("height", 1.0 if kind == "contacts" else 0.3))
                labels.append(common["name"])
                panel_tracks.append(track)
            canvas.add(obj)
        if not square or panel_tracks:
            canvas.add(TrackXAxis())
        fig = canvas.render()
        axis_tracks.reserve(fig)
        fig.axes[0].get_subplotspec().get_gridspec().update(hspace=0.24)
        panel_axes = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
        heatmap_ax = panel_axes[0]
        from .annotations import layers_for_view, draw_layers, feature_payload
        layers = layers_for_view(view)
        draw_layers(heatmap_ax, layers, region, region_y, triangle=not square)
        heatmap_ax.patch.set_gid("contact-map-bounds")
        heatmap_ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e6:g}"))
        heatmap_ax.set_xlabel(f"{region.chrom} (Mb)")
        if square:
            heatmap_ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e6:g}"))
            heatmap_ax.set_ylabel(f"{region_y.chrom} (Mb)")
            if not panel_tracks:
                heatmap_ax.tick_params(axis="x", labelbottom=True)
        for ax, name in zip(panel_axes[1:], labels[1:]):
            ax.set_title(name, loc="left", fontsize=8)
        for ax, track in zip(panel_axes[1:], panel_tracks):
            if track["kind"] == "contacts":
                ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1e6:g}"))
                ax.set_ylabel(f"{track.get('chrom_y') or region_y.chrom} (Mb)" if track.get("style") == "square" else "Distance (Mb)")
            if track["kind"] in {"genes", "tads", "loops"}:
                ax.set_yticks([])
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.spines["left"].set_visible(False)
        if square or any(t["kind"] == "contacts" and t.get("style") == "square" for t in panel_tracks):
            # Equal aspect can shrink the main axes; keep annotation rulers aligned.
            fig.canvas.draw()
            main_position = heatmap_ax.get_position()
            width = min(ax.get_position().width for ax in panel_axes)
            left = main_position.x0 + (main_position.width - width) / 2
            for ax in panel_axes:
                position = ax.get_position()
                ax.set_position([left, position.y0, width, position.height])
        axis_tracks.draw(fig, heatmap_ax, region.chrom, region.start, region.end, region_y)
        fig.savefig(output, bbox_inches="tight", dpi=180 if Path(output).suffix.lower() == ".png" else 140)
        files = {Path(output).suffix.lower().lstrip("."): output}
        if capture_formats:
            for extension in ("svg", "pdf", "png"):
                if extension not in files:
                    path = str(Path(output).with_suffix("." + extension))
                    fig.savefig(path, bbox_inches="tight", dpi=180)
                    files[extension] = path
        plt.close(fig)
        limits = heatmap_ax.images[0].get_clim()
        if view.get("log", True):
            limits = [float(np.sign(v)*np.expm1(abs(v))) for v in limits]
        map_data = None
        if square and "svg" in files:
            from .mapdata import save_map_data, svg_plot_bounds
            map_data = save_map_data(output, matrix, shown, control_matrix, region, region_y, resolution)
            map_data.update(svg_plot_bounds(files["svg"]) or {})
            map_data["features"] = feature_payload(layers, region, region_y)
        return {"output": output, "files": files, "binsize": resolution, "region": region.text(), "region_y": region_y.text(),
                "map_data": map_data,
                "cache": {"hits": state.hits, "misses": state.misses, "bytes": state.bytes} if state else {},
                "display": {"mode": view["mode"], "cmap": view.get("cmap", "Reds"),
                            "log": view.get("log", True), "norm": view.get("norm") or "raw",
                            "color_min": float(limits[0]), "color_max": float(limits[1]), "display_mode": display_mode}}
    finally:
        plt.close("all")
        if state is None and hasattr(reader, "close"):
            reader.close()
        if state is None:
            for panel_reader in secondary_readers:
                if hasattr(panel_reader, "close"):
                    panel_reader.close()


def main(job_path, result_path, state=None):
    try:
        job = json.loads(Path(job_path).read_text(encoding="utf-8"))
        if job["action"] == "sequence_inspect":
            from .sequenceio import inspect_sequence
            import struct
            tracks = []
            chromsizes = job.get("chromsizes")
            for path in job["paths"]:
                try:
                    track = inspect_sequence(path, chromsizes, job.get("indexes", {}).get(path))
                    tracks.append(track)
                    if track["kind"] == "reference" or not chromsizes:
                        chromsizes = track["chromsizes"]
                except (ValueError, OSError, struct.error) as exc:
                    tracks.append({"path": path, "error": str(exc)})
            result = {"tracks": tracks}
        elif job["action"] == "sequence_query":
            from .sequenceio import query_sequence
            result = query_sequence(job, state)
        elif job["action"] == "sequence_index":
            from .sequenceio import tool_output
            list(tool_output("samtools", ["index", "-b", "-o", job["output"], job["path"]]))
            result = {"index": job["output"]}
        elif job["action"] == "assembly_preview":
            from .assembly import render_assembly
            result = render_assembly(job, state)
        elif job["action"] == "assembly_export":
            from .assembly import export_assembly
            result = export_assembly(job, state)
        elif job["action"] == "plot":
            from .execution import execute_plot
            result = execute_plot(job)
        elif job["action"] == "inspect":
            result = inspect_file(job["input"])
        elif job["action"] == "genome":
            from .genome import genome_data, draw_genome
            meta = job["metadata"]
            data, layout, resolution = genome_data(job["input"], meta["chromsizes"], meta["resolutions"], job.get("resolution"), state)
            result = draw_genome(data, layout, resolution, job["output"], log=job.get("log", False), name=Path(job["input"]).name)
        elif job["action"] == "genome_export":
            from .genome import draw_genome
            from .mapdata import MapData
            previous = job["result"]
            data = MapData(previous["map_data"]).observed
            result = draw_genome(data, previous["chromosomes"], previous["binsize"], job["output"],
                                 log=previous["log"], cmap=previous["cmap"], name=previous["name"])
        elif job["action"] == "gene_search":
            from .genesearch import search_genes
            result = search_genes(job["paths"], job["query"], job["chromsizes"])
        else:
            output = job["output"] if job.get("export") else str(Path(job["output"]).with_suffix(".svg"))
            result = render(job["view"], output, state=state)
        Path(result_path).write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
        return 0
    except Exception:
        traceback.print_exc()
        return 1


def _own_process_group():
    """Become a session leader so the GUI can signal the whole worker process
    tree (including samtools/bcftools children) on cancel."""
    import os
    if os.name == "posix":
        try:
            os.setsid()
        except OSError:
            pass  # already a group leader


def serve():
    import signal
    _own_process_group()
    def stop(signum, frame):
        # Convert SIGTERM into an exception so `finally` blocks run and the
        # external tools started by the current job are terminated.
        raise SystemExit(0)
    try:
        signal.signal(signal.SIGTERM, stop)
    except (ValueError, OSError):
        pass  # not the main thread, or unsupported platform
    from .cache import RenderCache
    state = RenderCache()
    try:
        for line in sys.stdin:
            try:
                request = json.loads(line)
                code = main(request["job"], request["result"], state=state)
            except Exception:
                traceback.print_exc()
                code = 1
            print("\n@poreckit-result:" + json.dumps({"code": code}), flush=True)
    finally:
        state.close()


if __name__ == "__main__":
    _own_process_group()
    if sys.argv[1:] == ["--serve"]:
        serve()
    else:
        raise SystemExit(main(*sys.argv[1:]))
