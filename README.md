# poreckit GUI

A native **PySide6 / Qt Widgets** genome browser for the sibling
[`poreckit`](../poreckit) library. One Python UI for Windows, Linux and macOS.
It includes a focused genome browser, a track/alignment viewer, an assembly editor and a full poreckit plotting workbench.

## Run

Keep this project alongside `poreckit`:

```text
work/
  poreckit/
  poreckit_gui/
```

Python 3.11+ and a Rust toolchain are needed to install the local poreckit
source. `uv` creates the virtual environment in this directory:

```sh
uv sync
uv run poreckit-gui
```

On Windows, use a Rust MSVC toolchain and Visual Studio C++ build tools.
Qt and the scientific Python dependencies are installed by `uv sync`.
No separate Qt installation is needed.

## Integration with poreckit

`uv` builds the sibling source checkout and installs poreckit into this project's
`.venv`. The GUI imports that installed package: the browser uses the public
reader/Canvas/track APIs, and the workbench calls the original Click command
objects or Python plotting functions. Native parsing and computation still run
in poreckit's Rust extension. No shell command strings are used for plotting.

This is a regular installation, not an editable source import. After changing
the sibling poreckit source, rebuild its installation here:

```sh
uv sync --reinstall-package poreckit
```

## Browse

- Open or drop a `.cool`, `.mcool`, `.hic`, `.h5`, `.pairs[.gz]` or `.parquet` file.
- Use the chromosome selector and locus field, e.g. `chr1:0-2,000,000`.
  Coordinates are **zero-based, half-open**, matching poreckit. 
- Use `+` / `−` to zoom, arrows to pan, and back/forward for locus history.
  Click the chromosome overview to reposition the window.
- Switch triangle/square view, resolution, normalization and colormap.
  `auto` selects an available resolution targeting about 250 bins. Raw
  counts work without normalization weights. Other normalization choices
  require corresponding weights in the input file.
- Display changes apply automatically after a short debounce; disable **Auto
  update** to apply changes with **Update map** instead. Available resolutions
  and normalization names come from the selected file.
- The default is a square, equal-aspect, linear white–red observed-contact map.
  **Reset display settings** restores these display settings. **Color range** can
  use the 95th/99th percentile of positive off-diagonal contacts, the full data
  range, or a manual upper limit. Saturation changes colors without modifying
  counts. The header reports the applied view, normalization, scale and color
  maximum. The reset changes only the display mode, colormap, log scale and color range;
  it preserves the file, locus, resolution and normalization.
- Add contact-map panels, GTF/GFF genes, bigWig/bedGraph signals, BED TADs or BEDPE loops.
  Check to show/hide and drag to reorder; Auto update applies the changes.
  Double-click a track or click **Settings** to edit its name, mode, color,
  placement and height. Contact panels can independently be square or triangular;
  signals support filled, line, bar and point views. New annotation tracks default
  to **Follow map layout**: square maps put genes/signals above and to the left,
  TADs in diagonal boxes, and loops in mirrored anchor boxes. Triangle maps stack
  annotations below the contact map. **Placement** also supports above, left,
  both axes, separate panels, or in-map overlays as appropriate to each track.
  Existing explicitly configured placements are preserved; use **Settings →
  Placement → Follow map layout** to adopt the new defaults for an older track.
  All panels share the locus. Multiple files and repeated tracks are supported.
  **Source chromosome** handles naming differences such as `1` and `chr1`;
  it does not convert assemblies. Contact maps and bigWig readers also try the
  equivalent prefix when the exact chromosome is absent.
  Axis-side tracks share the matrix coordinates and keep signal values intact
  when rotated; the two orientations reuse the same underlying data query.
- Load a second contact file with **File → Load control map…** to unlock split, difference and ratio views against it, including observed/expected modes that use each file's own expected curves.
- **View → Whole-genome overview…** (`Ctrl+G`) aggregates contacts across the whole genome; click a chromosome pair to open it in the browser.
- The **Window** menu opens additional browser windows, duplicates the current view, and can link loci across windows so navigation stays synchronized.
- **File → Export matrix data…** saves the displayed square map as NPZ with the observed and displayed arrays, the optional control matrix and coordinate metadata.
- Export PNG, PDF or SVG. Save/open JSON sessions to retain settings and
  references to data files. Relative paths resolve against the session file;
  data files are not copied. Demo data is temporary and is not session-persistent.
- Try **File → Open synthetic demo** without external datasets. These counts are synthetic and have no biological meaning. **Open synthetic track demo** fills the Tracks / alignments tab with example files.

The UI follows the desktop's Qt style, fonts and palette, without a custom
stylesheet. Qt is cross-platform; it does not use GTK-native widgets on GNOME.
All application text is in English. Qt file dialogs keep their controls in
English even when the desktop uses another language.
The interface uses locus navigation, editable track lists and a shared genomic view.

## Scope and limitations

The **Genome Browser** tab supports a primary contact map and additional aligned contact/annotation panels. The **Tracks / alignments** tab browses signal, alignment, variant, reference and annotation files at its own locus. The **Assembly editor** tab reorders, cuts, joins and inverts scaffolds with automatic coordinate lifting and a live remapped contact-map preview. The **Plot Workbench** tab exposes every plotting command, public plotting function and exported track type in the installed poreckit version. The current catalog has
**30 commands, 84 functions and 38 track constructors** (including GC and
restriction-site helpers). [Full catalog](docs/plot-catalog.md).

The workbench includes inter-chromosome and genome maps, Pore-C/multiway views,
QC/reports, TADs, loops, stripes, compartments, saddles, pileups, differential
and aggregate figures, alignment/coverage/sashimi tracks, and the other
poreckit plot families. Non-plotting data-management commands are excluded.
Gene-name search uses local GTF/GFF/BED annotations. Genome data is not downloaded automatically.
Files lacking chromosome lengths need an explicit locus.

## Tracks and assemblies

**File → Open genomic tracks…** loads bigWig/bedGraph signals, BAM/CRAM/SAM alignments, FASTA references, VCF/BCF variants and GTF/GFF/BED files into the **Tracks / alignments** tab. The tab has its own locus field, per-track display options and track-session save/open. **File → Open remote track…** reads bigWig, BAM, VCF and FASTA directly over HTTP(S) with byte-range requests; an optional custom index URL (BAI/CSI/TBI) is supported.

The **Assembly editor** tab starts from the loaded contact map or an imported `.assembly` text file. Scaffolds can be dragged, moved, cut, split, joined and inverted with undo/redo; positions lift automatically and the preview re-queries the source contacts for the edited layout, bounded so very large layouts must zoom to a smaller region first. Projects save and reopen with their source mapping. **Export** writes the edited layout as `.assembly` text, AGP, corrected FASTA, a corrected contact map (`.cool`) or a map image, and never overwrites the source contact file.

## Plot Workbench

1. Find an entry in the searchable catalog. **Commands** compute plots from
   source files. **Functions** plot existing matrices/results. **Track stack**
   composes all available track types into one figure.
2. Fill the generated parameter form. Required parameters have `*`; **Show all
   optional parameters** exposes additional analysis/display settings. Output
   destinations are managed by the viewer. Blank fields retain
   poreckit's own defaults; tooltips and the description show the original
   API help. Click **Use browser data / region** to populate compatible inputs
   from the Genome Browser. **Apply locus** copies the locus into the selected
   form and every track in the stack; click Render to query the new region.
3. For track stacks, configure a track and click **Add / update track**. Repeat
   for additional tracks. Check to show/hide, drag to reorder, double-click to
   edit, then **Render track stack**. The renderer preserves poreckit's own
   overlay and independent-panel semantics. A standalone track can be rendered
   directly without first adding it to a stack.
4. Click **Render plot**. Multiple results, including figures embedded in QC
   reports, appear in the figure dropdown. All figures use the same viewer.
5. **Scroll to zoom around the pointer**, **drag to pan**, use **Zoom + / −**,
   **Fit** to reset, or **100%** for native size. Manual zoom remains unchanged
   when the window resizes. Ctrl+= / Ctrl+- / Ctrl+0 act on the active tab.
6. **Export** copies a complete SVG, PDF or PNG figure without rerunning the
   analysis. **Viewport** exports the currently visible, zoomed/panned area as
   PNG. **Export all** saves every figure plus generated tables/HTML/NPZ files
   in a ZIP archive. Results are temporary until exported.
7. **Save preset / Open preset** retain parameters, track order/visibility,
   locus and data folder. They reference inputs without copying data, and
   opening a preset does not run it. The File-menu session/export actions
   follow the active tab.

Image zoom/pan applies to all figures, including QC and aggregate plots that
have no genomic region. Genome Browser navigation changes the data window;
workbench locus edits affect region parameters and current-matrix references.
Zooming an image alone does not requery data or change its resolution.

### Function and track data

Functions that expect matrices or computed result dictionaries retain those
contracts; the GUI does not substitute a contact file for an analysis result.
Use the matching Command to compute a figure from a contact file, or import
prepared data using the function's form. Every function signature is exposed,
including extra keyword options where supported.

- Scalar strings, numeric and boolean parameters have typed controls.
- Data fields accept inline JSON. Enable **Array** for a numeric vector or
  matrix; leave it off for lists of names, dictionaries or collections.
- **Browse** imports JSON, NPY/NPZ or numeric CSV/TSV/TXT. NPY produces an array;
  NPZ produces a dictionary of arrays, or its only array if there is just one.
  Pick a specific NPZ array with `{"$file": "/path/result.npz", "key": "matrix"}`.
- Nested NumPy arrays use `{"$array": [[1, 2], [3, 4]]}` inside JSON structures.
  For example, `[{"$array": [[1, 2], [2, 3]]}]` is a collection of one matrix.
- **Current matrix** inserts `{"$matrix": true}` and queries the Genome
  Browser input using the workbench locus and current context. A fully explicit
  reference can supply `input`, `region`, `binsize` and `norm` under `$matrix`.
- Reader fields accept contact-file paths; structured references may also use
  `{"$reader": "/path/sample.mcool"}` or `{"$contacts": "/path/reads.pairs"}`.
- Repeated/tuple CLI parameters accept JSON arrays, e.g. `--viewpoint` takes
  `["chr1", 100000, 110000]`, while repeatable `--annotation` takes
  `[["genes", "/path/genes.bed"], ["peaks", "/path/peaks.bed"]]`.
- Relative paths resolve from **Data folder**, independent of the output folder.
  Python code and pickle/object-array inputs are not loaded.

[An executable example preset](examples/differential-matrix.json) uses explicitly
synthetic matrices and can be opened directly in Plot Workbench.

Rendering uses reusable, cancellable worker processes. Genome Browser retains
up to two readers and a shared 64 MiB LRU matrix/signal cache keyed by file identity/version,
locus, resolution and normalization. File edits invalidate cached reads. Preview
generates SVG only; PNG/PDF are rendered on demand from the displayed settings
and cached data. Cancelling a task stops its worker and clears its cache. Genome Browser and current-matrix
queries reject dense matrices exceeding 2500 bins per side. Workbench commands
also check explicit region/resolution pairs; genome-wide, reader-based track
and analysis-specific allocations follow poreckit's own limits. Small bins may
still require reading large source files. Original errors appear in the run log;
failed or cancelled workbench tasks retain previous results.

## Development and distribution

```sh
uv run pytest -q
```

Headless tests use Qt's offscreen platform and actual poreckit rendering. They
check catalog completeness against source definitions, structured data loading,
CLI argument conversion, representative commands and API functions, multi-figure
QC/reports, track composition, GUI zoom/pan/export, presets and cancellation.
Catalog coverage does not mean every plot has been run against biological data.
Only Linux execution has been verified locally. The CI workflow runs the
same smoke tests on Windows, Linux and macOS; those results require CI execution.

Build native installers on each target OS with PyInstaller:

```sh
uv sync --locked --group build
uv run --no-sync python scripts/build_release.py
```

Outputs in `dist/releases/`: Windows setup EXE and portable ZIP, macOS DMG
(Intel and Apple Silicon separately), and Linux DEB/Arch `.pkg.tar.zst`. The build runs
the frozen application's GUI/worker self-test before packaging. GitHub Actions
builds all five targets (including separate Debian/Arch environments) and publishes version tags only after all jobs pass.
Optional Windows signing and macOS signing/notarization are supported.
See [packaging and CI instructions](docs/packaging.md) for prerequisites,
private core-repository access, signing secrets, verification and release steps.
