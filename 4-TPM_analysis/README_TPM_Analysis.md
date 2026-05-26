# TPM Analysis
It converts single-particle tracking output from PALMTracer into per-pixel spring constant maps at one or more spatial sampling scales, in either single-label or dual-label mode.

## What It Does

- Loads `trcPALMTracer*.txt` trajectory files from one or more cell folders
- Computes a per-pixel spring constant `k = 1 / Rg²` at every zoom level in `ZOOM_VALUES`
- Saves a tab-delimited per-pixel table for every (cell, zoom) combination
- Optionally renders one TIFF map per zoom level using a discrete colour scale
- In single-label mode, writes a combined descriptive-statistics Excel file across all cells and zooms
- In dual-label mode, renders a three-panel TIFF per cell and zoom: protein 1 (red gradient), protein 2 (blue gradient), and a five-colour `log10(k1) − log10(k2)` difference overlay
- Writes a plain-text `config_snapshot.txt` recording every CONFIG value used in the run

The notebook is configured entirely from a single CONFIG cell. No other cell needs to be edited for routine use.

## Repository Files

- `TPM_Analysis.ipynb`: main pipeline notebook (CONFIG + functions + run cell + diagnostic cell)
- *input* directory (per cell, user-supplied): `trcPALMTracer*.txt` and an optional `.rgn` boundary file
- Generated outputs (under `OUTPUT_ROOT`, one sub-folder per cell):
  - `<Protein>_cell<N>_zoom<Z>x_spring_constants.txt` — per-pixel table
  - `<Protein>_cell<N>_zoom<Z>x_map.tiff` — single-label discrete-colour map (when `SAVE_IMAGES = True`)
  - `<cell>_zoom<Z>x_dual.tiff` — dual-label three-panel figure (when `SAVE_IMAGES = True`)
  - `spring_constant_stats_summary.xlsx` — single-label descriptive statistics
  - `config_snapshot.txt` — text record of every CONFIG value used

## Requirements

- Python 3.8+
- Jupyter Notebook or JupyterLab
- Packages:
  - `numpy` (≥1.21)
  - `pandas` (≥1.3)
  - `matplotlib` (≥3.5)
  - `scipy` (≥1.7) — `ndimage.zoom` for bilinear upsampling
  - `openpyxl` (≥3.0) — formatted Excel output
  - `tqdm` (≥4.0) — progress bar for parallel processing
  - `scikit-image` (optional) — fast polygon rasterisation for boundary masking (matplotlib `Path` fallback used if absent)

Install dependencies:

```bash
pip install numpy pandas matplotlib scipy openpyxl tqdm scikit-image
```

## Input Data Format

### Trajectory files

```text
<input_folder>/trcPALMTracer*.txt
```

Tab-delimited, exported directly from PALMTracer. Each file has:

- Row 0 — metadata: tab-separated key–value pairs including `Width`, `Height` (camera frame in pixels), `nb_Tracks`, etc.
- Row 1 — column names: `Track`, `Plane`, `CentroidX(px)`, `CentroidY(px)`, and additional columns.
- Rows 2+ — one localisation per row.

When multiple `trcPALMTracer*.txt` files exist in one folder (e.g. from multiple acquisition sessions), they are concatenated and `Track` IDs are made globally unique by adding a per-file offset, so that Track 1 from file A and Track 1 from file B are treated as different tracks.

### Boundary file (optional)

```text
<input_folder>/<dendrite>.rgn
```

A MetaMorph `.rgn` polygon file. Field 6 stores the polygon as a flat list of `(row, col)` integer coordinates in camera pixels (row → x axis, col → y axis; the parser swaps them to standard image convention). When provided, the boundary is overlaid on every image; with `MASK_TO_BOUNDARY = True`, out-of-boundary pixels are additionally set to white.

### Expected folder structure

For **single-label mode**, the protein name and cell number are auto-extracted from the folder path:

```text
.../ProteinName/cellN/loc and trc files/trcPALMTracer*.txt
```

For **dual-label mode**, the cell name is supplied explicitly as a dict key in `DUAL_INPUT_FOLDERS`; the two protein folders themselves can have any structure.

## Important Parameters

All parameters live in the CONFIG cell.

- `N_PROTEINS` (default `1`)
  - `1` — single-label mode; one discrete-colour map per zoom level.
  - `2` — dual-label mode; three-panel red / blue / difference figure per zoom level.
- `INPUT_FOLDERS` (single-label only)
  - Dict mapping each trc folder to its `.rgn` boundary file (or `None`).
- `DUAL_INPUT_FOLDERS` (dual-label only)
  - Dict mapping each cell name to a 3-tuple `(protein1_folder, protein2_folder, rgn_path_or_None)`. The third entry must be present; use `None` when no boundary is available.
- `PROTEIN1_NAME`, `PROTEIN2_NAME`
  - Display names for the red and blue channels in dual-label mode (e.g. `"PSD95"`, `"GluA2"`).
- `NM_PER_PIXEL` (default `160.0`)
  - Camera calibration in nm per pixel.
- `ZOOM_VALUES` (default `[1, 2, 4, 8, 16]`)
  - Pixel binning factors. `bin_size = NM_PER_PIXEL / zoom`. At `zoom = 4` and `NM_PER_PIXEL = 160`, the bin is 40 nm. Smaller bins resolve finer spatial structure but require more trajectories per pixel for reliable estimates.
- `MIN_POINTS` (default `6`)
  - Minimum localisations per track (total across all frames). Tracks shorter than this are discarded before pixel binning. Use `6` for ensemble SPT; `1`–`2` for sparse or single-trajectory data.
- `VMIN`, `VMAX`
  - `log10(k)` range for the colour scale. Leave both as `None` on a first single-label run to auto-derive from data; the printed values can then be pasted back to lock the scale across subsequent runs. Dual-label mode **requires** both to be set explicitly.
- `USE_LOC_PRECISION_FLOOR` (default `False`)
  - `False` — single-localisation segments (`Rg² = 0`) → NaN, rendered white, excluded from statistics. Use for ensemble SPT.
  - `True` — single-localisation segments → `k = 1 / (10 nm)²`, corresponding to typical SMLM localisation precision. Use for bead or immobile-control data.
- `BASE_COLORS`
  - 5-hex list for the single-label discrete colormap, low `k` → high `k`. Default: yellow → orange → purple → blue → green.
- `BASE_COLORS_RED`, `BASE_COLORS_BLUE`
  - 6-hex lists for the white → dark red and white → dark blue gradients used in dual-label mode.
- `THRESH_WEAK`, `THRESH_STRONG` (defaults `0.3`, `1.0`)
  - Difference-map thresholds on `|log10(k1) − log10(k2)|`. `|diff| ≤ THRESH_WEAK` → black (co-confined); `THRESH_WEAK < |diff| ≤ THRESH_STRONG` → light red / blue; `|diff| > THRESH_STRONG` → dark red / blue.
- `MASK_TO_BOUNDARY` (default `False`)
  - `True` — set pixels outside the `.rgn` polygon to white. `False` — show the full field of view.
- `BOUNDARY_LINEWIDTH` (default `0.3`)
  - Line width of the ROI boundary overlay.
- `MAX_WORKERS` (default `None`)
  - `None` — all CPU cores. `1` — serial (use if memory pressure occurs in parallel runs).
- `SAVE_IMAGES` (default `True`)
  - `False` — skip all TIFF rendering; still saves `.txt` tables and `.xlsx` stats. Useful for fast numerical-only runs.
- `DPI` (default `300`)
  - Output image resolution. Use `600` or `1200` for high-resolution journal submission.

## Run

Open the notebook in Jupyter:

```bash
jupyter notebook TPM_Analysis.ipynb
```

Then:

1. Run the **Imports** cell.
2. Edit the **CONFIG** cell — at minimum: set `N_PROTEINS`, fill in `INPUT_FOLDERS` or `DUAL_INPUT_FOLDERS`, set `OUTPUT_ROOT`.
3. Run all remaining function-definition cells in order.
4. Run the **Run** cell. It dispatches to `main()` (single-label) or `main_dual()` (dual-label) based on `N_PROTEINS`.

Outputs are written under `OUTPUT_ROOT` in sub-folders named by protein and cell. Typical execution time scales linearly with the number of cells and the number of zoom levels; on a modern workstation, one cell at all five default zooms takes about 1–3 minutes.

### First run: finding VMIN and VMAX

For dual-label mode, both `VMIN` and `VMAX` must be set before the run. The **Diagnostic** cell (last cell of the notebook) loads the first folder in `INPUT_FOLDERS`, computes `log10(k)` at each zoom level, and prints summary statistics:

```text
zoom=1x   n=   843  min=-3.21  P2=-2.89  median=-1.72  P98=0.44  max=1.01
zoom=2x   n=  3241  min=-2.88  P2=-2.61  median=-1.34  P98=0.89  max=1.55
zoom=4x   n= 11065  min=-2.55  P2=-2.30  median=-0.91  P98=1.23  max=2.10
zoom=8x   n= 38412  min=-2.22  P2=-1.98  median=-0.58  P98=1.61  max=2.74
```

Use these values to set `VMIN` and `VMAX` before re-running the main pipeline. Locking the scale ensures all cells and zoom levels share a consistent colormap.

## Methodology

- **Spring constant estimator.** `k = 1 / Rg²` with `Rg² = (1/n) Σ[(xⱼ − x̄)² + (yⱼ − ȳ)²]`. This is the finite-sample estimator of `⟨|r − r₀|²⟩ = 2 kBT / k*` from the overdamped Langevin equation. `kBT ≡ 1`, so the reported value is `k = k* / (2 kBT)` — a constant offset that does not affect any relative comparison.
- **Method** Each track is split by pixel bin and `Rg²` is computed independently within each segment. Multiple segments from different tracks that share the same pixel are averaged. This ensures each pixel reflects only the local confinement measured within that spatial region.
- **Vectorised implementation.** Pixel assignment, per-segment deviation, per-segment `Rg²` aggregation, and per-pixel averaging are all expressed as pandas `transform()` and `agg()` calls, keeping the inner loops in compiled code.
- **Boundary masking.** Uses `skimage.draw.polygon` for fast rasterisation, falling back to `matplotlib.path.Path.contains_points` when scikit-image is unavailable.
- **Rendering pipeline (display only).** Raw `log10(k)` is upsampled 4× with bilinear interpolation (`scipy.ndimage.zoom`, `order=1`) **before** colour mapping. `imshow` then receives the raw values directly and maps them through a `ListedColormap` + `BoundaryNorm`, so colour-bin boundaries remain discrete while the spatial appearance is smooth. NaN pixels (no data, boundary-masked, or excluded single-localisation segments) render as white via `cmap.set_bad`.
- **Dual-label difference panel.** Both proteins are processed independently. The two `log10(k)` grids are computed on the same camera pixel coordinate system (shared `frame_w`, `frame_h` from protein 1's trc header). `diff = log10(k1) − log10(k2)` is mapped through the five-colour scheme using `THRESH_WEAK` and `THRESH_STRONG`. Pixels where only one protein has data fall back to that protein's gradient colour.

## Reproducibility

- A `config_snapshot.txt` recording every CONFIG value (including every input folder path) is written to `OUTPUT_ROOT` at the start of every run, providing a complete audit trail.
- All algorithmic operations are deterministic given the same input. The only sources of non-determinism are file-system iteration order (handled by `sorted()` in `load_trc_files`) and parallel-execution order (output order is independent of run order, since each cell writes to its own sub-folder).

## Notes

- Display transforms (4× bilinear upsampling, colour-bin assignment) affect rendered TIFFs only. The `.txt` tables and the Excel statistics are computed from raw per-pixel `Rg²` values.
- Single-localisation segments (`Rg² = 0`) are handled by `USE_LOC_PRECISION_FLOOR` — see Important Parameters.
- `parse_spine_rois()` is included in the file I/O utilities for ad-hoc spine-level analysis but is not called by the main pipeline.
- `compute_spring_constants()` returns a **single DataFrame**. Unpacking as `df, _ = compute_spring_constants(...)` raises `TypeError`.

## Common Issues

- **No `trcPALMTracer*.txt` files found**
  - Filenames must start with `trcPALMTracer` and end with `.txt` exactly (the match is case-sensitive on Linux/macOS).
- **All pixels are white / no data**
  - Lower `MIN_POINTS` to `2` or `1` if your tracks are short. Verify the trc files contain the columns `Track`, `CentroidX(px)`, `CentroidY(px)`.
- **`ValueError: not enough values to unpack (expected 3, got 2)` in dual-label mode**
  - One of the `DUAL_INPUT_FOLDERS` values is a 2-tuple. Every entry must be a 3-tuple `(protein1_folder, protein2_folder, rgn_path_or_None)` — use `None` for the third element when no boundary file is available.
- **`ValueError: VMIN and VMAX must be set` in dual-label mode**
  - Auto-derive is not supported for dual-label runs. Run single-label mode (or the Diagnostic cell) first to find the data range, then set both values explicitly.
- **Memory error during parallel processing**
  - Set `MAX_WORKERS = 1` to force serial execution.
- **Boundary polygon looks rotated by 90°**
  - The parser swaps row → x and col → y from field 6. If the overlay appears rotated, the `.rgn` file may use a non-standard field layout; inspect with a text editor and adjust `parse_rgn()` if needed.
- **`TypeError: cannot unpack non-sequence`**
  - `compute_spring_constants()` returns one DataFrame, not a tuple. Use `df = compute_spring_constants(...)`.
- **Colormap looks different from previous runs**
  - `VMIN` / `VMAX` were auto-derived and differed between runs. Paste the printed values into CONFIG to lock the scale.
- **TIFF generation is slow**
  - Set `SAVE_IMAGES = False` to skip all rendering; the `.txt` tables and stats Excel are still produced. Or reduce `ZOOM_VALUES`. Or increase `MAX_WORKERS` if memory allows.
