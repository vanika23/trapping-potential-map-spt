# Plot Condinement panels multizoom

**Script:** `plot_confinement_panels_multizoom.py`

This variant builds a **four-panel** figure: three trajectory fields at different pixel sampling scales (1×, 2×, 4×), plus a cumulative frequency distribution of `log₁₀(k)` across those scales. Trajectory segments use the same dominant-cell / in-cell `Rg²` logic as the TPM zoom pipeline (`get_zoomed_RgMat`), so the displayed paths align with the CFD data source.


## What It Does

1. Loads Monte Carlo trajectories for three confinement regimes (`J = 0`, `−2`, `−3`) from `MC_simulation_traj/Uniform_J/`.
2. For each pixel scale (160, 80, 40 nm):
   - Draws a **640 × 640 nm** display window (4×4 cells at 160 nm) with the matching pixel grid; finer panels show the same field of view with more grid lines (8×8 at 80 nm, 16×16 at 40 nm). The underlying simulation lattice remains `L = 960` nm.
   - Plots decimated, in-cell trajectory segments for particles that pass the zoom-scale validity check (≥ `MIN_CELL_POINTS` steps in the dominant grid cell, finite in-cell `Rg²`).
3. Builds a KDE-smoothed **CFD** (cumulative frequency distribution, 0–100%) of `log₁₀(k) = log₁₀(1 / Rg²)` from pre-computed `Rg2Mat` maps at all three scales, pooled over `J ∈ {0, −2, −3}`.
4. Writes `pixel_scale_figure_FINAL_exp_1.png` (300 dpi) and `pixel_scale_figure_FINAL_exp_1.pdf`.

## Figure Layout

| Panel | Content | Pixel size |
|-------|---------|------------|
| 1 | Trajectories + grid | 160 nm (1×) |
| 2 | Trajectories + grid | 80 nm (2×) |
| 3 | Trajectories + grid | 40 nm (4×) |
| 4 | CFD of `log₁₀(k)` | 160, 80, 40 nm (overlaid) |

**Colours:** Brownian `#E8479A`, weakly confined `#2E5FA3`, strongly confined `#2A7D3F`.

One extra strongly confined trace (picked with seed `77`, gated at 40 nm) is drawn on top in green in **all three** trajectory panels.

The figure uses `GridSpec(1, 5)` with width ratios `[1, 1, 1, 0.18, 1.2]` (narrow spacer before the CFD panel). On-screen preview is `figsize=(16.0, 4.2)` at 150 dpi; saved PNG uses 300 dpi.

## Requirements

- Python 3.8+
- `numpy` (≥1.21)
- `matplotlib` (≥3.5)
- `scipy` (≥1.7) — `gaussian_filter1d`, `gaussian_kde`

```bash
pip install numpy matplotlib scipy
```

## Input Data

Default data directory (set in script):

```text
./MC_simulation_traj/Uniform_J
```

Run the script from the repository root so this relative path resolves correctly.

### Trajectory files

```text
MC_simulation_traj/Uniform_J/J={J}_seed_{seed}.npz
```

| Key | Description |
|-----|-------------|
| `xhist` | x-coordinates, shape `(n_steps, n_particles)` |
| `yhist` | y-coordinates, same shape |

| File (default `SEED = 10`) | Regime |
|------------------------------|--------|
| `J=0_seed_10.npz` | Brownian |
| `J=-2_seed_10.npz` | Weakly confined |
| `J=-3_seed_10.npz` | Strongly confined |

### Zoom / CFD files

```text
MC_simulation_traj/Uniform_J/J={J}_seed_{seed}_zoom_{1,2,4}x.npz
```

| Key | Description |
|-----|-------------|
| `Rg2Mat` | Per-pixel `Rg²` (nm²); flattened for the CFD |

| Zoom suffix | Pixel size |
|-------------|------------|
| `zoom_1x` | 160 nm |
| `zoom_2x` | 80 nm |
| `zoom_4x` | 40 nm |

The CFD loader pools `J = 0`, `−2`, and `−3` for each scale. Missing files are skipped without error.

### Example layout

```text
vanika/
  plot_confinement_panels_multizoom.py
  MC_simulation_traj/
    Uniform_J/
      J=0_seed_10.npz
      J=-2_seed_10.npz
      J=-3_seed_10.npz
      J=0_seed_10_zoom_1x.npz
      J=0_seed_10_zoom_2x.npz
      J=0_seed_10_zoom_4x.npz
      ... (same zoom files for J=-2, J=-3)
```

## Run

```bash
python plot_confinement_panels_multizoom.py
```

**Outputs:**

- `pixel_scale_figure_FINAL_exp_1.png`
- `pixel_scale_figure_FINAL_exp_1.pdf`

A preview window opens via `plt.show()`. Comment out that line for headless runs.

## Key Parameters

| Symbol / name | Default | Role |
|---------------|---------|------|
| `DATA_DIR` | `./MC_simulation_traj/Uniform_J` | Input folder |
| `SEED` | `10` | Simulation seed for trajectories and zoom `.npz` files |
| `L` | `960` | Lattice size (nm), periodic boundaries |
| `RES` | `5` | Sub-cell resolution (nm); matches `MC_Model2.py` |
| `MIN_CELL_POINTS` | `6` | Minimum in-cell steps for valid `Rg²` (matches `get_zoomed_RgMat`) |
| `GRID_CELLS_COARSEST` | `4` | Coarsest panel shows a 4×4 pixel block |
| `FIELD_NM` | `640` | Plot extent (nm): `4 × 160` |
| `DECIMATE` | `5` | Use every Nth step when drawing trajectories |
| `SMOOTH_SIG` | `2.0` | Gaussian smoothing σ on displayed traces (display only) |
| `TARGET_B`, `TARGET_W`, `TARGET_S` | `52`, `20`, `8` nm | Per-class display `Rg` targets (display only) |

**Pixel sizes** are defined in `ZOOM_PANELS`: `zoom_1x` → 160 nm, `zoom_2x` → 80 nm, `zoom_4x` → 40 nm.

**Particle counts** (after validity gating): 14 Brownian, 10 weak, 14 strong, plus one optional highlighted strong trace.

## Trajectory Panel Logic

Unlike `pixel_scale_figure_FINAL_refined.py`, this script does **not** randomly reassign trajectory centres across the field. For each particle and scale:

1. **Dominant cell** — grid cell with the most decimated steps.
2. **In-cell mask** — only steps inside that cell are drawn; others are `NaN`.
3. **Validity** — `pick_particles` keeps only particles with finite in-cell `Rg²` at the chosen zoom (same rules as the zoom map pipeline).
4. **Shared particle set** — particles are picked **once**, gated at the finest scale (40 nm). A particle valid at 40 nm is necessarily valid at 80 nm and 160 nm (its dense sub-cell sits inside the coarser cells), so the **same** trajectories appear in all three panels. A localisation cannot appear at one scale and vanish at another.
5. **Display scaling** — each confinement class is scaled about its segment centroid so Brownian, weak, and strong traces are visually separable (`TARGET_* / median(Rg)`). Does not affect the CFD.
6. **Periodic jumps** — large step discontinuities within the segment are masked; light `gaussian_filter1d` smoothing is applied for display.

## CFD Panel

- Values: `log₁₀(1 / Rg²)` from positive, finite entries in `Rg2Mat`.
- Smoothing: Gaussian KDE (`bw_method=0.15`), integrated and normalised to **0–100%** (cumulative frequency).
- Axis labels: `log k` (x), `Cumulative Frequency (%)` (y); title `CFD`.
- Axis limits: x ∈ [−2.8, 2.5], y ∈ [0, 105].

## Reproducibility

Fixed random seeds in `pick_particles`:

| Seed | Use |
|------|-----|
| `10` | Brownian particle IDs (14 particles) |
| `20` | Weakly confined particle IDs (10 particles) |
| `30` | Strongly confined particle IDs (14 particles) |
| `77` | Extra highlighted strong trace (all three trajectory panels) |

All picks are filtered with `zoom_nm=40` (finest scale) before assignment to any panel.


## Common Issues

- **`FileNotFoundError` on `J=..._seed_10.npz`**
  - Confirm `MC_simulation_traj/Uniform_J/` exists and you run from the repo root (or update `DATA_DIR`).

- **Sparse or empty trajectory panels**
  - Few particles may pass in-cell validity at 40 nm; this is expected when `MIN_CELL_POINTS` is not met in the dominant cell.

- **Empty CFD curve**
  - Missing `*_zoom_{1,2,4}x.npz` for the chosen `SEED` and `J` values; check filenames include the negative sign in `J=-2`.

- **Import errors**
  - `pip install numpy matplotlib scipy`
