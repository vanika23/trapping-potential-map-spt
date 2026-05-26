# Plot Confinement Panels Uniform J

This repository contains a Python analysis script for 2D lattice Monte Carlo trajectory data under multiple confinement regimes.  
It loads simulation files, computes per-particle confinement metrics, and generates high-resolution visual outputs.

## What It Does

- Loads trajectory data from `.npz` files in `MC_simulation_traj/Uniform_J_fast/`
- Computes per-particle `Rg` and `log(k) = log10(1 / Rg^2)`
- Builds a combined figure with:
  - trajectory panels for three confinement classes
  - cumulative frequency distribution (CFD) of `log(k)` on a 0–100% scale
- Builds a standalone CFD figure

## Repository Files

- `plot_confinement_panels_uniform_J.py`: main analysis and plotting script
- `MC_simulation_traj/Uniform_J_fast/`: expected input directory with trajectory `.npz` files
- Generated outputs:
  - `sim_trajectories.png`
  - `CFD_confinement.png`

## Requirements

- Python 3.9+ (recommended)
- Packages:
  - `numpy`
  - `matplotlib`
  - `scipy`

Install dependencies:

```bash
pip install numpy matplotlib scipy
```

## Input Data Format

The script expects files named:

```text
MC_simulation_traj/Uniform_J_fast/J={J}_seed_{seed}.npz
```

Typical `J` values used by default are `0`, `-2`, and `-3`.
In this dataset and plotting setup:

- `J = 0`: free diffusion (Brownian-like, no restoring confinement force)
- `J = -2`: weak/intermediate confinement (weakly confined)
- `J = -3`: strong confinement (strongly confined; deep trap regime)

As `|J|` increases, confinement strength increases and trajectories are more spatially restricted.

Each `.npz` file should include:

- `xhist`: x positions, shape `(n_steps+1, n_particles)`
- `yhist`: y positions, shape `(n_steps+1, n_particles)`
- Optional metadata keys such as `J`, `seed`, `nPart`, `nEpoch`, `L`

## Important Parameters

- `FIELD_R`
  - Shared trajectory-panel half-range (in lattice units) for all confinement classes.
  - A larger value shows more empty space; a smaller value zooms in.
  - Keeping this fixed across classes makes spatial spread visually comparable.
- `smooth_sigma`
  - Gaussian smoothing strength used when drawing trajectories (mainly confined classes).
  - Affects display aesthetics only; it does not change metric computation.
  - `0` means no smoothing (Brownian panel).
- `log(k)` computation
  - Computed per particle from raw trajectories using:
    - `Rg^2 = mean((x - mean(x))^2 + (y - mean(y))^2)`
    - `log(k) = log10(1 / Rg^2)`
  - CFD curves are generated from pooled per-particle `log(k)` values across selected seeds, integrated via KDE and normalised to 0–100%.

## Run

From the repository root:

```bash
python plot_confinement_panels_uniform_J.py
```

The script checks for required input files, then writes:

- `sim_trajectories.png` (combined trajectory + CFD figure)
- `CFD_confinement.png` (standalone CFD figure)

## Using Functions in Other Code

You can import helper functions into notebooks or other scripts:

```python
from MC_Sim_analysis_refined import compute_logk, plot_trajectory_panel, plot_cfd
```

## Notes

- **Strongly confined** trajectories and CFD use **`J = -3`** simulation files directly (`J=-3_seed_{seed}.npz`).
- Gaussian smoothing (`smooth_sigma = 2`) is applied when **drawing** weakly and strongly confined trajectories to soften lattice-step artefacts; `Rg`, `log(k)`, and CFD values are always computed from **unsmoothed** raw trajectories.
- Missing files are handled with warnings and skipped where possible.

## Common Issues

- **`FILE NOT FOUND` for a class**
  - Ensure `MC_simulation_traj/Uniform_J_fast/` exists and files match `J={J}_seed_{seed}.npz`.
- **Import errors**
  - Reinstall dependencies with `pip install numpy matplotlib scipy`.
- **No output images**
  - Check console warnings for missing input files or empty datasets.
- **CFD panel looks flat or empty in `sim_trajectories.png`**
  - The combined figure y-axis must span 0–100% to match CFD normalisation in `plot_cfd`. If you change scaling in `plot_cfd`, update `cfd_ax.set_ylim` in `make_figure()` accordingly.
