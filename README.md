# Confinement Analysis — Monte Carlo Simulation \& TPM Spring Constants

Tools for studying particle confinement on a 2D lattice. Lattice Monte Carlo simulations generate
trajectories under different interaction strengths (`J`), which are then visualised as confinement
panels, while a parallel pipeline turns experimental single-molecule tracking data into per-pixel
spring-constant maps. The shared metric throughout is `k = 1 / Rg²` (radius of gyration → confinement
strength).

## Folders

|#|Folder|Purpose|
|-|-|-|
|1|[`1-MC\_sim\_Uniform\_J`](1-MC_sim_Uniform_J/README.md)|Run lattice Monte Carlo simulations (`J ∈ {0, −2, −3}`, multiple seeds) and write trajectory + zoomed `Rg²` map `.npz` files. Includes a script to derive spring constants from experimental tracks.|
|2|[`2-plot\_confinement\_panels\_uniform\_J`](2-plot_confinement_panels_uniform_J/README_plot_confinement_panels_uniform_J.md)|Plot trajectory panels for the three confinement regimes plus a cumulative frequency distribution (CFD) of `log(k)`.|
|3|[`3-plot\_confinement\_panels\_multizoom\_uniform\_J`](3-plot_confinement_panels_multizoom_uniform_J/README_plot_confinement_panels_multizoom.md)|Multi-scale variant: trajectory fields at 160/80/40 nm pixel scales with an overlaid CFD across scales.|
|4|[`4-TPM\_analysis`](4-TPM_analysis/README_TPM_Analysis.md)|Notebook pipeline converting PALMTracer SPT data into per-pixel spring-constant maps (single- or dual-label), with TIFF maps and statistics.|

## Typical Workflow

1. **Simulate** (folder 1) → produces `.npz` trajectories and zoomed `Rg²` maps.
2. **Visualise** (folders 2 \& 3) → read those `.npz` files and render confinement figures + CFDs.
3. **Experimental analysis** (folder 4) → convert real PALMTracer tracks into spring-constant maps using the same `k = 1 / Rg²` estimator.

## Confinement Regimes (`J`)

|`J`|Regime|
|-|-|
|`0`|Free diffusion (Brownian)|
|`−2`|Weakly confined|
|`−3`|Strongly confined|

As `|J|` increases, confinement strengthens and trajectories become more spatially restricted.

## Requirements

Python 3.8+ with `numpy`, `matplotlib`, `scipy`. Folder 1 also needs `numba` (simulation) and `pandas`
(experimental script); folder 4 additionally needs `pandas`, `openpyxl`, `tqdm`, and optionally
`scikit-image`, run from Jupyter.

```bash
pip install numpy matplotlib scipy numba pandas openpyxl tqdm scikit-image
```

See each folder's README for detailed inputs, parameters, and outputs.

