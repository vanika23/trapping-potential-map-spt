# MC simulation trajectories

This folder contains a **lattice Monte Carlo** simulation for particle trajectories on a 2D periodic grid, plus a separate script that derives **pixel-level spring constants** from experimental single-molecule tracking data.

## Contents

| File | Role |
|------|------|
| [`MC_helper.py`](MC_helper.py) | Core routines: neighbor checks, Metropolis MC moves, radius of gyration (with periodic boundaries), and spatial maps of mean \(R_g^2\) on a zoomed grid. |
| [`MC_Model2.py`](MC_Model2.py) | Driver script: initializes particles, runs simulations for several interaction strengths and seeds, saves trajectories and zoomed \(R_g^2\) maps as `.npz` files. |
| [`generate_spring_constants.py`](generate_spring_constants.py) | Post-processes **PALM/tracer** tab-separated exports (not the MC output): aggregates tracks per pixel, computes \(R_g\), and writes per-pixel spring constant estimates to text files. |

## Dependencies

- **Python 3** with `numpy`, `matplotlib`, **`numba`** (required for `MC_helper` / `MC_Model2`).
- **`pandas`** is only needed for `generate_spring_constants.py`.

Install example:

```bash
pip install numpy matplotlib numba pandas
```

Run scripts from this directory so that `from MC_helper import *` resolves correctly:

```bash
cd MC_simulation_traj
python MC_Model2.py
```

## `MC_helper.py` — simulation primitives

- **`isempty(i, j, lat, L)`** — With periodic boundaries on an `L×L` lattice, returns which von Neumann neighbors are empty.
- **`run_mc(r, lat, J, nEpoch)`** — Picks a random particle, proposes a hop to a random empty neighbor. Energy change is modeled as `dE = -J` for an accepted move (Metropolis: always accept if `dE ≤ 0`, else accept with probability `exp(-dE)`). Returns full time histories `xhist`, `yhist` with shape `(nEpoch + 1, nPart)`.
- **`get_Rg(x, y, L)`** — Mean-squared radius of gyration for a single track using periodic images.
- **`get_zoomed_RgMat(res, zoom, xhist, yhist, L, nPart)`** — Coarse-grains the lattice into cells of size related to `zoom` and `res` (nm). For each cell, averages \(R_g^2\) over particles that have **at least 6** time points inside that cell (same threshold idea as in the experimental pipeline).

## `MC_Model2.py` — what it runs and what it writes

**Setup (defaults in code):**

- Output directory: **`Uniform_J_fast/`** (created if missing).
- **Seeds:** `1` … `10`.
- **Particles:** `nPart = 200` on a square lattice **`L = 960`** (no two particles on the same site initially; random placement).
- **Couplings:** `J ∈ {0, -2, -3}`.
- **MC length:** `nEpoch = 10000` steps per run.

**Outputs:**

1. **Full trajectories** — `Uniform_J_fast/J=<J>_seed_<seed>.npz`  
   Fields include at least: `xhist`, `yhist`, `L`, `nEpoch`, `J`, `seed`, `nPart`.

2. **Zoomed \(R_g^2\) maps** — `Uniform_J_fast/J=<J>_seed_<seed>_zoom_<n>x.npz`  
   For `zoom ∈ {160, 80, 40, 20}` (comment in code: 160 nm reference), with `res = 5` nm used in the grid logic. Contains `xc`, `yc`, `Rg2Mat` from `get_zoomed_RgMat`.

To change the experiment, edit the `prefix`, `nPart`, `L`, `J` list, `nEpoch`, or zoom list at the top of `MC_Model2.py`.

## `generate_spring_constants.py` — experimental spring constants

This script is **not** coupled to `MC_Model2.py` by default. It:

- Loops over **cells** `[1, 2, 3, 4, 6, 7, 8]` and protein label **`PF11`**.
- Expects tracer files under  
  `../Data-27-09-22/PF11/cell<cell>/loc and trc files/trcPALMTracer<id>.txt`  
  (merges `trcPALMTracer0` … `trcPALMTracer11` when present).
- Writes  
  `../Data-27-09-22/PF11/spring_constant/cell<cell>/timeAvgd_spring_constant_zoom_<zoom>x.txt`  
  for `zoom ∈ {1, 4, 8}`.

Each output line encodes pixel indices, physical coordinates (nm), \(R_g\), a spring constant proxy \(K \propto 1/R_g\), and `log10(K)`. Header in file:  
`<xid> <yid> <x in nm> <y in nm> <Rg> <Spring Const K> <log10(K)>`.

**Paths are relative** to where you run the script; adjust `datadir` / `saveDir` if your data live elsewhere. The script uses `pandas.DataFrame.append`; on **pandas 2.0+** that API was removed—you may need an older pandas or a small refactor to `pd.concat`.

## Typical workflow

1. **Simulate:** `python MC_Model2.py` → fills `Uniform_J_fast/` with `.npz`.
2. **Analyze:** load `.npz` in Python or in your project’s analysis notebooks (e.g. parent-folder scripts that plot pixel-scale figures).
3. **Experimental spring maps:** only if you have the matching `Data-27-09-22` layout—run `generate_spring_constants.py` after fixing paths and pandas compatibility.

## License / attribution

Add your lab’s preferred citation or license here if this code accompanies a publication.
