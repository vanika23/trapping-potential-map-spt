# Trapping Potential Mapping (TPM)

Spatially resolved maps of confinement strength from single-particle tracking data.

TPM treats the transient confinement of a membrane molecule as overdamped Brownian
motion in an effective harmonic potential, and converts the squared radius of
gyration of each trajectory segment into a local trapping stiffness
`k = 2·kBT / Rg²`. Assigning `log10(k)` to a defined pixel grid renders a
continuous, per-pixel map of confinement strength across a cell.

This repository contains the complete analysis code for:

> Dhingra V., Choudhary P., Sarkar S., Nair D.
> *Nanoscale Trapping Potential Mapping in Neurons by Single-Particle Tracking.*

A detailed description of the algorithm, including the derivation of the
estimator and all analysis parameters, is given in the **Methods** section of the
manuscript (sections *"Near a minimum, all potentials can be approximated as a
harmonic potential"*, *"Brownian motion of a particle in a harmonic trap"*,
*"Simulated trajectory validation"*, and *"Statistical analysis"*).

---

## Contents

| # | Folder | Purpose |
|---|--------|---------|
| 1 | [`1-MC_sim_Uniform_J/`](1-MC_sim_Uniform_J/) | Lattice Monte Carlo simulation of confined diffusion (`J ∈ {0, −2, −3}`, 10 seeds). Writes trajectory and zoomed `Rg²` map `.npz` files. |
| 2 | [`2-plot_confinement_panels_uniform_J/`](2-plot_confinement_panels_uniform_J/) | Trajectory panels for the three confinement regimes plus a cumulative frequency distribution (CFD) of `log10(k)`. Generates Supplementary Fig. 1A–D. |
| 3 | [`3-plot_confinement_panels_multizoom/`](3-plot_confinement_panels_multizoom/) | Multi-scale variant: trajectory fields at 160/80/40 nm pixel scales with an overlaid CFD across scales. Generates Supplementary Fig. 4. |
| 4 | [`4-TPM_analysis/`](4-TPM_analysis/) | Main pipeline. Converts PALMTracer SPT output into per-pixel `log10(k)` maps (single- or dual-label), TIFF renderings, per-pixel tables, and summary statistics. Generates Figs. 1C–D, 2, 3. |
| — | [`demo_data/`](demo_data/) | Real sptPALM demo dataset (PALMTracer exports from one hippocampal neuron). See §3.1. |

Each folder has its own README with full parameter documentation. **Start with
this file**, then consult the folder README for the component you need.

### Software not included here

Two upstream steps use established third-party tools and are not reimplemented
in this repository:

- **PALMTracer** — single-molecule detection, sub-pixel localisation and
  trajectory assembly (wavelet segmentation, 2-D Gaussian fitting, simulated
  annealing linking). MetaMorph plug-in, available from
  <https://neuro-intramuros.u-bordeaux.fr/displayresearchprojects/70/11>.
- **ThunderSTORM** — generation of the simulated immobilised-bead TIFF stacks.
  ImageJ plug-in, freely available at <https://zitmen.github.io/thunderstorm/>.

The TPM pipeline consumes PALMTracer `trcPALMTracer*.txt` exports as its input.

---

## 1. System requirements

### Operating systems

| OS | Status |
|----|--------|
| Windows 10 / 11 (64-bit) | Primary development and analysis platform; all manuscript figures produced here |
| Ubuntu 24.04 LTS (64-bit) | Verified: full demo runs end-to-end |
| macOS 12+ | Expected to work (pure Python, no OS-specific calls); not formally tested |

No operating-system-specific code is used. Note that trajectory filename
matching (`trcPALMTracer*.txt`) is **case-sensitive on Linux and macOS**.

### Software dependencies

Python **3.9 – 3.12**.

| Package | Minimum | Versions verified | Used by |
|---------|---------|-------------------|---------|
| `numpy` | 1.21 | 1.23.5, 2.4.4 | all |
| `scipy` | 1.7 | 1.9.3, 1.17.1 | folders 2, 3, 4 |
| `matplotlib` | 3.5 | 3.5.3, 3.10.8 | all |
| `pandas` | 1.3 | 1.4.4, 3.0.2 | folders 1, 4 |
| `numba` | 0.56 | 0.66.0 | folder 1 only |
| `openpyxl` | 3.0 | 3.1.5 | folder 4 |
| `tqdm` | 4.0 | 4.66 | folder 4 |
| `scikit-image` | 0.19 | 0.26.0 | folder 4 (optional; matplotlib `Path` fallback used if absent) |
| `jupyter` / `notebook` | — | 7.x | folder 4 |


Pinned versions for exact reproduction are in [`requirements.txt`](requirements.txt).

### Hardware

No non-standard hardware is required and no GPU is used. Any desktop or laptop
capable of running Python 3.9+ is sufficient.

**Memory is the binding constraint, and it depends on the finest zoom level you
render.** TIFF rendering upsamples the `log10(k)` grid 4× before colour mapping,
so peak memory scales roughly as `zoom²`. Measured on a 201 × 215 pixel camera
frame with 472,000 localisations:

| `ZOOM_VALUES` | Finest bin | Peak RAM | Wall time |
|---------------|-----------|----------|-----------|
| `[1, 2]` | 80 nm | 0.6 GB | 6 s |
| `[1, 2, 4]` | 40 nm | 1.4 GB | 9 s |
| `[1, 2, 4, 8]` | 20 nm | **~4 GB** | 30 s |
| `[1, 2, 4, 8]` with `SAVE_IMAGES = False` | 20 nm | 0.3 GB | 5 s |

Recommended:

- **8 GB RAM** to render the full four-scale series including the 20 nm maps.
  With 4 GB or less, either drop `8` from `ZOOM_VALUES` or set
  `SAVE_IMAGES = False` — the numerical results (`.txt` tables, statistics
  workbook) are identical either way, and cost 0.3 GB. Rendering, not
  computation, is what exhausts memory.
- 4+ CPU cores (folder 4 parallelises across cells with `ProcessPoolExecutor`;
  set `MAX_WORKERS = 1` to run serially or if memory is tight).
- **~1 GB free disk** if you re-run the full Monte Carlo simulation in folder 1
  (30 trajectory `.npz` files of ~32 MB each, plus 120 zoom-map files).

---

## 2. Installation guide

```bash
git clone https://github.com/vanika23/trapping-potential-map-spt.git
cd trapping-potential-map-spt
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

**Typical install time on a normal desktop computer: 2–4 minutes** on a standard
broadband connection (dominated by downloading `numba`, `scipy` and
`matplotlib` wheels). Cloning the repository itself takes a few seconds — it
contains no large binary files.

To verify the installation:

```bash
python -c "import numpy, scipy, pandas, matplotlib, numba, tqdm, openpyxl; print('ok')"
```

---

## 3. Demo

The demo runs the main TPM pipeline (folder 4) on a real sptPALM dataset and
reproduces, in miniature, the scale-dependence result of Fig. 1D and the
emergence of discrete high-`k` foci shown in Fig. 1C.

### 3.1 The demo dataset

`demo_data/` contains two PALMTracer demo datasets from live hippocampal
neurons: single-label mode (`demo_data/ProteinA`) and dual-label mode
(`demo_data/dual_label_mode`).

**Single-label demo**
| | |
|---|---|
| Files | `trcPALMTracer1.txt` … `trcPALMTracer5.txt` (5 acquisition sessions, concatenated by the pipeline) |
| Camera frame | 201 × 215 pixels at 160 nm/px (32.2 × 34.4 µm) |
| Frames per acquisition | 4,000 |
| Localisations | 472,062 |
| Tracks | 104,293 (about 19% pass the `MIN_POINTS = 6` filter) |
| Size on disk | 36 MB total |

Track IDs restart at 1 in each file; `load_trc_files()` offsets them so tracks
from different sessions are never merged.

If you prefer a lighter download, `trcPALMTracer1.txt` alone (8.6 MB, 112,947
localisations) exercises every code path and produces the same qualitative
result. Expected values for both configurations are given below.

**Dual-label demo**
Two channels: `cell01/ProteinX` and `cell01/ProteinY`, 
with a shared boundary ROI, `cell01/ProteinX_dendrite.rgn`, sitting alongside
both protein folders rather than inside either one, that's the path `DUAL_INPUT_FOLDERS` third tuple element points to.

||ProteinX|ProteinY|
|---|---|---|
|Files| `trcPALMTracer1.txt` … `trcPALMTracer5.txt`| `trcPALMTracer1.txt` … `trcPALMTracer5.txt`|
|Camera frame| 200 × 200 px at 160 nm/px (32.0 × 32.0 µm)| same |
|Frames per acquisition | 5,000 | same|
| Localisations | 120,528| 179,480 |
| Tracks | 5,090 (100% pass the `MIN_POINTS = 6` filter) | 4,395 (100% pass the `MIN_POINTS = 6` filter)|
| Size on disk | 9.0 MB | 14.0 MB|

Unlike the single-label demo, every track in both channels clears
`MIN_POINTS = 6` on the first pass, so no pixels are dropped for being too
short — useful for checking the dual-label difference-map thresholds without
track-length noise confounding the picture.

Track IDs restart at 1 in each file, in both datasets and both channels;
`load_trc_files()` offsets them so tracks from different sessions are never
merged.

### 3.2 Run the pipeline

```bash
jupyter notebook 4-TPM_analysis/TPM_Analysis.ipynb
```

In the **CONFIG** cell set:

```python
N_PROTEINS    = 1
INPUT_FOLDERS = {r"../demo_data/DemoProtein/cell1/loc_and_trc_files": None}
OUTPUT_ROOT   = r"../demo_output"
NM_PER_PIXEL  = 160.0
ZOOM_VALUES   = [1, 2, 4]      # add 8 for 20 nm maps; needs ~8 GB RAM (see §1)
MIN_POINTS    = 6
VMIN = None
VMAX = None
MASK_TO_BOUNDARY = False
```

Then run all cells top to bottom. The `ProteinName/cellN/` levels of the path
are parsed to name the outputs, so keep that folder layout.

### 3.3 Expected output

Console (all five files, `ZOOM_VALUES = [1, 2, 4, 8]`):

```
  -> Config snapshot saved: demo_output/config_snapshot.txt

Auto-deriving colormap scale...
  -> VMIN=-3.926, VMAX=4.455

=== demo_data/DemoProtein/cell1/loc_and_trc_files ===
  DemoProtein_cell1: 472,062 localisations, 104,293 unique tracks  [frame 201x215 px]
    zoom=1x (160nm/px):  11044 pixels  log10k median=-2.84
    zoom=2x (80nm/px):   30133 pixels  log10k median=-1.70
    zoom=4x (40nm/px):   69258 pixels  log10k median=-0.55
    zoom=8x (20nm/px):  122945 pixels  log10k median=+0.67
  -> Stats Excel saved: demo_output/spring_constant_stats_summary.xlsx

All done.
```

With `trcPALMTracer1.txt` only and `ZOOM_VALUES = [1, 2, 4]`:

```
  DemoProtein_cell1: 112,947 localisations, 23,474 unique tracks  [frame 201x215 px]
    zoom=1x (160nm/px):  6379 pixels  log10k median=-2.91
    zoom=2x (80nm/px):  14091 pixels  log10k median=-1.81
    zoom=4x (40nm/px):  25624 pixels  log10k median=-0.65
```

Files written under `demo_output/`:

```
demo_output/
├── config_snapshot.txt                          (complete record of every CONFIG value)
├── spring_constant_stats_summary.xlsx           (descriptive statistics, all zooms)
└── DemoProtein_cell1/
    ├── DemoProtein_cell1_zoom{1,2,4,8}x_spring_constants.txt
    └── DemoProtein_cell1_zoom{1,2,4,8}x_map.tiff
```

Each `*_spring_constants.txt` is a tab-delimited per-pixel table with columns
`xid`, `yid`, `x (nm)`, `y (nm)`, `Rg²`, `k`, `log10(k)`. Blank `Rg²`/`k` entries
are pixels whose only contributing segment was a single localisation
(`Rg² = 0`); these render white and are excluded from statistics when
`USE_LOC_PRECISION_FLOOR = False`.

The `*_map.tiff` files are discrete five-colour `log10(k)` maps (yellow = low `k`
→ green = high `k`), NaN pixels white. At 160 nm the map is spatially diffuse;
by 40 nm the dendritic arbor is clearly resolved and discrete high-`k` foci
appear at synaptic sites — the behaviour described in Fig. 1C.

Two things a reviewer should expect and not mistake for errors: most pixels are
blank (only about 19% of tracks reach the 6-localisation minimum), and the
median `log10(k)` rises as the pixel size falls, which is the scale-dependence
effect reported in Fig. 1D.

### 3.4 Expected run time

Measured on a single 2.10 GHz Intel Xeon core with 3 GB RAM under Ubuntu 24.04 —
a deliberately modest reference machine, so a normal desktop will be faster:

| Configuration | Wall time |
|---|---|
| All 5 files, `[1, 2, 4]`, images on | ~9 s |
| All 5 files, `[1, 2, 4, 8]`, images on | ~30 s (needs ~4 GB RAM) |
| All 5 files, `[1, 2, 4, 8]`, `SAVE_IMAGES = False` | ~5 s |
| `trcPALMTracer1.txt` only, `[1, 2, 4]`, images on | ~5 s |

## 4. Instructions for use

### 4.1 Running TPM on your own SPT data

1. **Export trajectories from PALMTracer.** Bead registration and localisation
   processing in PALMTracer write `trcPALMTracer*.txt` (one file per stack
   processed) into a `.PT` folder next to the raw image stack.

   **Don't point `INPUT_FOLDERS` at the `.PT` folder itself.** It keeps
   accumulating other PALMTracer output — drift-correction data, registration
   logs, state files — as you reprocess, so it's a moving target for a path
   baked into CONFIG. Instead, copy just the `trcPALMTracer*.txt` files out
   into their own folder, organised as:

   ```
   <root>/<ProteinName>/cell<N>/trcPALMTracer*.txt
   ```

   The protein name and cell number are parsed from the two directory levels
   above the trc files, so this layout is required in single-label mode.
   Multiple `trcPALMTracer*.txt` files in one folder are concatenated, with
   `Track` IDs offset so they remain globally unique — this is also why all
   trc files for one cell should sit together in the same folder.

2. **Optionally add a boundary ROI** — a MetaMorph `.rgn` polygon file in the
   same folder — to restrict the map to a dendrite or spine.

3. **Edit the CONFIG cell** of `4-TPM_analysis/TPM_Analysis.ipynb`. At minimum:
   `N_PROTEINS`, `INPUT_FOLDERS` (or `DUAL_INPUT_FOLDERS`), `OUTPUT_ROOT`,
   `NM_PER_PIXEL` (your camera calibration).

4. **Choose the colour scale.** On a first single-label run leave
   `VMIN = VMAX = None` to auto-derive, then paste the printed values back into
   CONFIG so the scale is locked across all cells and conditions. Dual-label
   mode requires both to be set explicitly.

5. **Run all cells.** `config_snapshot.txt` is written at the start of every run
   and records every parameter used, providing a complete audit trail.

Key analysis parameters and the values used in the manuscript:

| Parameter | Manuscript value | Meaning |
|-----------|------------------|---------|
| `NM_PER_PIXEL` | 160.0 | EMCCD camera calibration |
| `ZOOM_VALUES` | `[1, 2, 4, 8]` | Bin size = `NM_PER_PIXEL / zoom` → 160, 80, 40, 20 nm |
| `MIN_POINTS` | 6 | Minimum localisations per track; below this the Rg estimate is substantially biased by sampling error |
| `USE_LOC_PRECISION_FLOOR` | `False` for SPT, `True` for bead/immobile controls | Handling of single-localisation segments |
| `THRESH_WEAK`, `THRESH_STRONG` | 0.3, 1.0 | Dual-label difference-map thresholds on \|Δlog10(k)\| |

See [`4-TPM_analysis/README_TPM_Analysis.md`](4-TPM_analysis/README_TPM_Analysis.md)
for the complete parameter reference, methodology notes, and troubleshooting.

### 4.2 Running the simulations

```bash
cd 1-MC_sim_Uniform_J
python MC_Model2.py
```

Writes 30 trajectory files and 120 zoom-map files to `Uniform_J_fast/`
(~1 GB total). Run time: **approximately 3–5 minutes** on the reference machine
above (Numba JIT compilation adds ~4 s on the first call).

Then, from the repository root:

```bash
python 2-plot_confinement_panels_uniform_J/plot_confinement_panels_uniform_J.py
python 3-plot_confinement_panels_multizoom/plot_confinement_panels_multizoom.py
```

Both scripts read from `./MC_simulation_traj/`, so either run them from a
directory with that layout or edit `DATA_DIR` at the top of each script.


---

## 5. Reproduction instructions

To reproduce the quantitative results in the manuscript:

| Figure | Component | Notes |
|--------|-----------|-------|
| Fig. 1C, 1D | Folder 4, single-label, PSD95 datasets | `ZOOM_VALUES = [1, 2, 4, 8]`, `VMIN`/`VMAX` locked to the manuscript scale (−3.7, +3.7) |
| Fig. 2A–D | Folder 4, single-label, MB-PSD95 / Homer / GluA1 | Statistics exported to `spring_constant_stats_summary.xlsx`, then Kruskal–Wallis with Dunn's post-hoc in GraphPad Prism v8.4.2 |
| Fig. 3A–C | Folder 4, dual-label (`N_PROTEINS = 2`), PSD95 + GluA2 | `THRESH_WEAK = 0.3`, `THRESH_STRONG = 1.0` |
| Supp. Fig. 1A–D | Folders 1 → 2 | `J ∈ {0, −2, −3}`, seeds 1–10 |
| Supp. Fig. 3 | Folder 4 with `USE_LOC_PRECISION_FLOOR = True` | TetraSpeck beads and ThunderSTORM-simulated beads through an identical pipeline |
| Supp. Fig. 4 | Folders 1 → 3 | `SEED = 10` |

Statistical tests reported in the manuscript (Shapiro–Wilk, Kruskal–Wallis with
Dunn's post-hoc, Mann–Whitney U with Bonferroni correction) were performed in
**GraphPad Prism v8.4.2** on the `log10(k)` values exported by this pipeline,
not in Python.

### Determinism

- **Folder 4 (TPM pipeline)** is fully deterministic. File iteration order is
  fixed by `sorted()`, and each cell writes to its own output sub-folder, so
  parallel execution order does not affect results.
- **`make_demo_data.py`** uses a fixed seed and is byte-for-byte reproducible.
- **Folder 1 (Monte Carlo)** is *not* currently bit-reproducible. `MC_Model2.py`
  calls `np.random.seed(seed)` outside the Numba-jitted `run_mc()`; Numba
  maintains its own PRNG state, which NumPy's seed does not reach, so the
  Metropolis move sequence differs between runs. Ensemble statistics
  (`log10(k)` distributions pooled over 10 seeds) are unaffected, but individual
  trajectory files will not be identical.

---

## 6. Data availability

`demo_data/` contains a real sptPALM dataset from one live hippocampal neuron
(5 PALMTracer acquisition files, 472,062 localisations), included so that the
pipeline can be run and verified without access to the full experimental
dataset. See §3.1.


Source data for all main and supplementary figures are available from the
corresponding authors on reasonable request.

## 7. License

Released under the [MIT License](LICENSE), an
[OSI-approved](https://opensource.org/licenses/MIT) open source license.
It permits unrestricted use, modification and redistribution, including for
commercial purposes, provided the copyright notice and permission notice are
retained. No warranty is given.

---

## 8. Citation

```bibtex
@article{dhingra_tpm,
  title   = {Nanoscale Trapping Potential Mapping in Neurons by Single-Particle Tracking},
  author  = {Dhingra, Vanika and Choudhary, Prakhar and Sarkar, Sumantra and Nair, Deepak},
  journal = {},
  year    = {},
  doi     = {}
}
```

## 9. Contact

Vanika Dhingra, Centre for Neuroscience, Indian Institute of Science, Bangalore.

Corresponding authors: Deepak Nair (deepak@iisc.ac.in)
