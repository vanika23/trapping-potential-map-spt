"""
plot_confinement_panels_uniform_J.py
==================
Analysis and visualisation script for the 2D lattice Monte Carlo (MC) simulation
used to validate the Trapping Potential Mapping (TPM) pipeline.

The simulation generates trajectories of particles on a periodic 960 x 960 nm² lattice
under three confinement regimes (Brownian, Weakly Confined, Strongly Confined),
controlled by the harmonic coupling parameter J. This script:

    1. Loads pre-generated trajectory files from the Uniform_J/ folder
    2. Computes log(k) = log10(1/Rg²) per particle directly from raw trajectories
    3. Plots trajectory panels for each confinement class (shared field of view)
    4. Plots cumulative density functions (CFDs) of log(k) for each class
    5. Saves all outputs as high-resolution PNG files

Usage:
    python plot_confinement_panels_uniform_J.py

Requirements:
    numpy, matplotlib, scipy
    Uniform_J/ folder containing J=0, J=-2, J=-3 trajectory .npz files

Output files:
    sim_trajectories.png   — trajectory panels + CFD (600 dpi)
    CFD_confinement.png    — CFD only (600 dpi)
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')   # non-interactive backend for script use
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d

# ── GLOBAL STYLE ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'DejaVu Sans',
    'font.size':        11,
    'axes.labelsize':   12,
    'xtick.labelsize':  10,
    'ytick.labelsize':  10,
    'pdf.fonttype':     42,
    'ps.fonttype':      42,
    'svg.fonttype':     'none',
})

# ── COLOURS AND LABELS ────────────────────────────────────────────────────────
# Three confinement classes:
#   J = 0  → Brownian (free diffusion, no restoring force)
#   J = -2 → Weakly confined (intermediate harmonic trap)
#   J = -3 → Strongly confined (deep harmonic trap)
COLORS = {
    0:  '#E8479A',   # Brownian       — hot pink
    -2: '#2E5FA3',   # Weakly conf.   — blue
    -3: '#2A7D3F',   # Strongly conf. — green
}
LABELS = {
    0:  'Brownian',
    -2: 'Weakly Confined',
    -3: 'Strongly Confined',
}

# Data directory — folder containing J=*_seed_*.npz files
DATA_DIR = './MC_simulation_traj/Uniform_J_fast'

# Shared field of view for all three trajectory panels (lattice units).
# All panels use the same axis range so trajectory sizes are visually comparable.
FIELD_R = 120


# ── DATA LOADING ──────────────────────────────────────────────────────────────

def load_trajectory(J, seed=1):
    """
    Load raw trajectory arrays for a given J value and random seed.

    Parameters
    ----------
    J    : int   — coupling parameter (0, -2, -3, ...)
    seed : int   — random seed index (1-10)

    Returns
    -------
    xhist, yhist : ndarray (n_steps+1, n_particles), float
                   Returns (None, None) if file not found.
    """
    fname = os.path.join(DATA_DIR, f'J={J}_seed_{seed}.npz')
    if not os.path.exists(fname):
        print(f'  Warning: {fname} not found')
        return None, None
    d = np.load(fname)
    return d['xhist'].astype(float), d['yhist'].astype(float)


def compute_rg_per_particle(xhist, yhist, L=960):
    """
    Compute the radius of gyration Rg for each particle trajectory.

    Rg² = (1/n) Σⱼ [(xⱼ - x₀)² + (yⱼ - y₀)²]

    where x₀ = mean(x) and y₀ = mean(y) are the trajectory centroids,
    corresponding to the equilibrium position r₀ = <r> from the theory.

    Particles whose trajectory crosses the periodic boundary (identified by
    a single step larger than L/4) are excluded and assigned Rg = inf.

    Parameters
    ----------
    xhist, yhist : ndarray (n_steps, n_particles)
    L            : int — lattice size (nm), default 960

    Returns
    -------
    rg : ndarray (n_particles,) — Rg in lattice units; inf for excluded particles
    """
    rg = []
    for i in range(xhist.shape[1]):
        x = xhist[:, i]
        y = yhist[:, i]
        # Exclude particles that cross periodic boundaries
        if np.any(np.abs(np.diff(x)) > L // 4) or \
           np.any(np.abs(np.diff(y)) > L // 4):
            rg.append(np.inf)
        else:
            rg2 = np.mean((x - x.mean())**2 + (y - y.mean())**2)
            rg.append(np.sqrt(rg2))
    return np.array(rg)


def compute_logk(J, seeds):
    """
    Compute log(k) = log10(1/Rg²) for all particles across multiple seeds.

    This uses per-particle Rg computed directly from raw trajectories —
    not the pre-binned zoom files — to ensure the CFD reflects the true
    per-trajectory spring constant distribution independently of pixel size.

    k = 1/Rg² is the effective spring constant (k_BT set to 1).
    The true spring constant k* = 2·k_BT/Rg² under the isotropic harmonic model;
    k = k*/2 is used consistently throughout for all comparisons.

    Parameters
    ----------
    J     : int        — coupling parameter
    seeds : list[int]  — list of seed indices to pool

    Returns
    -------
    logk : ndarray — log(k) values for all valid particles across seeds
    """
    logk = []
    for seed in seeds:
        xh, yh = load_trajectory(J, seed)
        if xh is None:
            continue
        for i in range(xh.shape[1]):
            x = xh[:, i]; y = yh[:, i]
            rg2 = np.mean((x - x.mean())**2 + (y - y.mean())**2)
            if rg2 > 0:
                logk.append(np.log10(1.0 / rg2))
    return np.array(logk)


# ── TRAJECTORY PANEL HELPERS ──────────────────────────────────────────────────

def place_non_overlapping(n, min_sep, seed=13):
    """
    Place n trajectory centres randomly within FIELD_R with minimum separation.

    Uses a simple rejection sampler. Returns up to n positions;
    may return fewer if packing is too dense within the attempt budget.

    Parameters
    ----------
    n       : int   — number of positions
    min_sep : float — minimum centre-to-centre separation (lattice units)
    seed    : int   — random seed for reproducibility

    Returns
    -------
    positions : list of (x, y) tuples
    """
    rng = np.random.default_rng(seed)
    positions = []
    bound = FIELD_R * 0.70
    attempts = 0
    while len(positions) < n and attempts < 5000:
        attempts += 1
        px = rng.uniform(-bound, bound)
        py = rng.uniform(-bound, bound)
        if all(np.sqrt((px - ex)**2 + (py - ey)**2) > min_sep
               for ex, ey in positions):
            positions.append((px, py))
    return positions


def select_particles(rg, J, n_show):
    """
    Select particle indices for display based on confinement class.

    - Brownian   (J=0):  particles with Rg > 80% of median (large, fills panel)
    - Weakly     (J=-2): 1 high-spread particle + (n_show-1) near-median
    - Strongly   (J=-3): particles with Rg within 0.4×–2.5× of median

    Parameters
    ----------
    rg     : ndarray — per-particle Rg values (inf = excluded)
    J      : int     — coupling parameter
    n_show : int     — number of trajectories to display

    Returns
    -------
    chosen : ndarray of int — particle indices
    """
    fin = np.where(np.isfinite(rg))[0]
    med = np.median(rg[fin])
    rng = np.random.default_rng(42 + abs(J))

    if J == 0:
        # Brownian: pick large-Rg particles so they fill the shared field
        good   = fin[rg[fin] > med * 0.8]
        chosen = rng.choice(good, size=min(n_show, len(good)), replace=False)

    elif J == -2:
        # Weakly confined: mix one high-spread and several typical
        high   = fin[rg[fin] > med * 1.4]
        mid    = fin[np.abs(rg[fin] - med) < med * 0.6]
        pick_h = rng.choice(high, size=min(1, len(high)), replace=False)
        pick_m = rng.choice(mid,  size=min(n_show - 1, len(mid)), replace=False)
        chosen = np.concatenate([pick_h, pick_m])

    elif J == -3:
        # Strongly confined: typical particles near median Rg (exclude extremes)
        good   = fin[(rg[fin] > med * 0.4) & (rg[fin] < med * 2.5)]
        chosen = rng.choice(good, size=min(n_show, len(good)), replace=False)

    else:
        chosen = rng.choice(fin, size=min(n_show, len(fin)), replace=False)

    return chosen


def plot_trajectory_panel(ax, J, seed=1, n_show=3, decimate=4, smooth_sigma=0):
    """
    Plot trajectory panel for one confinement class.

    Trajectories are centred on their mean position and placed at random
    non-overlapping positions within the shared field (±FIELD_R lattice units).

    For confined classes (J=-2, J=-3), a Gaussian filter is applied to the
    raw lattice trajectories to produce smooth, curved steps that visually
    resemble continuous diffusion — matching the appearance of experimental
    sptPALM data. The underlying Rg statistics are unaffected (Rg is computed
    before smoothing).

    Parameters
    ----------
    ax           : matplotlib Axes
    J            : int   — confinement class (0, -2, -3)
    seed         : int   — trajectory file seed
    n_show       : int   — number of trajectories to show
    decimate     : int   — plot every nth step (reduces file size)
    smooth_sigma : float — Gaussian filter sigma (0 = no smoothing)
    """
    col = COLORS[J]

    xhist, yhist = load_trajectory(J, seed)
    if xhist is None:
        ax.text(0.5, 0.5, 'data not found', ha='center',
                va='center', transform=ax.transAxes, color='gray')
        ax.axis('off')
        return

    # Compute Rg for particle selection
    rg = compute_rg_per_particle(xhist, yhist)
    chosen = select_particles(rg, J, n_show)

    # Placement
    traj_rg   = np.median(rg[chosen])
    min_sep   = traj_rg * 2.8   # manually tuned heuristic visual separation factor
    pos_seed  = 99 if J == 0 else 13 + abs(J) # seed selection for reproducibility
    positions = place_non_overlapping(n_show, min_sep, seed=pos_seed)

    for i, (xoff, yoff) in zip(chosen, positions):
        x = xhist[::decimate, i]
        y = yhist[::decimate, i]

        # Apply Gaussian smoothing BEFORE centring.
        # This converts axial-only MC steps into curved paths resembling
        # continuous Brownian motion, without changing the overall trajectory shape.
        if smooth_sigma > 0:
            x = gaussian_filter1d(x.astype(float), sigma=smooth_sigma)
            y = gaussian_filter1d(y.astype(float), sigma=smooth_sigma)

        # Centre trajectory on its mean position, offset to placement location
        xcm = np.nanmean(x); ycm = np.nanmean(y)
        xc  = x - xcm + xoff
        yc  = y - ycm + yoff

        # Clip any remaining large jumps (periodic boundary artefacts)
        th  = rg[i] * 6
        bad = (np.abs(np.diff(xc, prepend=xc[0])) > th) | \
              (np.abs(np.diff(yc, prepend=yc[0])) > th)
        xc[bad] = np.nan
        yc[bad] = np.nan

        ax.plot(xc, yc, '-', color=col, lw=0.8, alpha=0.85, rasterized=True)

    # All panels share the same axis limits — size differences are meaningful
    ax.set_xlim(-FIELD_R, FIELD_R)
    ax.set_ylim(-FIELD_R, FIELD_R)
    ax.set_aspect('equal')

    # Minimal border
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#BBBBBB')
        spine.set_linewidth(0.8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(LABELS[J], fontsize=12, fontweight='bold', color=col, pad=5)


def plot_cfd(ax, J, seeds, bins):
    """
    Compute and plot the CFD of log(k) for one confinement class.

    CFD is computed via KDE on per-particle log(k) values pooled across seeds.
    KDE bandwidth is set to 0.25 to smooth over the discrete lattice structure
    without over-smoothing the separation between classes.

    Parameters
    ----------
    ax    : matplotlib Axes
    J     : int        — coupling parameter
    seeds : list[int]  — seeds to pool
    bins  : ndarray    — x-axis evaluation points
    """
    logk = compute_logk(J, seeds)
    if len(logk) == 0:
        print(f'  Warning: no log(k) values for J={J}')
        return

    kde  = gaussian_kde(logk, bw_method=0.25)
    dens = kde(bins)
    dx   = bins[1] - bins[0]
    cfd  = np.cumsum(dens * dx)
    cfd  = (cfd / cfd[-1]) * 100  # normalise to 0–100%

    ax.plot(bins, cfd, color=COLORS[J], lw=2.2, solid_capstyle='round')


# ── MAIN FIGURE ───────────────────────────────────────────────────────────────

def make_figure():
    """
    Generate the combined trajectory + CFD figure.

    Layout: 3 trajectory panels (Brownian | Weakly Confined | Strongly Confined)
    + 1 CFD panel, all on a single row.

    Trajectory panels share a fixed field of view (±FIELD_R lattice units)
    so spatial spread is directly comparable across confinement classes.
    """

    # Per-class display settings
    # n_show       : number of trajectories to display
    # decimate     : plot every nth step (lower = denser lines)
    # smooth_sigma : Gaussian filter sigma applied to confined trajectories
    # cfd_seeds    : which seeds to pool for the CFD
    class_config = {
        0:  dict(n_show=2, decimate=4, smooth_sigma=0, cfd_seeds=[1]),
        -2: dict(n_show=4, decimate=2, smooth_sigma=2, cfd_seeds=[1, 2, 3]),
        -3: dict(n_show=6, decimate=1, smooth_sigma=2, cfd_seeds=[1]),
    }
    J_values = [0, -2, -3]

    fig = plt.figure(figsize=(15, 4.0))
    gs  = GridSpec(1, 4, figure=fig,
                   width_ratios=[1, 1, 1, 1.3],
                   wspace=0.18,
                   left=0.02, right=0.98,
                   top=0.88, bottom=0.12)

    traj_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    cfd_ax    = fig.add_subplot(gs[0, 3])

    # ── Trajectory panels ────────────────────────────────────────────────
    for ax, J in zip(traj_axes, J_values):
        cfg = class_config[J]
        print(f'Plotting trajectory panel: J={J} ({LABELS[J]})')
        plot_trajectory_panel(
            ax, J,
            seed         = 1,
            n_show       = cfg['n_show'],
            decimate     = cfg['decimate'],
            smooth_sigma = cfg['smooth_sigma'],
        )

    # ── CFD panel ────────────────────────────────────────────────────────
    bins = np.linspace(-5, 3.5, 500)
    for J in J_values:
        cfg = class_config[J]
        print(f'Plotting CFD: J={J} ({LABELS[J]}), seeds={cfg["cfd_seeds"]}')
        plot_cfd(cfd_ax, J, cfg['cfd_seeds'], bins)

    cfd_ax.set_xlabel(r'$\log(k)$', fontsize=12)
    cfd_ax.set_ylabel('Cumulative Frequency Distribution (%)', fontsize=10)
    cfd_ax.set_xlim(-5, 3.5)
    cfd_ax.set_ylim(-0.02, 101)
    cfd_ax.axhline(0, color='#CCCCCC', lw=0.6, zorder=0)
    cfd_ax.spines['top'].set_visible(False)
    cfd_ax.spines['right'].set_visible(False)
    cfd_ax.tick_params(direction='out', length=3)
    # Legend omitted — add manually in figure editing software

    return fig


# ── STANDALONE CFD FIGURE ─────────────────────────────────────────────────────

def make_cfd_figure():
    """
    Generate a standalone CFD figure (without trajectory panels).
    Useful for supplementary figures or presentations.
    """
    J_values   = [0, -2, -3]
    cfd_seeds  = {0: [1], -2: [1, 2, 3], -3: [1]}

    fig, ax = plt.subplots(figsize=(5, 5))
    bins = np.linspace(-5, 3.5, 500)

    for J in J_values:
        plot_cfd(ax, J, cfd_seeds[J], bins)

    ax.set_xlabel(r'$\log(k)$', fontsize=12)
    ax.set_ylabel('Cumulative Frequency Distribution (%)', fontsize=12)
    ax.set_xlim(-5, 3.5)
    ax.set_ylim(-0.02, 101)
    ax.axhline(0, color='#CCCCCC', lw=0.6, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out', length=3)

    return fig


# ── ENTRY POINT ───────────────────────────────────────────────────────────────

if __name__ == '__main__':

    print('=' * 60)
    print('MC Simulation Analysis — Trapping Potential Mapping (TPM)')
    print('=' * 60)
    print(f'Data directory: {DATA_DIR}/')
    print()

    # Check data is present
    for J in [0, -2, -3]:
        ok = os.path.exists(os.path.join(DATA_DIR, f'J={J}_seed_1.npz'))
        print(f'  J={J:3d} ({LABELS[J]:20s}): {"OK" if ok else "FILE NOT FOUND"}')
    print()

    # Combined trajectory + CFD figure
    print('Generating combined figure (trajectories + CFD)...')
    fig = make_figure()
    out = 'sim_trajectories.png'
    fig.savefig(out, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  Saved: {out}')
    print()

    # Standalone CFD
    print('Generating standalone CFD figure...')
    fig_cfd = make_cfd_figure()
    out_cfd = 'CFD_confinement.png'
    fig_cfd.savefig(out_cfd, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.close(fig_cfd)
    print(f'  Saved: {out_cfd}')
    print()
    print('Done.')
