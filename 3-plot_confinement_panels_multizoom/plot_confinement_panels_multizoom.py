# =============================================================================
# plot_confinement_panels_multizoom.py
#
# Generates the pixel sampling scale demonstration figure for:
#   "Trapping Potential Mapping (TPM)"
#
# Figure consists of four panels:
#   1–3 — 960×960 nm fields at 1×, 2×, 4× pixel sampling (160, 80, 40 nm),
#         showing cell-local trajectory segments derived from the same zoom
#         pipeline as the CDF (get_zoomed_RgMat / *_zoom_{1,2,4}x.npz).
#   4   — Cumulative Frequency Distribution (CFD) of log10(k) at all three scales.
#
# DATA SOURCE:
#   Trajectories: Uniform_J_fast/J={J}_seed_{seed}.npz  (xhist, yhist)
#   Zoom maps:    Uniform_J_fast/J={J}_seed_{seed}_zoom_{1,2,4}x.npz (Rg2Mat)
#
# DEPENDENCIES: numpy, matplotlib, scipy
# =============================================================================

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde

DATA_DIR = './MC_simulation_traj/Uniform_J'
SEED = 10
L = 960
RES = 5                    # nm sub-cell resolution (matches MC_Model2.py)
MIN_CELL_POINTS = 6        # minimum in-cell steps (matches get_zoomed_RgMat)

COL_BROWN  = '#E8479A'
COL_WEAK   = '#2E5FA3'
COL_STRONG = '#2A7D3F'

DECIMATE   = 5
SMOOTH_SIG = 2.0

TARGET_B = 52.0
TARGET_W = 20.0
TARGET_S =  8.0

# zoom_key, pixel size (nm), panel title
ZOOM_PANELS = [
    ('zoom_1x', 160.0, '1× (160 nm)'),
    ('zoom_2x',  80.0, '2× (80 nm)'),
    ('zoom_4x',  40.0, '4× (40 nm)'),
]

# Display FOV: the coarsest panel shows a 4×4 pixel grid, so the window is
# 4 × 160 = 640 nm. Finer panels share this FOV and scale up automatically
# (80 nm → 8×8, 40 nm → 16×16). The true simulation field is L = 960 nm.
GRID_CELLS_COARSEST = 4
FIELD_NM = GRID_CELLS_COARSEST * max(p for _, p, _ in ZOOM_PANELS)


# =============================================================================
# Data loading
# =============================================================================

def load_traj(J, seed=SEED):
    fname = os.path.join(DATA_DIR, f'J={J}_seed_{seed}.npz')
    d = np.load(fname)
    return d['xhist'].astype(float), d['yhist'].astype(float)


def pixel_dx(zoom_nm):
    return int(zoom_nm // RES)


def compute_rg_per_particle(xhist, yhist):
    n_particles = xhist.shape[1]
    return np.array([
        np.sqrt(np.nanmean(
            (xhist[:, i] - np.nanmean(xhist[:, i]))**2 +
            (yhist[:, i] - np.nanmean(yhist[:, i]))**2
        ))
        for i in range(n_particles)
    ])


def load_logk_for_scale(zoom_key):
    logk_all = []
    for J in [0, -2, -3]:
        fname = os.path.join(DATA_DIR, f'J={J}_seed_{SEED}_{zoom_key}.npz')
        try:
            rg2 = np.load(fname)['Rg2Mat'].flatten()
            rg2 = rg2[(rg2 > 0) & np.isfinite(rg2)]
            logk_all.append(np.log10(1.0 / rg2))
        except FileNotFoundError:
            pass
    return np.concatenate(logk_all) if logk_all else np.array([])


# =============================================================================
# Zoom-scale trajectory helpers (match MC_helper.get_zoomed_RgMat logic)
# =============================================================================

def dominant_cell(x, y, dx, gx, gy):
    """
    Return grid indices and nm bounds for the cell containing the most steps.
    x, y are 1D decimated coordinate arrays (may contain NaN).
    """
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() == 0:
        return None

    xi = (x[valid] // dx).astype(int)
    yi = (y[valid] // dx).astype(int)
    ok = (xi >= 0) & (xi < gx) & (yi >= 0) & (yi < gy)
    if not ok.any():
        return None

    xi, yi = xi[ok], yi[ok]
    # mode cell
    pairs = xi * gy + yi
    uniq, counts = np.unique(pairs, return_counts=True)
    best = uniq[counts.argmax()]
    ci = best // gy
    cj = best % gy

    x_lo = ci * dx
    x_hi = (ci + 1) * dx
    y_lo = cj * dx
    y_hi = (cj + 1) * dx
    return ci, cj, x_lo, x_hi, y_lo, y_hi


def in_cell_rg2(xcol, ycol, x_lo, x_hi, y_lo, y_hi):
    """Rg² using only in-cell points (≥ MIN_CELL_POINTS), periodic-aware."""
    mask = (xcol >= x_lo) & (xcol < x_hi) & (ycol >= y_lo) & (ycol < y_hi)
    n = mask.sum()
    if n < MIN_CELL_POINTS:
        return np.nan
    xc = xcol[mask].mean()
    yc = ycol[mask].mean()
    xic = xcol[mask] - xc
    xic = xic - np.round(xic / L) * L
    yic = ycol[mask] - yc
    yic = yic - np.round(yic / L) * L
    return np.mean(xic**2 + yic**2)


def particle_valid_at_zoom(xhist, yhist, particle_idx, zoom_nm):
    """True if particle has a valid in-cell Rg² at its dominant cell (zoom pipeline)."""
    dx = pixel_dx(zoom_nm)
    gx = gy = int(np.ceil(L / dx))
    x = xhist[::DECIMATE, particle_idx]
    y = yhist[::DECIMATE, particle_idx]
    cell = dominant_cell(x, y, dx, gx, gy)
    if cell is None:
        return False
    _, _, x_lo, x_hi, y_lo, y_hi = cell
    return np.isfinite(in_cell_rg2(x, y, x_lo, x_hi, y_lo, y_hi))


def prepare_cell_trajectory(xhist, yhist, particle_idx, zoom_nm, scale):
    """
    Extract in-cell trajectory segment at true spatial coordinates for a given
    pixel sampling scale. Applies per-class visual scaling about the segment centroid.
    """
    dx = pixel_dx(zoom_nm)
    gx = gy = int(np.ceil(L / dx))

    x = xhist[::DECIMATE, particle_idx].astype(float).copy()
    y = yhist[::DECIMATE, particle_idx].astype(float).copy()

    cell = dominant_cell(x, y, dx, gx, gy)
    if cell is None:
        return np.array([np.nan]), np.array([np.nan])
    _, _, x_lo, x_hi, y_lo, y_hi = cell

    outside = ~((x >= x_lo) & (x < x_hi) & (y >= y_lo) & (y < y_hi))
    x[outside] = np.nan
    y[outside] = np.nan

    if np.sum(np.isfinite(x)) < MIN_CELL_POINTS:
        return x, y

    # Visual scaling about segment centroid (display only)
    cx = np.nanmean(x)
    cy = np.nanmean(y)
    x = (x - cx) * scale + cx
    y = (y - cy) * scale + cy

    # Mask periodic jumps within segment
    rg_est = np.sqrt(np.nanmean((x - cx)**2 + (y - cy)**2))
    thresh = max(rg_est * 5.0, 5.0)
    for arr in (x, y):
        jumps = np.abs(np.diff(arr, prepend=arr[0])) > thresh
        arr[jumps] = np.nan

    valid = np.isfinite(x)
    if valid.sum() > 4:
        x[valid] = gaussian_filter1d(x[valid], sigma=SMOOTH_SIG)
        y[valid] = gaussian_filter1d(y[valid], sigma=SMOOTH_SIG)

    return x, y


def pick_particles(rg_vals, n, seed, rg_max=150, zoom_nm=None,
                   xhist=None, yhist=None):
    med = np.median(rg_vals[(rg_vals > 0) & (rg_vals < rg_max)])
    good = np.where((rg_vals > med * 0.4) & (rg_vals < med * 2.5))[0]
    if zoom_nm is not None and xhist is not None:
        good = np.array([
            i for i in good
            if particle_valid_at_zoom(xhist, yhist, i, zoom_nm)
        ])
    np.random.seed(seed)
    if len(good) == 0:
        return np.array([], dtype=int)
    return np.random.choice(good, size=min(n, len(good)), replace=False)


# =============================================================================
# Plotting helpers
# =============================================================================

def draw_pixel_grid(ax, pixel_nm, xlo=0.0, xhi=FIELD_NM, ylo=0.0, yhi=FIELD_NM):
    lw = 0.45 if pixel_nm >= 80 else 0.25
    col = '#CCCCCC' if pixel_nm >= 80 else '#DDDDDD'
    for xv in np.arange(xlo, xhi + 0.1, pixel_nm):
        ax.axvline(xv, color=col, lw=lw, zorder=0)
    for yv in np.arange(ylo, yhi + 0.1, pixel_nm):
        ax.axhline(yv, color=col, lw=lw, zorder=0)


def draw_trajectories(ax, xhist, yhist, particle_indices, zoom_nm, scale,
                      color, lw, alpha, xlo=0.0, xhi=FIELD_NM,
                      ylo=0.0, yhi=FIELD_NM):
    for pi in particle_indices:
        tx, ty = prepare_cell_trajectory(xhist, yhist, pi, zoom_nm, scale)
        outside = (tx < xlo) | (tx > xhi) | (ty < ylo) | (ty > yhi)
        tx[outside] = np.nan
        ty[outside] = np.nan
        ax.plot(tx, ty, '-', color=color, lw=lw, alpha=alpha,
                rasterized=True, zorder=2)


def format_traj_axis(ax):
    ax.set_xlim(0.0, FIELD_NM)
    ax.set_ylim(0.0, FIELD_NM)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color('#888888')
        sp.set_linewidth(0.8)


def plot_smooth_cdf(ax, data, color, label, lw=2.0):
    d = data[np.isfinite(data)]
    if d.size == 0:
        return
    kde = gaussian_kde(d, bw_method=0.15)
    x_fine = np.linspace(d.min() - 0.3, d.max() + 0.3, 2000)
    pdf = kde(x_fine)
    cdf = np.cumsum(pdf) * (x_fine[1] - x_fine[0])
    cdf = 100.0 * cdf / cdf[-1]          # scale to percent for CFD
    ax.plot(x_fine, cdf, '-', color=color, lw=lw, label=label)


def draw_zoom_panel(ax, zoom_nm, title,
                    traj_sets, extra_s=None):
    """
    Draw one trajectory panel at a given pixel sampling scale.

    traj_sets: list of (xhist, yhist, indices, scale, color, lw, alpha)
    extra_s: optional (xhist, yhist, idx, scale) for highlighted trajectory
    """
    draw_pixel_grid(ax, zoom_nm)
    for xhist, yhist, indices, scale, color, lw, alpha in traj_sets:
        draw_trajectories(ax, xhist, yhist, indices, zoom_nm, scale,
                          color, lw, alpha)
    if extra_s is not None:
        xh, yh, pi, sc = extra_s
        tx, ty = prepare_cell_trajectory(xh, yh, pi, zoom_nm, sc)
        outside = (tx < 0) | (tx > FIELD_NM) | (ty < 0) | (ty > FIELD_NM)
        tx[outside] = np.nan
        ty[outside] = np.nan
        ax.plot(tx, ty, '-', color=COL_STRONG, lw=0.70, alpha=0.92,
                rasterized=True, zorder=3)
    format_traj_axis(ax)
    ax.set_title(title, fontsize=10, pad=4)


# =============================================================================
# Load trajectories and compute display scales
# =============================================================================

xb, yb = load_traj(J=0)
xw, yw = load_traj(J=-2)
xs, ys = load_traj(J=-3)

rg_b = compute_rg_per_particle(xb, yb)
rg_w = compute_rg_per_particle(xw, yw)
rg_s = compute_rg_per_particle(xs, ys)

med_b = np.median(rg_b[(rg_b > 0) & (rg_b < 200)])
med_w = np.median(rg_w[(rg_w > 0) & (rg_w < 100)])
med_s = np.median(rg_s[(rg_s > 0) & (rg_s < 50)])

scale_b = TARGET_B / med_b
scale_w = TARGET_W / med_w
scale_s = TARGET_S / med_s

# Pick particles ONCE, gated at the finest scale (40 nm = most restrictive).
# A particle valid at 40 nm is necessarily valid at 80/160 nm (its dense
# sub-cell sits inside the coarser cells), so the identical set renders in all
# three panels — a localisation cannot appear at one scale and vanish at another.
FINEST_NM = min(p for _, p, _ in ZOOM_PANELS)
shared_picks = {
    'b': pick_particles(rg_b, 14, seed=10, zoom_nm=FINEST_NM, xhist=xb, yhist=yb),
    'w': pick_particles(rg_w, 10, seed=20, zoom_nm=FINEST_NM, xhist=xw, yhist=yw),
    's': pick_particles(rg_s, 14, seed=30, zoom_nm=FINEST_NM, xhist=xs, yhist=ys),
}
ZOOM_PICKS = {pixel_nm: shared_picks for _, pixel_nm, _ in ZOOM_PANELS}

extra_s_candidates = pick_particles(rg_s, 5, seed=77, zoom_nm=FINEST_NM, xhist=xs, yhist=ys)
extra_s_idx = extra_s_candidates[0] if len(extra_s_candidates) else None

lk_160 = load_logk_for_scale('zoom_1x')
lk_80  = load_logk_for_scale('zoom_2x')
lk_40  = load_logk_for_scale('zoom_4x')


# =============================================================================
# Build figure: 3 zoom trajectory panels + CFD
# =============================================================================

fig = plt.figure(figsize=(16.0, 4.2), dpi=150)
gs = GridSpec(1, 5, figure=fig,
              width_ratios=[1, 1, 1, 0.18, 1.2],
              wspace=0.18)

for col, (_, pixel_nm, title) in enumerate(ZOOM_PANELS):
    ax = fig.add_subplot(gs[col])
    picks = ZOOM_PICKS[pixel_nm]
    traj_sets = [
        (xh, yh, picks[key], sc, col, lw, al)
        for xh, yh, key, sc, col, lw, al in [
            (xb, yb, 'b', scale_b, COL_BROWN,  0.50, 0.72),
            (xw, yw, 'w', scale_w, COL_WEAK,   0.55, 0.75),
            (xs, ys, 's', scale_s, COL_STRONG, 0.60, 0.80),
        ]
    ]
    extra = None
    if extra_s_idx is not None:
        extra = (xs, ys, extra_s_idx, scale_s)
    draw_zoom_panel(ax, pixel_nm, title, traj_sets, extra_s=extra)

# CFD panel
ax_cdf = fig.add_subplot(gs[4])
plot_smooth_cdf(ax_cdf, lk_160, '#D94F7A', '160 nm')
plot_smooth_cdf(ax_cdf, lk_80,  '#5B9BD5',  '80 nm')
plot_smooth_cdf(ax_cdf, lk_40,  '#4CAF50',  '40 nm')

ax_cdf.set_xlabel('log k', fontsize=11)
ax_cdf.set_ylabel('Cumulative Frequency (%)', fontsize=11)
ax_cdf.set_title('CFD', fontsize=10, pad=4)
ax_cdf.set_xlim(-2.8, 2.5)
ax_cdf.set_ylim(0, 105)
ax_cdf.legend(fontsize=9, frameon=False, loc='lower right')
ax_cdf.spines['top'].set_visible(False)
ax_cdf.spines['right'].set_visible(False)
ax_cdf.tick_params(labelsize=9)

fig.patch.set_facecolor('white')

out_png = 'pixel_scale_figure_FINAL_exp_1.png'
out_pdf = 'pixel_scale_figure_FINAL_exp_1.pdf'

plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.show()

print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
