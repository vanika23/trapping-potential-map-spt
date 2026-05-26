import numpy as np 
import matplotlib.pyplot as plt 
from numba import njit

# from scipy.signal import periodogram, find_peaks, find_peaks_cwt
import os

@njit
def isempty(i,j,lat,L):
    nbrs = [
        [(i+1)%L,j],
        [i,(j+1)%L],
        [(i-1)%L,j],
        [i,(j-1)%L]
            ]
    empty = False 
    empty_nbrs = []
    for nbr in nbrs:
       x,y = nbr 
       if lat[x,y] == 0:
           empty = True 
           empty_nbrs.append(nbr)
    return empty_nbrs,empty 
       
@njit
def run_mc(r, lat, J, nEpoch):
    nPart = r.shape[0]
    L = lat.shape[0]
    nT = nEpoch + 1
    xhist = np.empty((nT, nPart), dtype=np.int64)
    yhist = np.empty((nT, nPart), dtype=np.int64)

    for epoch in range(nT):
        xhist[epoch, :] = r[:, 0]
        yhist[epoch, :] = r[:, 1]
        for iTrial in range(nPart):
            id = np.random.randint(nPart)
            i = r[id,0]
            j = r[id,1]

            empty_nbrs,empty = isempty(i,j,lat,L)

            if empty: 
                
                dE = -J 
                nbrID = np.random.randint(len(empty_nbrs))
                nbr = empty_nbrs[nbrID]
                x,y = nbr 
                lat[i,j] = 0
                lat[x,y] = 1 
                if dE <= 0:
                    r[id,:] = np.array([x,y])
                    # print(id,'before: ', [i,j], 'after: ', r[id,:] )
                else:
                    p = np.random.rand()
                    if p < np.exp(-dE): 
                        r[id,:] = np.array([x,y])
                        # print(id,'before: ', [i,j], 'after: ', r[id,:] )
                    else:
                        lat[i,j] = 1
                        lat[x,y] = 0 
    
    return xhist, yhist


@njit
def get_Rg(x,y,L):
    xc = np.nanmean(x)
    yc = np.nanmean(y)

    xic = x - xc 
    xic = xic - np.round(xic/L)*L
    yic = y - yc 
    yic = yic - np.round(yic/L)*L  

    rg2 = np.nanmean(xic**2 + yic**2)
    return rg2 

@njit
def _rg2_in_cell(xcol, ycol, x_lo, x_hi, y_lo, y_hi, L):
    """Rg^2 for one track, using only points inside the cell (min 6 points)."""
    nT = xcol.shape[0]
    sx = 0.0
    sy = 0.0
    n = 0
    for t in range(nT):
        x = xcol[t]
        y = ycol[t]
        if x >= x_lo and x < x_hi and y >= y_lo and y < y_hi:
            sx += x
            sy += y
            n += 1
    if n < 6:
        return np.nan

    xc = sx / n
    yc = sy / n
    rg2_sum = 0.0
    for t in range(nT):
        x = xcol[t]
        y = ycol[t]
        if x >= x_lo and x < x_hi and y >= y_lo and y < y_hi:
            xic = x - xc
            xic = xic - np.round(xic / L) * L
            yic = y - yc
            yic = yic - np.round(yic / L) * L
            rg2_sum += xic * xic + yic * yic
    return rg2_sum / n


@njit
def _grid_n_cells(L, dx):
    """Match len(np.arange(0, L + 1, dx)) - 1 from the original code."""
    n_xg = 0
    v = 0
    while v <= L:
        n_xg += 1
        v += dx
    return n_xg - 1


@njit
def get_zoomed_RgMat(res, zoom, xhist, yhist, L, nPart):
    dx = zoom // res
    gx = _grid_n_cells(L, dx)
    gy = gx

    xc = np.empty(gx)
    yc = np.empty(gy)
    for i in range(gx):
        xc[i] = i * dx + dx / 2.0
        yc[i] = i * dx + dx / 2.0

    Rg2Mat = np.full((gx, gy), np.nan)
    nT = xhist.shape[0]
    counts = np.zeros((gx, gy, nPart), dtype=np.int32)

    for t in range(nT):
        for k in range(nPart):
            xi = xhist[t, k] // dx
            yi = yhist[t, k] // dx
            if 0 <= xi < gx and 0 <= yi < gy:
                counts[xi, yi, k] += 1

    for i in range(gx):
        x_lo = i * dx
        x_hi = (i + 1) * dx
        for j in range(gy):
            y_lo = j * dx
            y_hi = (j + 1) * dx
            rg_sum = 0.0
            rg_n = 0
            for k in range(nPart):
                if counts[i, j, k] < 6:
                    continue
                rg2 = _rg2_in_cell(
                    xhist[:, k], yhist[:, k], x_lo, x_hi, y_lo, y_hi, L
                )
                if not np.isnan(rg2):
                    rg_sum += rg2
                    rg_n += 1
            if rg_n > 0:
                Rg2Mat[i, j] = rg_sum / rg_n

    return xc, yc, Rg2Mat