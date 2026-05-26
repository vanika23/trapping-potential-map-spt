import numpy as np 
import matplotlib.pyplot as plt 
from numba import njit

# from scipy.signal import periodogram, find_peaks, find_peaks_cwt
import os

hfont = {'fontname':'Helvetica'}
import matplotlib.font_manager as font_manager

from MC_helper import *

# Initialize in clusters 

prefix ='Uniform_J_fast/'
if not os.path.exists(prefix):
    os.makedirs(prefix)

# @njit
def main(): 
    for seed in range(1,11):
        
        # plt.close('all')

        nPart = 200
        # nCluster = 4 
        # pPerCluster = 9
        # inCluster = nCluster*pPerCluster 

        np.random.seed(seed)
        L = 960 
        r = np.zeros((nPart,2)).astype(int)

        lat = np.zeros((L,L))

        for i in range(nPart):
            success = False
            while not success:
                x = int(np.random.rand()*L)
                y = int(np.random.rand()*L)

                if lat[x,y] == 0:
                    lat[x,y] = 1
                    r[i,:] = np.array([x,y])
                    success = True 
                else:
                    success = False 
        # fig,ax = plt.subplots()
        # ax.imshow(lat)


        ############ run MC moves ############
        for J in [0,-2,-3]:

            nEpoch = 10000 
            xhist, yhist = run_mc(r, lat, J, nEpoch)

            # fig,ax = plt.subplots(figsize=(5,5))
            # for i in range(nPart):
            #     ax.plot(xhist[:,i],yhist[:,i],'.-',ms=1,lw=1,alpha=0.4)

            print (xhist.shape)

            exptParam ='J=%d_seed_%d'%(J,seed)
            filename = prefix + exptParam +'.npz'
            np.savez(filename,xhist=xhist,yhist=yhist,L=L,nEpoch=nEpoch,J=J,seed=seed,nPart=nPart)


            # Grids of 160 nm 

            for zoom in [160, 80, 40, 20]:
                print(seed, J, zoom)
                res = 5  # nm
                xc, yc, Rg2Mat = get_zoomed_RgMat(res, zoom, xhist, yhist, L, nPart)
                filename = prefix + exptParam + '_zoom_%dx.npz'%(160/zoom)
                np.savez(filename,xc=xc,yc=yc,Rg2Mat=Rg2Mat)
            

main()