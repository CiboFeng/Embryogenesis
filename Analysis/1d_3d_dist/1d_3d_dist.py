"""Import Modules"""
import numpy as np
from scipy.spatial.distance import pdist,squareform
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from cpmt_sgnl import cpmt_sgnl
from pdbs import pdbs
from sklearn.decomposition import PCA

"""Set Arguments"""
outname='1d_3d_dist'
cl_=[['zp','zm'],['esc','esc']]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
nhist=50
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Calculate for Cells"""
d0_mean=np.zeros((nsw,ncl,nat))
d0_std=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        nfr0=np.shape(r)[0]
        D=np.zeros((nfr0,nat,nat))
        for k in range(nfr0):
            D[k]=squareform(pdist(r[k],'euclidean'))
        for k in range(nat):
            d_temp=np.zeros((nfr0,nat-k))
            for l in range(k,nat):
                d_temp[:,l-k]=D[:,l,l-k]
            d0_mean[i,j,k]=np.mean(d_temp)
            d0_std[i,j,k]=np.std(d_temp)

"""Calculate"""
d_mean=np.zeros((nsw,nfr,nat))
d_std=np.zeros((nsw,nfr,nat))
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        ntr=np.shape(r)[0]
        D=np.zeros((ntr,nat,nat))
        for k in range(ntr):
            D[k]=squareform(pdist(r[k],'euclidean'))
        for k in range(nat):
            d_temp=np.zeros((ntr,nat-k))
            for l in range(k,nat):
                d_temp[:,l-k]=D[:,l,l-k]
            d_mean[i,j,k]=np.mean(d_temp)
            d_std[i,j,k]=np.std(d_temp)

"""Output Data"""
np.savez(f'{outname}.npz',d0_mean=d0_mean,d0_std=d0_std,d_mean=d_mean,d_std=d_std)