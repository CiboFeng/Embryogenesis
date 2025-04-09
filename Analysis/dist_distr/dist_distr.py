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
reffile=[[],[]]
reffile[0]+=['../cont_prob/ref/zp.dat']
reffile[1]+=['../cont_prob/ref/zm.dat']
reffile[0]+=['../cont_prob/ref/e2cp.dat']
reffile[1]+=['../cont_prob/ref/e2cm.dat']
reffile[0]+=['../cont_prob/ref/l2cp.dat']
reffile[1]+=['../cont_prob/ref/l2cm.dat']
reffile[0]+=['../cont_prob/ref/8cp.dat']
reffile[1]+=['../cont_prob/ref/8cm.dat']
reffile[0]+=['../cont_prob/ref/esc.dat']
reffile[1]+=['../cont_prob/ref/esc.dat']
outname='dist_distr'
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
d0=np.zeros((nsw,ncl,nat,nhist))
P_d0=np.zeros((nsw,ncl,nat,nhist))
d0_mean=np.zeros((nsw,ncl,nat))
d0_std=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        r=r[-100000:]
        nfr0=np.shape(r)[0]
        D=np.zeros((nfr0,nat,nat))
        for k in range(nfr0):
            D[k]=squareform(pdist(r[k],'euclidean'))
        for k in range(nat):
            d_temp=np.zeros((nfr0,nat-k))
            for l in range(k,nat):
                d_temp[:,l-k]=D[:,l,l-k]
            hist,edge=np.histogram(d_temp,bins=nhist)
            d0[i,j,k]=(edge[:-1]+edge[1:])/2
            P_d0[i,j,k]=hist/np.sum(hist*(np.max(d0[i,j,k])-np.min(d0[i,j,k]))/nhist)
            d0_mean[i,j,k]=np.mean(d_temp)
            d0_std[i,j,k]=np.std(d_temp)

"""Output Data"""
np.savez(f'{outname}.npz',d0=d0,P_d0=P_d0,d0_mean=d0_mean,d0_std=d0_std)