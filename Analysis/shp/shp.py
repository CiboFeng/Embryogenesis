"""Import Modules"""
import numpy as np
import multiprocessing
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from pdbs import pdbs
from sklearn.decomposition import PCA

"""Set Arguments"""
outname='shp'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
ntr=[1214,974]
nfr=56
nat=600

"""Calculate for Cells"""
PA0=np.zeros((nsw,ncl,3))
Rg0=np.zeros((nsw,ncl))
Dlt0=np.zeros((nsw,ncl))
S0=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        nfr0=np.shape(r)[0]
        r-=np.mean(r,axis=1,keepdims=True)
        IT=np.sum(np.einsum('ijk,ijl->ijkl',r,r),axis=1)/nat
        tr=np.trace(IT,axis1=1,axis2=2)
        lmd=np.array([np.sort(np.linalg.eigvals(IT[k]))[::-1] for k in range(nfr0)])
        PA0[i,j]=np.mean(lmd**0.5,axis=0)
        Rg0[i,j]=np.mean(tr**0.5)
        Dlt0[i,j]=np.mean(9/2*np.var(lmd,axis=1)/tr**2)
        S0[i,j]=np.mean(27*np.product(lmd-np.mean(lmd,axis=1,keepdims=True),axis=1)/tr**3)

"""Read .pdb Files and Calculate"""
PA=[]
Rg=[]
Dlt=[]
S=[]
for i in range(nsw):
    PA+=[np.zeros((nfr,ntr[i],3))]
    Rg+=[np.zeros((nfr,ntr[i]))]
    Dlt+=[np.zeros((nfr,ntr[i]))]
    S+=[np.zeros((nfr,ntr[i]))]
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        r-=np.mean(r,axis=1,keepdims=True)
        IT=np.sum(np.einsum('ijk,ijl->ijkl',r,r),axis=1)/nat
        tr=np.trace(IT,axis1=1,axis2=2)
        lmd=np.array([np.sort(np.linalg.eigvals(IT[k]))[::-1] for k in range(ntr[i])])
        PA[i][j]=lmd**0.5
        Rg[i][j]=tr**0.5
        Dlt[i][j]=9/2*np.var(lmd,axis=1)/tr**2
        S[i][j]=27*np.product(lmd-np.mean(lmd,axis=1,keepdims=True),axis=1)/tr**3

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
for i in range(nsw):
    pyw.write(f'The Index of Cells, the Length of Principal Axis 1, 2, 3, Radius of Gyration, Aspherical Quantity, and '
              f'S ({sw[i]}):\n')
    for j in range(ncl):
        pyw.write(f'{j} {PA0[i,j,0]} {PA0[i,j,1]} {PA0[i,j,2]} {Rg0[i,j]} {Dlt0[i,j]} {S0[i,j]}\n')
    pyw.write('\n')
    pyw.write(f'The Index of Frames, The Index of Trajectories, the Length of Principal Axis 1, 2, 3, Radius of '
              f'Gyration, Aspherical Quantity, and S ({sw[i]}):\n')
    for j in range(nfr):
        for k in range(ntr[i]):
            pyw.write(f'{j} {k} {PA[i][j,k,0]} {PA[i][j,k,1]} {PA[i][j,k,2]} {Rg[i][j,k]} {Dlt[i][j,k]} {S[i][j,k]}\n')
    pyw.write('\n')



