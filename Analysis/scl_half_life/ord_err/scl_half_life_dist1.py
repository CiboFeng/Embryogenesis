"""Import Modules"""
import numpy as np
from scipy.spatial.distance import pdist,squareform
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs

"""Set Arguments"""
npyfile='../cont_prob/cont_prob.npy'
outname=f'scl_half_life_dist'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
nfr0=400100
ntr=[1214,974]
dlt=0.2
scl=[[0,20],[20,60],[60,140],[140,300]]
nsc=len(scl)
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Calculate the Reference Distance Matrix"""
_,r=pdbs(f'../trj/esc/trj.pdb')
D0=np.zeros((nat,nat))
for i in range(nfr0):
    D0+=squareform(pdist(r[i],'euclidean'))
D0/=nfr0

"""Calculate for Cells"""
D=np.zeros((nsw,ncl,nat,nat))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        for k in range(nfr0):
            D[i,j]+=squareform(pdist(r[k],'euclidean'))
        D[i,j]/=nfr0

Q0=np.zeros((nsw,ncl,nsc,nat))
for i in range(nsc):
    D0_temp=np.zeros((scl[i][1]-scl[i][0],nat))
    D_temp=np.zeros((nsw,ncl,scl[i][1]-scl[i][0],nat))
    for j in range(nat-scl[i][1]+1):
        D0_temp[:,j]=D0[j+scl[i][0]:j+scl[i][1],j]
        D_temp[:,:,:,j]=D[:,:,j+scl[i][0]:j+scl[i][1],j]
    for j in range(nat-scl[i][1]+1,nat):
        D0_temp[:,j]=D0[j-scl[i][1]:j-scl[i][0],j]
        D_temp[:,:,:,j]=D[:,:,j-scl[i][1]:j-scl[i][0],j]
    Q0[:,:,i,:]=np.mean(np.exp(-(D_temp-D0_temp[np.newaxis,np.newaxis,:,:])**2/(2*(10*dlt)**2)),axis=-2)

"""Calculate"""
D=np.zeros((nsw,nfr,nat,nat))
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        for k in range(ntr[i]):
            D[i,j]+=squareform(pdist(r[k],'euclidean'))
        D[i,j]/=ntr[i]

Q=np.zeros((nsw,nfr,nsc,nat))
for i in range(nsc):
    D0_temp=np.zeros((scl[i][1]-scl[i][0],nat))
    D_temp=np.zeros((nsw,nfr,scl[i][1]-scl[i][0],nat))
    for j in range(nat-scl[i][1]+1):
        D0_temp[:,j]=D0[j+scl[i][0]:j+scl[i][1],j]
        D_temp[:,:,:,j]=D[:,:,j+scl[i][0]:j+scl[i][1],j]
    for j in range(nat-scl[i][1]+1,nat):
        D0_temp[:,j]=D0[j-scl[i][1]:j-scl[i][0],j]
        D_temp[:,:,:,j]=D[:,:,j-scl[i][1]:j-scl[i][0],j]
    Q[:,:,i,:]=np.mean(np.exp(-(D_temp-D0_temp[np.newaxis,np.newaxis,:,:])**2/(2*(10*dlt)**2)),axis=-2)

"""Output Data"""
np.savez(f'{outname}.npz',Q0=Q0,Q=Q)
