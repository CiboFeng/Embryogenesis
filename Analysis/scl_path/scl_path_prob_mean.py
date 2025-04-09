"""Import Modules"""
import numpy as np
import datetime
import multiprocessing
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs
from cont import cont

"""Set Arguments"""
npyfile='../cont_prob/cont_prob.npy'
outname=f'scl_path_prob_mean'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
ntr=[1214,974]
d0=0.6
beta=3
# scl0=20
# scl=[]
# i=0
# while scl0*(2**(i+1)-1)<=nat:
#     scl.append([scl0*(2**i-1),scl0*(2**(i+1)-1)])
#     i+=1
# scl.append([scl0*(2**i-1),nat])
scl=[[0,20],[20,60],[60,140],[140,300],[300,600]]
nsc=len(scl)

"""Calculate for Cells"""
P0=np.load(npyfile)[:,[0,-1]]
Q0=np.zeros((nsw,ncl,nsc))
for i in range(nsc):
    P_temp=[]
    for j in range(nat):
        for k in range(j):
            if abs(j-k) in range(scl[i][0],scl[i][1]):
                P_temp+=[P0[:,:,j,k]]
    Q0[:,:,i]=np.mean(np.array(P_temp),axis=0)

"""Calculate"""
Q=[]
for i in range(nsw):
    Q+=[np.zeros((nfr,ntr[i],nsc))]

    P=np.zeros((nfr,ntr[i],nat,nat))
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        # D+np.linalg.norm(r[:,:,np.newaxis,:]-r[:,np.newaxis,:,:],axis=-1)
        D=np.zeros((ntr[i],nat,nat))
        for k in range(ntr[i]):
            D[k]=squareform(pdist(r[k],'euclidean'))
        P[j]=0.5*(1-np.tanh(beta*(D/10-d0)))

    for j in range(nsc):
        P_temp=[]
        for k in range(nat):
            for l in range(k):
                if abs(k-l) in range(scl[j][0],scl[j][1]):
                    P_temp+=[P[:,:,k,l]]
        Q[i][:,:,j]=np.mean(np.array(P_temp),axis=0)

"""Output Data"""
np.savez(f'{outname}.npz',Q0=Q0,**{f'Q{i+1}':Q[i] for i in range(nsw)})
