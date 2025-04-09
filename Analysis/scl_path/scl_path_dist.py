"""Import Modules"""
import numpy as np
from scipy.spatial.distance import pdist,squareform
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs

"""Set Arguments"""
outname=f'scl_path_dist'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
ntr=[1214,974]
nfr0=40010
dlt=0.2
# scl0=20
# scl=[]
# i=0
# while scl0*(2**(i+1)-1)<=nat:
#     scl.append([scl0*(2**i-1),scl0*(2**(i+1)-1)])
#     i+=1
# scl.append([scl0*(2**i-1),nat])
scl=[[0,20],[20,60],[60,140],[140,300],[300,600]]
nsc=len(scl)

"""Calculate the Reference Distance Matrix"""
_,r=pdbs(f'../trj/esc/trj.pdb')
D0=np.zeros((nat,nat))
for i in range(nfr0):
    D0+=squareform(pdist(r[i],'euclidean'))
D0/=nfr0

"""Calculate for Cells"""
Q0=np.zeros((nsw,ncl,nsc))
for i in range(nsw):
    D=np.zeros((ncl,nfr0,nat,nat))
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        for k in range(nfr0):
            D[j,k]=squareform(pdist(r[k],'euclidean'))

    for j in range(nsc):
        D0_temp=[]
        D_temp=[]
        for k in range(nat):
            for l in range(k):
                if abs(k-l) in range(scl[j][0],scl[j][1]):
                    D0_temp+=[D0[k,l]]
                    D_temp+=[D[:,:,k,l]]
        D0_temp=np.array(D0_temp)
        D_temp=np.array(D_temp)
        Q0[i,:,j]=np.mean(np.exp(-(D_temp-D0_temp[:,np.newaxis,np.newaxis])**2/(2*(10*dlt)**2)),axis=(0,2))
print('OK')
"""Calculate"""
Q=[]
for i in range(nsw):
    Q+=[np.zeros((nfr,ntr[i],nsc))]

    D=np.zeros((nfr,ntr[i],nat,nat))
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        for k in range(ntr[i]):
            D[j,k]=squareform(pdist(r[k],'euclidean'))

    for j in range(nsc):
        D0_temp=[]
        D_temp=[]
        for k in range(nat):
            for l in range(k):
                if abs(k-l) in range(scl[j][0],scl[j][1]):
                    D0_temp+=[D0[k,l]]
                    D_temp+=[D[:,:,k,l]]
        D0_temp=np.array(D0_temp)
        D_temp=np.array(D_temp)
        Q[i][:,:,j]=np.mean(np.exp(-(D_temp-D0_temp[:,np.newaxis,np.newaxis])**2/(2*(10*dlt)**2)),axis=0)

"""Output Data"""
np.savez(f'{outname}.npz',Q0=Q0,**{f'Q{i+1}':Q[i] for i in range(nsw)})
