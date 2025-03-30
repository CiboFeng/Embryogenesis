"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from dist_prob import dist_prob

"""Set Arguments"""
npyfile='../cont_prob/cont_prob.npy'
outname='r2_cont_prob_cross'
split_label='/'
filename=npyfile.split(split_label)[-1].split('.')[0]
cl=['ZP','ZM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
idx_cl=[(0,0),(1,0),(0,-1)]

"""Read Data and Calculate"""
P=np.load(npyfile)
for i in range(nsw):
    for j in range(nfr+2):
        dp=dist_prob(P[i,j])
        norm=np.zeros((nat,nat))
        for k in range(nat):
            for l in range(nat):
                norm[k,l]=dp[abs(k-l)]
        P[i,j]-=norm
P=P.reshape(nsw,nfr+2,-1)
Pl=np.zeros((nsw,nfr,nat**2))
for i in range(nsw):
    Pl[i]=P[i,1:nfr+1]

r2=np.zeros((nsw,nsw,nfr,nfr))
for i in range(nsw):
    for j in range(nsw):
        for k in range(nfr):
            for l in range(nfr):
                r2[i,j,k,l]=r2_score(Pl[i,k],Pl[j,l])

"""Output Data"""
np.save(f'{outname}.npy',r2)


