"""Import Modules"""
import numpy as np
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from dist_prob import dist_prob

"""Set Arguments"""
npyfile='../cont_prob/cont_prob.npy'
reffile=[]
reffile+=['../cont_prob/ref/zp.dat']
reffile+=['../cont_prob/ref/zm.dat']
reffile+=['../cont_prob/ref/e2cp.dat']
reffile+=['../cont_prob/ref/e2cm.dat']
reffile+=['../cont_prob/ref/l2cp.dat']
reffile+=['../cont_prob/ref/l2cm.dat']
reffile+=['../cont_prob/ref/8cp.dat']
reffile+=['../cont_prob/ref/8cm.dat']
reffile+=['../cont_prob/ref/esc.dat']
outname='r_cont_prob'
split_label='/'
filename=npyfile.split(split_label)[-1].split('.')[0]
# cl=['ZP','ZM','ESC']
cl=['ZP','ZM','E2CP','E2CM','L2CP','L2CM','8CP','8CM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
# idx_cl=[(0,0),(1,0),(0,-1)]
idx_cl=[-9,-8,-7,-6,-5,-4,-3,-2,-1]

"""Read Data and Calculate"""
P=np.load(npyfile)
# for i in range(nsw):
#     for j in range(nfr+2):
#         dp=dist_prob(P[i,j])
#         norm=np.zeros((nat,nat))
#         for k in range(nat):
#             for l in range(nat):
#                 norm[k,l]=dp[abs(k-l)]
#         P[i,j]-=norm
P=P.reshape(nsw,nfr+2,-1)

P_ref=np.zeros((len(reffile),nat,nat))
for i in range(len(reffile)):
    f=open(reffile[i],'r')
    lsf=f.readlines()
    for j in range(len(lsf)):
        ls_f=lsf[j].strip('\n').split()
        P_ref[i,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
# for i in range(len(reffile)):
#     dp=dist_prob(P_ref[i])
#     norm=np.zeros((nat,nat))
#     for j in range(nat):
#         for k in range(nat):
#             norm[j,k]=dp[abs(j-k)]
#     P_ref[i]-=norm
P_ref=P_ref.reshape(len(reffile),-1)

Pl0=np.zeros((ncl,nat**2))
for i in range(ncl):
    Pl0[i]=P_ref[idx_cl[i]]
Pl=np.zeros((nsw,nfr,nat**2))
for i in range(nsw):
    Pl[i]=P[i,1:nfr+1]

r=np.zeros((ncl,nsw,nfr))
for i in range(ncl):
    for j in range(nsw):
        for k in range(nfr):
            r[i,j,k]=np.corrcoef(Pl[j,k],Pl0[i])[0,1]

"""Output Data"""
np.save(f'{outname}.npy',r)


