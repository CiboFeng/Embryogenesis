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
outname='dist_cont_prob'
split_label='/'
filename=npyfile.split(split_label)[-1].split('.')[0]
# cl=['ZP','ZM','ESC']
cl=['ZP','ZM','E2CP','E2CM','L2CP','L2CM','8CP','8CM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600

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
P=P[:,1:nfr+1,:,:]

P_ref=np.zeros((ncl,nat,nat))
for i in range(ncl):
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

d=np.zeros((ncl,nsw,nfr))
for i in range(ncl):
    for j in range(nsw):
        for k in range(nfr):
            d[i,j,k]=(np.mean((P[j,k]-P_ref[i])**2))**0.5

"""Output Data"""
np.save(f'{outname}.npy',d)


