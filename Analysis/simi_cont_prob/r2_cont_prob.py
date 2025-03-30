"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from dist_prob import dist_prob

"""Set Arguments"""
npyfile='../cont_prob/cont_prob.npy'
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
outname='r2_cont_prob'
split_label='/'
filename=npyfile.split(split_label)[-1].split('.')[0]
cl=['Z','E2C','L2C','8C','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600

"""Read Data and Calculate"""
P_ref=np.zeros((nsw,ncl,nat,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile[i][j],'r')
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            P_ref[i,j,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
# for i in range(nsw):
#     for j in range(ncl):
#         dp=dist_prob(P_ref[i,j])
#         norm=np.zeros((nat,nat))
#         for k in range(nat):
#             for l in range(nat):
#                 norm[k,l]=dp[abs(k-l)]
#         P_ref[i,j]-=norm

P=np.load(npyfile)[:,1:nfr+1,:,:]
# for i in range(nsw):
#     for j in range(nfr+2):
#         dp=dist_prob(P[i,j])
#         norm=np.zeros((nat,nat))
#         for k in range(nat):
#             for l in range(nat):
#                 norm[k,l]=dp[abs(k-l)]
#         P[i,j]-=norm

P_ref=P_ref.reshape(nsw,ncl,-1)
P=P.reshape(nsw,nfr,-1)

r2=np.zeros((nsw,ncl,nsw,nfr))
for i in range(nsw):
    for j in range(ncl):
        for k in range(nsw):
            for l in range(nfr):
                r2[i,j,k,l]=r2_score(P[k,l],P_ref[i,j])

"""Output Data"""
np.save(f'{outname}.npy',r2)


