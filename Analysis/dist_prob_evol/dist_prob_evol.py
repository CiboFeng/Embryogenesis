"""Import Modules"""
import numpy as np
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
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
outname='dist_prob'
cl=['Z','E2C','L2C','8C','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600

"""Calculate"""
P_ref=np.zeros((nsw,ncl,nat,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile[i][j],'r')
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            P_ref[i,j,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
DP_ref=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    for j in range(ncl):
        DP_ref[i,j]=dist_prob(P_ref[i,j])

P=np.load(npyfile)
DP=np.zeros((nsw,nfr,nat))
for i in range(nsw):
    for j in range(1,nfr+1):
        DP[i,j-1]=dist_prob(P[i,j])

"""Output Data"""
np.savez(f'{outname}.npz',DP_ref=DP_ref,DP=DP)