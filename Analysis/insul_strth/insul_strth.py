"""Import Modules"""
import numpy as np
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from insul_score import insul_score

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
outname='insul_strth'
cl=['Z','E2C','L2C','8C','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
scp=5

"""Read Data and Calculate"""
P0=np.zeros((nat,nat))
f=open(reffile[0][-1],'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P0[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
pos=insul_score(P0,bin_size=10**5)[1]
k=1
while k<len(pos):
    if pos[k]-pos[k-1]<=scp or min(pos[k],nat-pos[k])<=scp:
        pos.pop(k)
        k-=1
    k+=1

P_ref=np.zeros((nsw,ncl,nat,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile[i][j],'r')
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            P_ref[i,j,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
ISs_ref=np.zeros((nsw,ncl,2*scp+1))
for i in range(nsw):
    for j in range(ncl):
        IS=insul_score(P_ref[i,j],bin_size=10**5)[0]
        for k in range(len(pos)):
            ISs_ref[i,j]+=IS[pos[k]-scp:pos[k]+scp+1]-np.mean(IS[pos[k]-scp:pos[k]+scp+1])
        ISs_ref[i,j]/=len(pos)

P=np.load(npyfile)
ISs=np.zeros((nsw,nfr,2*scp+1))
for i in range(nsw):
    for j in range(1,nfr+1):
        IS=insul_score(P[i,j],bin_size=10**5)[0]
        for k in range(len(pos)):
            ISs[i,j-1]+=IS[pos[k]-scp:pos[k]+scp+1]-np.mean(IS[pos[k]-scp:pos[k]+scp+1])
        ISs[i,j-1]/=len(pos)

"""Output Data"""
np.savez(f'{outname}.npz',ISs_ref=ISs_ref,ISs=ISs)