"""Import Modules"""
import numpy as np
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from insul_score import insul_score

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
outname='r_insul_score'
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
IS=np.zeros((nsw,nfr,nat))
for i in range(nsw):
    for j in range(1,nfr+1):
        IS[i,j-1]=insul_score(P[i,j],10**5)[0]

P_ref=np.zeros((ncl,nat,nat))
for i in range(ncl):
    f=open(reffile[i],'r')
    lsf=f.readlines()
    for j in range(len(lsf)):
        ls_f=lsf[j].strip('\n').split()
        P_ref[i,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
IS_ref=np.zeros((ncl,nat))
for i in range(ncl):
    IS_ref[i]=insul_score(P_ref[i],10**5)[0]

r2=np.zeros((ncl,nsw,nfr))
for i in range(ncl):
    for j in range(nsw):
        for k in range(nfr):
            r2[i,j,k]=np.corrcoef(IS[j,k],IS_ref[i])[0,1]

"""Output Data"""
np.save(f'{outname}.npy',r2)


