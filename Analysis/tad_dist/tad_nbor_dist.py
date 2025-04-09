"""Import Modules"""
import numpy as np
import datetime
import multiprocessing
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs
from insul_score import insul_score

"""Set Arguments"""
reffile='../cont_prob/ref/esc.dat'
npyfile='../cont_prob/cont_prob.npy'
outname=f'tad_nbor_dist'
cl_=[['zp','zm'],['esc','esc']]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
nhist=50

"""Define TAD"""
P_ref=np.zeros((nat,nat))
f=open(reffile,'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P_ref[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])

pos=insul_score(P_ref,bin_size=10**5)[1]
k=1
while k<len(pos):
    if pos[k]-pos[k-1]<=1:
        pos.pop(k)
        k-=1
    k+=1
ntad=len(pos)-1

"""Calculate for Cells"""
d0=np.zeros((nsw,ncl,nhist))
P_d0=np.zeros((nsw,ncl,nhist))
d0_mean=np.zeros((nsw,ncl))
d0_std=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        nfr0=np.shape(r)[0]
        r0=np.zeros((nfr0,ntad,3))
        rg=np.zeros((nfr0,ntad))
        for k in range(len(pos)-1):
            r_tad=r[:,pos[k]:pos[k+1],:]
            r0[:,k,:]=np.mean(r_tad,axis=1)
            rg[:,k]=(np.mean((np.linalg.norm(r_tad-r0[:,k:k+1,:],axis=-1))**2,axis=1))**0.5
        d_temp=np.linalg.norm(r0[:,1:,:]-r0[:,:-1,:],axis=-1)-(rg[:,1:]+rg[:,:-1])
        hist,edge=np.histogram(d_temp,bins=nhist)
        d0[i,j]=(edge[:-1]+edge[1:])/2
        P_d0[i,j]=hist/np.sum(hist*(np.max(d0[i,j])-np.min(d0[i,j]))/nhist)
        d0_mean[i,j]=np.mean(d_temp)
        d0_std[i,j]=np.std(d_temp)

"""Read .pdb Files and Calculate the Radius of Gyration"""
d=np.zeros((nsw,nfr,nhist))
P_d=np.zeros((nsw,nfr,nhist))
d_mean=np.zeros((nsw,nfr))
d_std=np.zeros((nsw,nfr))
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        ntr=np.shape(r)[0]
        r0=np.zeros((ntr,ntad,3))
        rg=np.zeros((ntr,ntad))
        for k in range(len(pos)-1):
            r_tad=r[:,pos[k]:pos[k+1],:]
            r0[:,k,:]=np.mean(r_tad,axis=1)
            rg[:,k]=(np.mean((np.linalg.norm(r_tad-r0[:,k:k+1,:],axis=-1))**2,axis=1))**0.5
        d_temp=np.linalg.norm(r0[:,1:,:]-r0[:,:-1,:],axis=-1)-(rg[:,1:]+rg[:,:-1])
        hist,edge=np.histogram(d_temp,bins=nhist)
        d[i,j]=(edge[:-1]+edge[1:])/2
        P_d[i,j]=hist/np.sum(hist*(np.max(d[i,j])-np.min(d[i,j]))/nhist)
        d_mean[i,j]=np.mean(d_temp)
        d_std[i,j]=np.std(d_temp)

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
for i in range(nsw):
    pyw.write(f'The Index of Cells, Neighbor Distance between TAD, Its Probability Density (Reference Value, {sw[i]}):\n')
    for j in range(ncl):
        for k in range(nhist):
            pyw.write(f'{j} {d0[i,j,k]} {P_d0[i,j,k]}\n')
    pyw.write('\n')
    pyw.write(f'The Index of Cells, Mean, and Standard Deviation of Neighbor Distance between TAD (Reference Value, {sw[i]}):\n')
    for j in range(ncl):
        pyw.write(f'{j} {d0_mean[i,j]} {d0_std[i,j]}\n')
    pyw.write('\n')
    pyw.write(f'Time, Neighbor Distance between TAD, Its Probability Density ({sw[i]}):\n')
    for j in range(nfr):
        for k in range(nhist):
            pyw.write(f'{t[j]} {d[i,j,k]} {P_d[i,j,k]}\n')
    pyw.write('\n')
    pyw.write(f'Time, Mean, and Standard Deviation of Neighbor Distance between TAD ({sw[i]}):\n')
    for j in range(nfr):
        pyw.write(f'{t[j]} {d_mean[i,j]} {d_std[i,j]}\n')
    pyw.write('\n')
pyw.close()
