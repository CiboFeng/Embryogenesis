"""Import Modules"""
import numpy as np
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from cpmt_sgnl import cpmt_sgnl
from pdbs import pdbs
from sklearn.decomposition import PCA

"""Set Arguments"""
csvfile='../loci_distr/mm9.csv'
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
outname='cpmt_size'
cl_=[['zp','zm'],['esc','esc']]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
nhist=20
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read and Calculate the Gene Density"""
f=open(csvfile)
lsf=f.readlines()
Pg=[]
for i in range(len(lsf)):
    ls_f=lsf[i].strip('\n').split(',')
    if ls_f[1]=='chr15' or ls_f[1]=='"chr15"':
        pg=(float(ls_f[2])+float(ls_f[3]))/2/1e5
        if 250<pg<=850:
            Pg+=[pg-251]

bins=np.arange(0-1/2,nat,1)
hist,edge=np.histogram(Pg,bins=bins)
Pg=(edge[:-1]+edge[1:])/2
P_Pg=hist/np.sum(hist)

"""Define Compartments"""
P0=np.zeros((nat,nat))
f=open(reffile[0][-1],'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P0[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
CS=cpmt_sgnl(P0)[0]
if np.corrcoef(CS,P_Pg)[0,1]<0:
    CS*=-1

"""Calculate for Cells"""
Rg0=np.zeros((nsw,ncl,2,nhist))
P_Rg0=np.zeros((nsw,ncl,2,nhist))
Rg0_mean=np.zeros((nsw,ncl,2))
Rg0_std=np.zeros((nsw,ncl,2))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        nfr0=np.shape(r)[0]
        rg=np.zeros((2,nfr0))
        ra=r[:,CS>0,:]
        ra0=np.mean(ra,axis=1,keepdims=True)
        rg[0]=(np.mean((np.linalg.norm(ra-ra0,axis=-1))**2,axis=1))**0.5
        rb=r[:,CS<0,:]
        rb0=np.mean(rb,axis=1,keepdims=True)
        rg[1]=(np.mean((np.linalg.norm(rb-rb0,axis=-1))**2,axis=1))**0.5
        for k in range(2):
            hist,edge=np.histogram(rg[k],bins=nhist)
            Rg0[i,j,k]=(edge[:-1]+edge[1:])/2
            P_Rg0[i,j,k]=hist/np.sum(hist*(np.max(Rg0[i,j,k])-np.min(Rg0[i,j,k]))/nhist)
            Rg0_mean[i,j,k]=np.mean(rg[k])
            Rg0_std[i,j,k]=np.std(rg[k])

"""Read .pdb Files and Calculate the Radius of Gyration"""
Rg=np.zeros((nsw,nfr,2,nhist))
P_Rg=np.zeros((nsw,nfr,2,nhist))
Rg_mean=np.zeros((nsw,nfr,2))
Rg_std=np.zeros((nsw,nfr,2))
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        ntr=np.shape(r)[0]
        rg=np.zeros((2,ntr))
        ra=r[:,CS>0,:]
        ra0=np.mean(ra,axis=1,keepdims=True)
        rg[0]=(np.mean((np.linalg.norm(ra-ra0,axis=-1))**2,axis=1))**0.5
        rb=r[:,CS<0,:]
        rb0=np.mean(rb,axis=1,keepdims=True)
        rg[1]=(np.mean((np.linalg.norm(rb-rb0,axis=-1))**2,axis=1))**0.5
        for k in range(2):
            hist,edge=np.histogram(rg[k],bins=nhist)
            Rg[i,j,k]=(edge[:-1]+edge[1:])/2
            P_Rg[i,j,k]=hist/np.sum(hist*(np.max(Rg[i,j,k])-np.min(Rg[i,j,k]))/nhist)
            Rg_mean[i,j,k]=np.mean(rg[k])
            Rg_std[i,j,k]=np.std(rg[k])

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
for i in range(nsw):
    pyw.write(f'The Index of Cells, Rg of Compartment A, Its Probability Density, Rg of Compartment B, Its Probability Density (Reference Value, {sw[i]}):\n')
    for j in range(ncl):
        for k in range(nhist):
            pyw.write(f'{j} {Rg0[i,j,0,k]} {P_Rg0[i,j,0,k]} {Rg0[i,j,1,k]} {P_Rg0[i,j,1,k]}\n')
    pyw.write('\n')
    pyw.write(f'The Index of Cells, Mean, and Standard Deviation of Rg of Compartment A, and B (Reference Value, {sw[i]}):\n')
    for j in range(ncl):
        pyw.write(f'{j} {Rg0_mean[i,j,0]} {Rg0_std[i,j,0]} {Rg0_mean[i,j,1]} {Rg0_std[i,j,1]}\n')
    pyw.write('\n')
    pyw.write(f'Time, Rg of Compartment A, Its Probability Density, Rg of Compartment B, Its Probability Density ({sw[i]}):\n')
    for j in range(nfr):
        for k in range(nhist):
            pyw.write(f'{t[j]} {Rg[i,j,0,k]} {P_Rg[i,j,0,k]} {Rg[i,j,1,k]} {P_Rg[i,j,1,k]}\n')
    pyw.write('\n')
    pyw.write(f'Time, Mean, and Standard Deviation of Rg of Compartment A, and B ({sw[i]}):\n')
    for j in range(nfr):
        pyw.write(f'{t[j]} {Rg_mean[i,j,0]} {Rg_std[i,j,0]} {Rg_mean[i,j,1]} {Rg_std[i,j,1]}\n')
    pyw.write('\n')
pyw.close()