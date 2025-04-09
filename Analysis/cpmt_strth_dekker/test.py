"""Import Modules"""
import numpy as np
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from cpmt_sgnl import cpmt_sgnl
from pca import pca
from sklearn.decomposition import PCA

"""Set Arguments"""
gtffile='../loci_distr/mm9.knownGene.gtf'
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
outname='test_esc_cor'
cl=['Z','E2C','L2C','8C','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600

"""Read and Calculate the Gene Density"""
f=open(gtffile)
lsf=f.readlines()
Pg=[]
for i in range(len(lsf)):
    ls_f=lsf[i].strip('\n').split()
    if ls_f[0]=='chr15':
        pg=(float(ls_f[3])+float(ls_f[4]))/2/1e5
        if 250<pg<=850:
            Pg+=[pg-251]

bins=np.arange(0-1/2,nat,1)
hist,edge=np.histogram(Pg,bins=bins)
Pg=(edge[:-1]+edge[1:])/2
P_Pg=hist/np.sum(hist)

"""Read Data and Do PCA"""
P0=np.zeros((nat,nat))
f=open(reffile[0][-1],'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P0[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
CS=cpmt_sgnl(P0)[0]
if np.corrcoef(CS,P_Pg)[0,1]<0:
    CS*=-1

P=np.zeros((nat,nat))
f=open(reffile[0][-1],'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
P_nor=cpmt_sgnl(P)[1]

CS=CS[::10]
P_nor=P_nor[::10,::10]
P_nor[P_nor==0]=np.min(P_nor[P_nor>0])/2
P_nor=np.log2(P_nor)
P_nor=np.ma.corrcoef(P_nor)
sort_idx=np.argsort(CS)
CS=CS[sort_idx]
P_nor=P_nor[sort_idx][:,sort_idx]
# P0_nor=P0_nor[sort_idx]

"""Output Data"""
np.savez(f'{outname}.npz',P0_nor=P_nor,CS=CS)