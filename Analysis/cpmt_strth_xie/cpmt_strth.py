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
outname='cpmt_strth'
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

"""Read Data and Calculate"""
P0=np.zeros((nat,nat))
f=open(reffile[0][-1],'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P0[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
CS=cpmt_sgnl(P0)[0]
if np.corrcoef(CS,P_Pg)[0,1]<0:
    CS*=-1
srt_idx=np.argsort(CS)
n1=np.sum(CS<0)
n2=np.sum(CS>0)

P_ref=np.zeros((nsw,ncl,nat,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile[i][j],'r')
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            P_ref[i,j,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
for i in range(nsw):
    for j in range(ncl):
        for k in range(nat):
            for l in range(nat):
                if abs(k-l)<20:
                    P_ref[i,j,k,l]=np.nan
CSs_ref=np.zeros((nsw,ncl,3))
for i in range(nsw):
    for j in range(ncl):
        P_srt=P_ref[i,j][srt_idx][:,srt_idx]
        CSs_ref[i,j,0]=np.nanmean(P_srt[:n1,:n1])
        CSs_ref[i,j,1]=np.nanmean(P_srt[-n2:,-n2:])
        CSs_ref[i,j,2]=np.nanmean(P_srt[-n2:,:n1])

P=np.load(npyfile)
for i in range(nsw):
    for j in range(nfr):
        for k in range(nat):
            for l in range(nat):
                if abs(k-l)<20:
                    P[i,j,k,l]=np.nan
CSs=np.zeros((nsw,nfr,3))
for i in range(nsw):
    for j in range(1,nfr+1):
        P_srt=P[i,j][srt_idx][:,srt_idx]
        CSs[i,j-1,0]=np.nanmean(P_srt[:n1,:n1])
        CSs[i,j-1,1]=np.nanmean(P_srt[-n2:,-n2:])
        CSs[i,j-1,2]=np.nanmean(P_srt[-n2:,:n1])

"""Output Data"""
np.savez(f'{outname}.npz',CSs_ref=CSs_ref,CSs=CSs)
