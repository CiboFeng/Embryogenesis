"""Import Modules"""
import numpy as np
from scipy.spatial.distance import pdist,squareform
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
outname='cpmt_size_each1'
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

idxa=[[],[]]
idxb=[[],[]]
if CS[0]>=0:
    idxa[0]+=[0]
else:
    idxb[0]+=[0]
for i in range(1,nat):
    if CS[i-1]<0 and CS[i]>0:
        idxa[0]+=[i]
        idxb[1]+=[i]
    if CS[i-1]>0 and CS[i]<0:
        idxa[1]+=[i]
        idxb[0]+=[i]
if CS[-1]>=0:
    idxa[1]+=[nat]
else:
    idxb[1]+=[nat]
idxa=np.array(idxa).T.tolist()
idxb=np.array(idxb).T.tolist()

lmax=[]
for i in range(len(idxa)):
    lmax+=[idxa[i][1]-idxa[i][0]]
for i in range(len(idxb)):
    lmax+=[idxb[i][1]-idxb[i][0]]
lmax=max(lmax)

"""Calculate for Cells"""
Rga0=np.zeros((nsw,ncl,len(idxa),nhist))
P_Rga0=np.zeros((nsw,ncl,len(idxa),nhist))
Rga0_mean=np.zeros((nsw,ncl,len(idxa)))
Rga0_std=np.zeros((nsw,ncl,len(idxa)))
Rgb0=np.zeros((nsw,ncl,len(idxb),nhist))
P_Rgb0=np.zeros((nsw,ncl,len(idxb),nhist))
Rgb0_mean=np.zeros((nsw,ncl,len(idxb)))
Rgb0_std=np.zeros((nsw,ncl,len(idxb)))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        for k in range(len(idxa)):
            ra=r[:,idxa[k][0]:idxa[k][1],:]
            ra0=np.mean(ra,axis=1,keepdims=True)
            rga=(np.mean((np.linalg.norm(ra-ra0,axis=-1))**2,axis=1))**0.5
            hist,edge=np.histogram(rga,bins=nhist)
            Rga0[i,j,k]=(edge[:-1]+edge[1:])/2
            P_Rga0[i,j,k]=hist/np.sum(hist*(np.max(Rga0[i,j,k])-np.min(Rga0[i,j,k]))/nhist)
            Rga0_mean[i,j,k]=np.mean(rga)
            Rga0_std[i,j,k]=np.std(rga)
        for k in range(len(idxb)):
            rb=r[:,idxb[k][0]:idxb[k][1],:]
            rb0=np.mean(rb,axis=1,keepdims=True)
            rgb=(np.mean((np.linalg.norm(rb-rb0,axis=-1))**2,axis=1))**0.5
            hist,edge=np.histogram(rgb,bins=nhist)
            Rgb0[i,j,k]=(edge[:-1]+edge[1:])/2
            P_Rgb0[i,j,k]=hist/np.sum(hist*(np.max(Rgb0[i,j,k])-np.min(Rgb0[i,j,k]))/nhist)
            Rgb0_mean[i,j,k]=np.mean(rgb)
            Rgb0_std[i,j,k]=np.std(rgb)

"""Read .pdb Files and Calculate the Radius of Gyration"""
Rga=np.zeros((nsw,nfr,len(idxa),nhist))
P_Rga=np.zeros((nsw,nfr,len(idxa),nhist))
Rga_mean=np.zeros((nsw,nfr,len(idxa)))
Rga_std=np.zeros((nsw,nfr,len(idxa)))
Rgb=np.zeros((nsw,nfr,len(idxb),nhist))
P_Rgb=np.zeros((nsw,nfr,len(idxb),nhist))
Rgb_mean=np.zeros((nsw,nfr,len(idxb)))
Rgb_std=np.zeros((nsw,nfr,len(idxb)))
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        for k in range(len(idxa)):
            ra=r[:,idxa[k][0]:idxa[k][1],:]
            ra0=np.mean(ra,axis=1,keepdims=True)
            rga=(np.mean((np.linalg.norm(ra-ra0,axis=-1))**2,axis=1))**0.5
            hist,edge=np.histogram(rga,bins=nhist)
            Rga[i,j,k]=(edge[:-1]+edge[1:])/2
            P_Rga[i,j,k]=hist/np.sum(hist*(np.max(Rga[i,j,k])-np.min(Rga[i,j,k]))/nhist)
            Rga_mean[i,j,k]=np.mean(rga)
            Rga_std[i,j,k]=np.std(rga)
        for k in range(len(idxb)):
            rb=r[:,idxb[k][0]:idxb[k][1],:]
            rb0=np.mean(rb,axis=1,keepdims=True)
            rgb=(np.mean((np.linalg.norm(rb-rb0,axis=-1))**2,axis=1))**0.5
            hist,edge=np.histogram(rgb,bins=nhist)
            Rgb[i,j,k]=(edge[:-1]+edge[1:])/2
            P_Rgb[i,j,k]=hist/np.sum(hist*(np.max(Rgb[i,j,k])-np.min(Rgb[i,j,k]))/nhist)
            Rgb_mean[i,j,k]=np.mean(rgb)
            Rgb_std[i,j,k]=np.std(rgb)

"""Output Data"""
np.savez(f'{outname}.npz',Rga0=Rga0,P_Rga0=P_Rga0,Rgb0=Rgb0,P_Rgb0=P_Rgb0,
         Rga0_mean=Rga0_mean,Rga0_std=Rga0_std,Rgb0_mean=Rgb0_mean,Rgb0_std=Rgb0_std,
         Rga=Rga,P_Rga=P_Rga,Rgb=Rgb,P_Rgb=P_Rgb,
         Rga_mean=Rga_mean,Rga_std=Rga_std,Rgb_mean=Rgb_mean,Rgb_std=Rgb_std)