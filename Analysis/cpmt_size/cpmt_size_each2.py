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
outname='cpmt_size_each2'
cl_=[['zp','zm'],['esc','esc']]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
nhist=50
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
da0_mean=np.zeros((nsw,ncl,len(idxa),lmax))
da0_std=np.zeros((nsw,ncl,len(idxa),lmax))
db0_mean=np.zeros((nsw,ncl,len(idxb),lmax))
db0_std=np.zeros((nsw,ncl,len(idxb),lmax))
for i in range(nsw):
    for j in range(ncl):
        _,r=pdbs(f'../trj/{cl_[j][i]}/trj.pdb')
        nfr0=np.shape(r)[0]
        for k in range(len(idxa)):
            r_=r[:,idxa[k][0]:idxa[k][1],:]
            nat_=idxa[k][1]-idxa[k][0]
            D=np.zeros((nfr0,nat_,nat_))
            for l in range(nfr0):
                D[l]=squareform(pdist(r_[l],'euclidean'))
            for l in range(nat_):
                d_temp=np.zeros((nfr0,nat_-l))
                for m in range(l,nat_):
                    d_temp[:,m-l]=D[:,m,m-l]
                da0_mean[i,j,k,l]=np.mean(d_temp)
                da0_std[i,j,k,l]=np.std(d_temp)
        for k in range(len(idxb)):
            r_=r[:,idxb[k][0]:idxb[k][1],:]
            nat_=idxb[k][1]-idxb[k][0]
            D=np.zeros((nfr0,nat_,nat_))
            for l in range(nfr0):
                D[l]=squareform(pdist(r_[l],'euclidean'))
            for l in range(nat_):
                d_temp=np.zeros((nfr0,nat_-l))
                for m in range(l,nat_):
                    d_temp[:,m-l]=D[:,m,m-l]
                db0_mean[i,j,k,l]=np.mean(d_temp)
                db0_std[i,j,k,l]=np.std(d_temp)

"""Read .pdb Files and Calculate the Radius of Gyration"""
da_mean=np.zeros((nsw,nfr,len(idxa),lmax))
da_std=np.zeros((nsw,nfr,len(idxa),lmax))
db_mean=np.zeros((nsw,nfr,len(idxb),lmax))
db_std=np.zeros((nsw,nfr,len(idxb),lmax))
for i in range(nsw):
    for j in range(nfr):
        _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
        ntr=np.shape(r)[0]
        for k in range(len(idxa)):
            r_=r[:,idxa[k][0]:idxa[k][1],:]
            nat_=idxa[k][1]-idxa[k][0]
            D=np.zeros((ntr,nat_,nat_))
            for l in range(ntr):
                D[l]=squareform(pdist(r_[l],'euclidean'))
            for l in range(nat_):
                d_temp=np.zeros((ntr,nat_-l))
                for m in range(l,nat_):
                    d_temp[:,m-l]=D[:,m,m-l]
                da_mean[i,j,k,l]=np.mean(d_temp)
                da_std[i,j,k,l]=np.std(d_temp)
        for k in range(len(idxb)):
            r_=r[:,idxb[k][0]:idxb[k][1],:]
            nat_=idxb[k][1]-idxb[k][0]
            D=np.zeros((ntr,nat_,nat_))
            for l in range(ntr):
                D[l]=squareform(pdist(r_[l],'euclidean'))
            for l in range(nat_):
                d_temp=np.zeros((ntr,nat_-l))
                for m in range(l,nat_):
                    d_temp[:,m-l]=D[:,m,m-l]
                db_mean[i,j,k,l]=np.mean(d_temp)
                db_std[i,j,k,l]=np.std(d_temp)

"""Output Data"""
np.savez(f'{outname}.npz',da0_mean=da0_mean,da0_std=da0_std,db0_mean=db0_mean,db0_std=db0_std,
         da_mean=da_mean,da_std=da_std,db_mean=db_mean,db_std=db_std)