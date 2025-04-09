"""Import Modules"""
import numpy as np
import multiprocessing
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs
from cpmt_sgnl import cpmt_sgnl

"""Set Arguments"""
# gtffile='mm9.knownGene.gtf'
csvfile='mm9.csv'
npyfile='../cont_prob/cont_prob.npy'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
nsw=len(sw)
ncl=len(cl)
outname=f'loci_distr'

nat=600
nfr=56
ntr=[1214,974]
nhist=50
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
rmax=52.3424

"""Read and Calculate the Gene Density"""
# f=open(gtffile)
# lsf=f.readlines()
# Pg=[]
# for i in range(len(lsf)):
#     ls_f=lsf[i].strip('\n').split()
#     if ls_f[0]=='chr15':
#         pg=(float(ls_f[3])+float(ls_f[4]))/2/1e5
#         if 250<pg<=850:
#             Pg+=[pg-251]

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

"""Define Compartment State"""
P0=np.load(npyfile)[0,-1]
cs=cpmt_sgnl(P0)[0]
if np.corrcoef(cs,P_Pg)[0,1]<0:
    cs*=-1
state=np.heaviside(cs,1)

"""Read .pdb Files for Cells and Calculate Distance to Center"""
r=[]
for i in range(nsw):
    r.append([])
    for j in range(ncl):
        r[-1]+=[pdbs(f'../trj/{cl_[j][i]}/trj.pdb')[1]]
r=np.array(r)

r0=np.mean(r,axis=-2,keepdims=True)
r=np.linalg.norm(r-r0,axis=-1)

"""Collect and Analysis"""
idx_a=np.where(state==1)[0]
idx_b=np.where(state==0)[0]
r_a=r[:,:,:,idx_a]
r_b=r[:,:,:,idx_b]

dr=rmax/(nhist-1)
bins=np.arange(0-dr/2,rmax+dr,dr)
R0_all=np.zeros((nsw,ncl,nhist))
P0_all=np.zeros((nsw,ncl,nhist))
R0_mean_all=np.zeros((nsw,ncl))
R0_std_all=np.zeros((nsw,ncl))
R0_a=np.zeros((nsw,ncl,nhist))
P0_a=np.zeros((nsw,ncl,nhist))
R0_mean_a=np.zeros((nsw,ncl))
R0_std_a=np.zeros((nsw,ncl))
R0_a_=np.zeros((nsw,ncl,nhist))
P0_a_=np.zeros((nsw,ncl,nhist))
R0_mean_a_=np.zeros((nsw,ncl))
R0_std_a_=np.zeros((nsw,ncl))
R0_b=np.zeros((nsw,ncl,nhist))
P0_b=np.zeros((nsw,ncl,nhist))
R0_mean_b=np.zeros((nsw,ncl))
R0_std_b=np.zeros((nsw,ncl))
R0_b_=np.zeros((nsw,ncl,nhist))
P0_b_=np.zeros((nsw,ncl,nhist))
R0_mean_b_=np.zeros((nsw,ncl))
R0_std_b_=np.zeros((nsw,ncl))
R0_g=np.zeros((nsw,ncl,nhist))
P0_g=np.zeros((nsw,ncl,nhist))
R0_mean_g=np.zeros((nsw,ncl))
R0_std_g=np.zeros((nsw,ncl))
R0_g_=np.zeros((nsw,ncl,nhist))
P0_g_=np.zeros((nsw,ncl,nhist))
R0_mean_g_=np.zeros((nsw,ncl))
R0_std_g_=np.zeros((nsw,ncl))
w=np.zeros((nsw,ncl,np.shape(r)[2],nat))
for i in range(nsw):
    for j in range(ncl):
        for k in range(np.shape(r)[2]):
            w[i,j,k]=P_Pg/(np.shape(r)[2])
for i in range(nsw):
    for j in range(ncl):
        hist,edge=np.histogram(r[i,j],bins=bins)
        R0_all[i,j]=(edge[:-1]+edge[1:])/2
        P0_all[i,j]=hist/np.sum(hist*dr)
        R0_mean_all[i,j]=np.mean(r[i,j])
        R0_std_all[i,j]=np.std(r[i,j])
        hist,edge=np.histogram(r_a[i,j],bins=bins)
        R0_a[i,j]=(edge[:-1]+edge[1:])/2
        P0_a[i,j]=hist/np.sum(hist*dr)
        R0_mean_a[i,j]=np.mean(r_a[i,j])
        R0_std_a[i,j]=np.std(r_a[i,j])
        hist,edge=np.histogram(r_a[i,j]/R0_mean_all[i,j],bins=nhist)
        R0_a_[i,j]=(edge[:-1]+edge[1:])/2
        P0_a_[i,j]=hist/np.sum(hist*(np.max(R0_a_[i,j])-np.min(R0_a_[i,j]))/(nhist-1))
        R0_mean_a_[i,j]=np.mean(r_a[i,j]/R0_mean_all[i,j])
        R0_std_a_[i,j]=np.std(r_a[i,j]/R0_mean_all[i,j])
        hist,edge=np.histogram(r_b[i,j],bins=bins)
        R0_b[i,j]=(edge[:-1]+edge[1:])/2
        P0_b[i,j]=hist/np.sum(hist*dr)
        R0_mean_b[i,j]=np.mean(r_b[i,j])
        R0_std_b[i,j]=np.std(r_b[i,j])
        hist,edge=np.histogram(r_b[i,j]/R0_mean_all[i,j],bins=nhist)
        R0_b_[i,j]=(edge[:-1]+edge[1:])/2
        P0_b_[i,j]=hist/np.sum(hist*(np.max(R0_b_[i,j])-np.min(R0_b_[i,j]))/(nhist-1))
        R0_mean_b_[i,j]=np.mean(r_b[i,j]/R0_mean_all[i,j])
        R0_std_b_[i,j]=np.std(r_b[i,j]/R0_mean_all[i,j])
        hist,edge=np.histogram(r[i,j],bins=bins,weights=w[i,j])
        R0_g[i,j]=(edge[:-1]+edge[1:])/2
        P0_g[i,j]=hist/np.sum(hist*dr)
        R0_mean_g[i,j]=np.sum(r[i,j]*w[i,j])
        R0_std_g[i,j]=(np.sum((r[i,j]-R0_mean_g[i,j])**2*w[i,j]))**0.5
        hist,edge=np.histogram(r[i,j]/R0_mean_all[i,j],bins=nhist,weights=w[i,j])
        R0_g_[i,j]=(edge[:-1]+edge[1:])/2
        P0_g_[i,j]=hist/np.sum(hist*(np.max(R0_g_[i,j])-np.min(R0_g_[i,j]))/(nhist-1))
        R0_mean_g_[i,j]=np.sum(r[i,j]*w[i,j]/R0_mean_all[i,j])
        R0_std_g_[i,j]=(np.sum((r[i,j]/R0_mean_all[i,j]-R0_mean_g_[i,j])**2*w[i,j]))**0.5

"""Read .pdb Files and Calculate Distance to Center"""
r=[]
for i in range(nsw):
    r+=[np.zeros((nfr,ntr[i],nat,3))]
for i in range(nsw):
    for j in range(nfr):
        r[i][j,:,:,:]=(pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')[1])

R=[]
for i in range(nsw):
    r0=np.mean(r[i],axis=-2,keepdims=True)
    R+=[np.linalg.norm(r[i]-r0,axis=-1)]
r=R

"""Collect and Analysis"""
idx_a=np.where(state==1)[0]
idx_b=np.where(state==0)[0]
r_a=[]
r_b=[]
for i in range(nsw):
    r_a+=[r[i][:,:,idx_a]]
    r_b+=[r[i][:,:,idx_b]]

dr=rmax/(nhist-1)
bins=np.arange(0-dr/2,rmax+dr,dr)
R_all=np.zeros((nsw,nfr,nhist))
P_all=np.zeros((nsw,nfr,nhist))
R_mean_all=np.zeros((nsw,nfr))
R_std_all=np.zeros((nsw,nfr))
R_a=np.zeros((nsw,nfr,nhist))
P_a=np.zeros((nsw,nfr,nhist))
R_mean_a=np.zeros((nsw,nfr))
R_std_a=np.zeros((nsw,nfr))
R_a_=np.zeros((nsw,nfr,nhist))
P_a_=np.zeros((nsw,nfr,nhist))
R_mean_a_=np.zeros((nsw,nfr))
R_std_a_=np.zeros((nsw,nfr))
R_b=np.zeros((nsw,nfr,nhist))
P_b=np.zeros((nsw,nfr,nhist))
R_mean_b=np.zeros((nsw,nfr))
R_std_b=np.zeros((nsw,nfr))
R_b_=np.zeros((nsw,nfr,nhist))
P_b_=np.zeros((nsw,nfr,nhist))
R_mean_b_=np.zeros((nsw,nfr))
R_std_b_=np.zeros((nsw,nfr))
R_g=np.zeros((nsw,nfr,nhist))
P_g=np.zeros((nsw,nfr,nhist))
R_mean_g=np.zeros((nsw,nfr))
R_std_g=np.zeros((nsw,nfr))
R_g_=np.zeros((nsw,nfr,nhist))
P_g_=np.zeros((nsw,nfr,nhist))
R_mean_g_=np.zeros((nsw,nfr))
R_std_g_=np.zeros((nsw,nfr))
w=[]
for i in range(nsw):
    w+=[np.zeros((ntr[i],nat))]
    for j in range(ntr[i]):
        w[i][j]=P_Pg/ntr[i]
for i in range(nsw):
    for j in range(nfr):
        hist,edge=np.histogram(r[i][j],bins=bins)
        R_all[i,j]=(edge[:-1]+edge[1:])/2
        P_all[i,j]=hist/np.sum(hist*dr)
        R_mean_all[i,j]=np.mean(r[i][j])
        R_std_all[i,j]=np.std(r[i][j])
        hist,edge=np.histogram(r_a[i][j],bins=bins)
        R_a[i,j]=(edge[:-1]+edge[1:])/2
        P_a[i,j]=hist/np.sum(hist*dr)
        R_mean_a[i,j]=np.mean(r_a[i][j])
        R_std_a[i,j]=np.std(r_a[i][j])
        hist,edge=np.histogram(r_a[i][j]/R_mean_all[i,j],bins=nhist)
        R_a_[i,j]=(edge[:-1]+edge[1:])/2
        P_a_[i,j]=hist/np.sum(hist*(np.max(R_a_[i,j])-np.min(R_a_[i,j]))/(nhist-1))
        R_mean_a_[i,j]=np.mean(r_a[i][j]/R_mean_all[i,j])
        R_std_a_[i,j]=np.std(r_a[i][j]/R_mean_all[i,j])
        hist,edge=np.histogram(r_b[i][j],bins=bins)
        R_b[i,j]=(edge[:-1]+edge[1:])/2
        P_b[i,j]=hist/np.sum(hist*dr)
        R_mean_b[i,j]=np.mean(r_b[i][j])
        R_std_b[i,j]=np.std(r_b[i][j])
        hist,edge=np.histogram(r_b[i][j]/R_mean_all[i,j],bins=nhist)
        R_b_[i,j]=(edge[:-1]+edge[1:])/2
        P_b_[i,j]=hist/np.sum(hist*(np.max(R_b_[i,j])-np.min(R_b_[i,j]))/(nhist-1))
        R_mean_b_[i,j]=np.mean(r_b[i][j]/R_mean_all[i,j])
        R_std_b_[i,j]=np.std(r_b[i][j]/R_mean_all[i,j])
        hist,edge=np.histogram(r[i][j],bins=bins,weights=w[i])
        R_g[i,j]=(edge[:-1]+edge[1:])/2
        P_g[i,j]=hist/np.sum(hist*dr)
        R_mean_g[i,j]=np.sum(r[i][j]*w[i])
        R_std_g[i,j]=(np.sum((r[i][j]-R_mean_g[i,j])**2*w[i]))**0.5
        hist,edge=np.histogram(r[i][j]/R_mean_all[i,j],bins=nhist,weights=w[i])
        R_g_[i,j]=(edge[:-1]+edge[1:])/2
        P_g_[i,j]=hist/np.sum(hist*(np.max(R_g_[i,j])-np.min(R_g_[i,j]))/(nhist-1))
        R_mean_g_[i,j]=np.sum(r[i][j]*w[i]/R_mean_all[i,j])
        R_std_g_[i,j]=(np.sum((r[i][j]/R_mean_all[i,j]-R_mean_g_[i,j])**2*w[i]))**0.5

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
for i in range(nsw):
    pyw.write(f'The Index of Cells, Location of All Loci, Its Probability Density, Location of Loci in Compartment A, '
              f'Its Probability Density, Relative Location of Loci in Compartment A, Its Probability Density, Similar '
              f'for Loci in Compartment B, and Weighted by Gene Density ({sw[i]}):\n')
    for j in range(ncl):
        for k in range(nhist):
            pyw.write(f'{j} {R0_all[i,j,k]} {P0_all[i,j,k]} '
                      f'{R0_a[i,j,k]} {P0_a[i,j,k]} {R0_a_[i,j,k]} {P0_a_[i,j,k]} '
                      f'{R0_b[i,j,k]} {P0_b[i,j,k]} {R0_b_[i,j,k]} {P0_b_[i,j,k]} '
                      f'{R0_g[i,j,k]} {P0_g[i,j,k]} {R0_g_[i,j,k]} {P0_g_[i,j,k]}\n')
    pyw.write('\n')
    pyw.write(f'The Index of Cells, Mean, and Standard Deviation of Location of All Loci, of Location of Loci in '
              f'Compartment A, of Relative Location of Loci in Compartment A, Similar for Loci in Compartment B, and '
              f'Weighted by Gene Density ({sw[i]}):\n')
    for j in range(ncl):
        pyw.write(f'{j} {R0_mean_all[i,j]} {R0_std_all[i,j]} '
                  f'{R0_mean_a[i,j]} {R0_std_a[i,j]} {R0_mean_a_[i,j]} {R0_std_a_[i,j]} '
                  f'{R0_mean_b[i,j]} {R0_std_b[i,j]} {R0_mean_b_[i,j]} {R0_std_b_[i,j]} '
                  f'{R0_mean_g[i,j]} {R0_std_g[i,j]} {R0_mean_g_[i,j]} {R0_std_g_[i,j]}\n')
    pyw.write('\n')
    pyw.write(f'Time, Location of All Loci, Its Probability Density, Location of Loci in Compartment A, Its Probability '
              f'Density, Relative Location of Loci in Compartment A, Its Probability Density, Similar for Loci in '
              f'Compartment B, and Weighted by Gene Density ({sw[i]}):\n')
    for j in range(nfr):
        for k in range(nhist):
            pyw.write(f'{t[j]} {R_all[i,j,k]} {P_all[i,j,k]} '
                      f'{R_a[i,j,k]} {P_a[i,j,k]} {R_a_[i,j,k]} {P_a_[i,j,k]} '
                      f'{R_b[i,j,k]} {P_b[i,j,k]} {R_b_[i,j,k]} {P_b_[i,j,k]} '
                      f'{R_g[i,j,k]} {P_g[i,j,k]} {R_g_[i,j,k]} {P_g_[i,j,k]}\n')
    pyw.write('\n')
    pyw.write(f'Time, Mean, and Standard Deviation of Location of All Loci, of Location of Loci in Compartment A, of '
              f'Relative Location of Loci in Compartment A, Similar for Loci in Compartment B, and Weighted by Gene '
              f'Density ({sw[i]}):\n')
    for j in range(nfr):
        pyw.write(f'{t[j]} {R_mean_all[i,j]} {R_std_all[i,j]} '
                  f'{R_mean_a[i,j]} {R_std_a[i,j]} {R_mean_a_[i,j]} {R_std_a_[i,j]} '
                  f'{R_mean_b[i,j]} {R_std_b[i,j]} {R_mean_b_[i,j]} {R_std_b_[i,j]} '
                  f'{R_mean_g[i,j]} {R_std_g[i,j]} {R_mean_g_[i,j]} {R_std_g_[i,j]}\n')
    pyw.write('\n')
pyw.close()










