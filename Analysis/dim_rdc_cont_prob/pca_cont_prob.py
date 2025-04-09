"""Import Modules"""
import numpy as np
from sklearn.decomposition import PCA
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from dist_prob import dist_prob
from pca import pca

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
outname='pca_cont_prob'
cl=['ZP_','ZM_','ESC_','ZP','ZM','E2CP','E2CM','L2CP','L2CM','8CP','8CM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
idx_cl=[0,nfr+2,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1]

"""Read Data and Do PCA"""
P_ref=np.zeros((len(reffile),nat,nat))
for i in range(len(reffile)):
    f=open(reffile[i],'r')
    lsf=f.readlines()
    for j in range(len(lsf)):
        ls_f=lsf[j].strip('\n').split()
        P_ref[i,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
for i in range(len(reffile)):
    dp=dist_prob(P_ref[i])
    norm=np.zeros((nat,nat))
    for j in range(nat):
        for k in range(nat):
            norm[j,k]=dp[abs(j-k)]
    P_ref[i]-=norm
P_ref=P_ref.reshape(len(reffile),-1)

P=np.load(npyfile)
for i in range(nsw):
    for j in range(nfr+2):
        dp=dist_prob(P[i,j])
        norm=np.zeros((nat,nat))
        for k in range(nat):
            for l in range(nat):
                norm[k,l]=dp[abs(k-l)]
        P[i,j]-=norm
P=P.reshape(nsw*(nfr+2),-1)

P=np.concatenate((P,P_ref),axis=0)

pca=PCA(n_components=2)
pca.fit(P)
pc_nor=pca.transform(P)
perct=pca.explained_variance_ratio_

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
pyw.write('The Percentage of the First Component:\n')
pyw.write(f'{perct[0]}\n')
pyw.write('\n')
pyw.write('The Percentage of the Second Component:\n')
pyw.write(f'{perct[1]}\n')
pyw.write('\n')
for i in range(ncl):
    pyw.write(f'The First, Second Component of {cl[i]}:\n')
    pyw.write(f'{pc_nor[idx_cl[i],0]} {pc_nor[idx_cl[i],1]}\n')
    pyw.write('\n')
for i in range(nsw):
    pyw.write(f'The First, Second Component of {sw[i]}:\n')
    for j in range(i*(nfr+2)+1,(i+1)*(nfr+2)-1):
        pyw.write(f'{pc_nor[j,0]} {pc_nor[j,1]}\n')
    pyw.write('\n')
pyw.close()


