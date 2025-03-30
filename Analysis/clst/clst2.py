"""Import Modules"""
import argparse
import numpy as np
from scipy.spatial.distance import squareform,pdist
from scipy.cluster.hierarchy import linkage,fcluster
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

"""Set Arguments"""
filep1='samp'
file2='clst1.pyw'
cls={'ZP':0.0,'ZM':1.0,'ESC':2.0}
ncl=len(cls)
outname='clst2'
sgm=4
cut=np.array([6,6,4])*sgm
nat=600
nhist=20
beta=3 ### In nm^-1.
d0=0.6 ### In nm.

"""Read Data and Do Cluster"""
D_rdc=[]
for i in range(ncl):
    D_rdc+=[np.load(f'{filep1}/{filep1}_{cls[list(cls.keys())[i]]}.npy')]
D_rdc=np.array(D_rdc)

q=pyw(file2,dir=True,mat=False)
Z=[]
for i in range(ncl):
    z=np.zeros((int(len(q[i])/4),4))
    for j in range(int(len(q[i])/4)):
        for k in range(4):
            z[j,k]=q[i][j*4+k][0]
    Z+=[z]

cent_idx=[]
for i in range(ncl):
    clst=fcluster(Z[i],cut[i],criterion='distance')
    cent_idx.append([])
    for j in np.unique(clst):
        p=np.where(clst==j)[0]
        o=np.mean(D_rdc[i,p],axis=0)
        d=np.linalg.norm(D_rdc[i,p]-o,axis=1)
        idx=p[np.argmin(d)]
        cent_idx[i]+=[idx]

"""Check the Distance Matrix of the Selected Samplings"""
DD=np.zeros((ncl,nhist))
PP=np.zeros((ncl,nhist))
for i in range(ncl):
    dd=squareform(pdist(D_rdc[i,cent_idx[i]],'euclidean'))[np.triu_indices(len(cent_idx[i]),k=1)]*nat/np.shape(D_rdc)[2]
    hist,edge=np.histogram(dd,bins=nhist)
    dd=(edge[:-1]+edge[1:])/2
    DD[i]=dd
    PP[i]=hist

"""Check the Reproducibility of Contact Probability Matrix"""
P=np.zeros((ncl,int((nat-1)*(nat-2)/2)))
for i in range(ncl):
    P[i]=np.mean(0.5*(1-np.tanh(beta*(D_rdc[i,cent_idx[i]]/10-d0))),axis=0)

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
for i in range(ncl):
    pyw.write(f'Cluster Center Index of {list(cls.keys())[i]}:\n')
    for j in range(len(cent_idx[i])):
        pyw.write(f'{200*cent_idx[i][j]}\n')
    pyw.write('\n')
for i in range(ncl):
    pyw.write(f'Distance among Selected Conformations of {list(cls.keys())[i]}, Population:\n')
    for j in range(nhist):
        pyw.write(f'{DD[i,j]} {PP[i,j]}\n')
    pyw.write('\n')
for i in range(ncl):
    pyw.write(f'Contact Probability Matrix of Selected Conformations of {list(cls.keys())[i]}:\n')
    for j in range(int((nat-1)*(nat-2)/2)):
        pyw.write(f'{P[i,j]}\n')
    pyw.write('\n')
pyw.close()
