"""Import Modules"""
import numpy as np
from sklearn.decomposition import PCA
from scipy import interpolate
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw

"""Set Arguments"""
file='loci_distr.pyw'
outname='pca_loci_distr'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl=len(cl)
nfr=56
nat=600
nhist=50
idx=[2,4,6]

"""Read Data and Do PCA"""
r0=np.zeros((nsw,ncl,nhist,len(idx)))
P0=np.zeros((nsw,ncl,nhist,len(idx)))
for i in range(nsw):
    q=pyw(file,f'The Index of Cells, Location of All Loci, Its Probability Density, Location of Loci in Compartment A, '
               f'Its Probability Density, Relative Location of Loci in Compartment A, Its Probability Density, Similar '
               f'for Loci in Compartment B, and Weighted by Gene Density ({sw[i]}):')
    for j in range(ncl):
        for k in range(nhist):
            for l in range(len(idx)):
                r0[i,j,k,l]=q[2*idx[l]+1][j*nhist+k]
                P0[i,j,k,l]=q[2*idx[l]+2][j*nhist+k]
# rho0=P0/(4*np.pi*r0**2)
rho0=P0

t=np.zeros(nfr)
r=np.zeros((nsw,nfr,nhist,len(idx)))
P=np.zeros((nsw,nfr,nhist,len(idx)))
for i in range(nsw):
    q=pyw(file,f'Time, Location of All Loci, Its Probability Density, Location of Loci in Compartment A, Its Probability '
               f'Density, Relative Location of Loci in Compartment A, Its Probability Density, Similar for Loci in '
               f'Compartment B, and Weighted by Gene Density ({sw[i]}):')
    for j in range(nfr):
        t[j]=q[0][j*nhist]
        for k in range(nhist):
            for l in range(len(idx)):
                r[i,j,k,l]=q[2*idx[l]+1][j*nhist+k]
                P[i,j,k,l]=q[2*idx[l]+2][j*nhist+k]
# rho=P/(4*np.pi*r**2)
rho=P

r=np.concatenate((r0[:,0][:,np.newaxis],r),axis=1)
r=np.concatenate((r,r0[:,-1][:,np.newaxis]),axis=1)
rho=np.concatenate((rho0[:,0][:,np.newaxis],rho),axis=1)
rho=np.concatenate((rho,rho0[:,-1][:,np.newaxis]),axis=1)

r_=np.linspace(np.max(np.min(r,axis=2)),np.min(np.max(r,axis=2)),nhist)
rho_=np.zeros((nsw,nfr+2,nhist,len(idx)))
for i in range(nsw):
    for j in range(nfr+2):
        for k in range(len(idx)):
            f=interpolate.interp1d(r[i,j,:,k],rho[i,j,:,k],kind='cubic')
            rho_[i,j,:,k]=f(r_)

rho=rho_.reshape(nsw*(nfr+2),nhist,len(idx))
rho_nor=np.zeros((nsw*(nfr+2),2,len(idx)))
perct=np.zeros((2,len(idx)))
for i in range(len(idx)):
    pca=PCA(n_components=2)
    pca.fit(rho[:,:,i])
    rho_nor[:,:,i]=pca.transform(rho[:,:,i])
    perct[:,i]=pca.explained_variance_ratio_

import matplotlib.pyplot as plt
h=np.linalg.norm(rho_[:,1:-1,:,:]-(rho_[:,2:,:,:]+rho_[:,:-2,:,:])/2,axis=2)
w=np.linalg.norm(rho_[:,2:,:,:]-rho_[:,:-2,:,:],axis=2)
de=h/w
m=np.sum(rho_*r_[np.newaxis,np.newaxis,:,np.newaxis]*(np.max(r_)-np.min(r_)),axis=2)
h=np.abs(m[:,1:-1,:]-(m[:,2:,:]+m[:,:-2,:])/2)
w=np.abs(m[:,2:,:]-m[:,:-2,:])
dm=h/w
plt.plot(t[2:-9],de[0,1:-10,0],'r')
plt.plot(t[2:-9],de[0,1:-10,0],'r.',label='Euclidean')
plt.plot(t[2:-9],dm[0,1:-10,0],'b')
plt.plot(t[2:-9],dm[0,1:-10,0],'b.',label='Mean')
plt.xscale('log')
plt.legend()
plt.savefig('test.png')
plt.show()

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
pyw.write('The Percentage of the First, Second Component of Compartment A, B, and Gene:\n')
pyw.write(f'{perct[0,0]} {perct[1,0]} {perct[0,1]} {perct[1,1]} {perct[0,2]} {perct[1,2]}\n')
pyw.write('\n')
for i in range(nsw):
    pyw.write(f'The First, Second Component of Compartment A, B, and Gene ({sw[i]}, Reference):\n')
    for j in [i*(nfr+2),(i+1)*(nfr+2)-1]:
        pyw.write(f'{rho_nor[j,0,0]} {rho_nor[j,1,0]} {rho_nor[j,0,1]} {rho_nor[j,1,1]} {rho_nor[j,0,2]} {rho_nor[j,1,2]}\n')
    pyw.write('\n')
    pyw.write(f'The First, Second Component of Compartment A, B, and Gene ({sw[i]}):\n')
    for j in range(i*(nfr+2)+1,(i+1)*(nfr+2)-1):
        pyw.write(f'{rho_nor[j,0,0]} {rho_nor[j,1,0]} {rho_nor[j,0,1]} {rho_nor[j,1,1]} {rho_nor[j,0,2]} {rho_nor[j,1,2]}\n')
    pyw.write('\n')
pyw.close()


