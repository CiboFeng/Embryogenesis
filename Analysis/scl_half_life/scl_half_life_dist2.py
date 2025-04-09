"""Import Modules"""
import numpy as np
from scipy.spatial.distance import pdist,squareform
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs

"""Set Arguments"""
npzfile='scl_half_life_dist_.npz'
outname=npzfile[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
nfr0=400100
ntr=[1214,974]
dlt=0.2
scl=[[0,20],[20,60],[60,140],[140,300]]
nsc=len(scl)
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read Data"""
data=np.load(npzfile)
Q0=data['Q0']
Q=data['Q']

"""Normalize the Reaction Coordinate"""
Q0=np.zeros((nsw,nfr,nsc,nat))
Q_=np.zeros((nsw,nfr,nsc,nat))
Q1=np.zeros((nsw,nfr,nsc,nat))
Q0+=Q
Q_+=Q
Q1+=Q
Q0[:,:,:,:]=Q0[:,[1],:,:]
Q1[:,:,:,:]=Q1[:,[-10],:,:]
### Only assign in this way can it be ensured without confusion, otherwise the "Q0[:,:,:]=Q0[:,:,[0]]" .etc will affect the initial variables.
Q=(Q_-Q0)/(Q1-Q0)
# Q=(Q-Q0[:,0:1])/(Q0[:,1:2]-Q0[:,0:1])

"""Collect the Half Lives"""
hl=np.zeros((nsw,nsc,nat))
for i in range(nsw):
    for j in range(nsc):
       for k in range(nat):
          for l in range(nfr-1):
             if Q[i,l,j,k]<=0.5 and Q[i,l+1,j,k]>=0.5:
                hl[i,j,k]=((0.5-Q[i,l,j,k])*t[l]+(Q[i,l+1,j,k]-0.5)*t[l+1])/(Q[i,l+1,j,k]-Q[i,l,j,k])
                break

t1mean=np.mean(hl,axis=-1)
t2mean=np.mean(hl**2,axis=-1)

"""Output Data"""
np.save(f'{outname}.npy',Q)
pyw=open(f'{outname}.pyw','w')
for i in range(nsw):
    pyw.write(f'The First, Second Moment of Half Lives ({sw[i]}):\n')
    for j in range(nsc):
       pyw.write(f'{t1mean[i,j]} {t2mean[i,j]}\n')
    pyw.write('\n')
pyw.close()
