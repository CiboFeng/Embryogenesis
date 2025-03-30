"""Import Modules"""
import argparse
import numpy as np
import os
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from pdbs import pdbs
from sklearn.decomposition import PCA

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-i",type=str,help="The switching index")
parser.add_argument("-j",type=str,help="The frame index")
parser.add_argument("-k",type=str,help="The row index")
parser.add_argument("-l",type=str,help="The column index")
args=parser.parse_args()
i=int(args.i)
j=int(args.j)
k=int(args.k)
l=int(args.l)
outname='shp_seg'
cl=['ZP','ZM','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=['zp','zm','esc']
sw_=['zp_esc','zm_esc']
nsw=len(sw)
ncl=len(cl)
nfr=56
nat=600

"""Read .pdb Files and Calculate"""
if abs(k-l)>0:
    _,r=pdbs(f'../trj/{sw_[i]}/reshp_{j}/trj.pdb')
    if k>l:
        r=r[:,l:k,:]
    else:
        r=r[:,k:l,:]
    r-=np.mean(r,axis=1,keepdims=True)
    IT=np.sum(np.einsum('ijk,ijl->ijkl',r,r),axis=1)/nat
    tr=np.trace(IT,axis1=1,axis2=2)
    lmd=np.array([np.sort(np.linalg.eigvals(IT[m]))[::-1] for m in range(np.shape(IT)[0])])
    PA=np.mean(lmd**0.5,axis=0)
    Rg=np.mean(tr**0.5)
    Dlt=np.mean(9/2*np.var(lmd,axis=1)/tr**2)
    S=np.mean(27*np.product(lmd-np.mean(lmd,axis=1,keepdims=True),axis=1)/tr**3)

    if abs(k-l)==1:
        PA=np.zeros(3)
        Rg=0.0
        Dlt=0.0
        S=0.0
    if abs(k-l)==2:
        PA=np.abs(PA)
        PA[1]=0.0
        PA[2]=0.0
        Dlt=1.0
        S=2.0
    if abs(k-l)==3:
        PA=np.abs(PA)
        PA[2]=0.0
else:
    PA=np.zeros(3)
    Rg=0.0
    Dlt=0.0
    S=0.0

"""Output Data"""
if not os.path.exists(f'{outname}.pyw'):
    try:
        os.mkdir(f'{outname}.pyw')
    except IOError:
        pass
f=open(f'{outname}.pyw/_','w')
f.write('The Index of Switching Process, Frame, Row, Column, the Length of Principal Axis 1, 2, 3, Radius of Gyration, Aspherical Quantity, and S:\n')
f.close()
f=open(f'{outname}.pyw/{i}.0_{j}.0_{k}.0_{l}.0','w')
f.write(f'{PA[0]} {PA[1]} {PA[2]} {Rg} {Dlt} {S}\n')
f.close()

