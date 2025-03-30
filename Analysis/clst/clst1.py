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
from pdbs import pdbs

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-f",type=str,help="The .pdb files")
args=parser.parse_args()
pdbfile=args.f
cls={'zp':0.0,'zm':1.0,'esc':2.0}
cl=cls[pdbfile.split('/')[-2].split('_')[-1]]
outname1='samp'
outname2='dim_rdc'
outname3='clst1'

"""Read Data and Do Cluster"""
_,r=pdbs(pdbfile)
nfr,nat,_=np.shape(r)
D=np.zeros((nfr,int((nat-1)*(nat-2)/2)))
for i in range(nfr):
    d=squareform(pdist(r[i],'euclidean'))
    D[i]=d[np.triu_indices(nat,k=2)]
D_rdc=D[::200]
Z=linkage(D_rdc,method='ward')
Z[:,2]*=nat/((nat-1)*(nat-2)/2)

"""Dimensionality Reduction for Visualization"""
# DD=squareform(pdist(D,'euclidean'))
# mds=MDS(n_components=2,dissimilarity='precomputed',random_state=0)
# d=mds.fit_transform(DD)
### MDS is too slow.
pca=PCA(n_components=2)
d=pca.fit_transform(D)

"""Output Data"""
if not os.path.exists(f'{outname1}'):
    try:
        os.mkdir(f'{outname1}')
    except IOError:
        pass
np.save(f'{outname1}/{outname1}_{cl}.npy',D_rdc)

if not os.path.exists(f'{outname2}.pyw'):
    try:
        os.mkdir(f'{outname2}.pyw')
    except IOError:
        pass
f=open(f'{outname2}.pyw/_','w')
f.write('Cell Type Index, Reduced Coordinate 1, Reduced Coordinate 2:\n')
f.close()
f=open(f'{outname2}.pyw/{cl}','w')
for i in range(nfr):
    f.write(f'{d[i,0]} {d[i,1]}\n')
f.close()

if not os.path.exists(f'{outname3}.pyw'):
    try:
        os.mkdir(f'{outname3}.pyw')
    except IOError:
        pass
f=open(f'{outname3}.pyw/_','w')
f.write('Linkage Matrix:\n')
f.close()
f=open(f'{outname3}.pyw/{cl}','w')
for i in range(np.shape(Z)[0]):
    for j in range(np.shape(Z)[1]):
        f.write(f'{Z[i,j]}\n')
f.close()

