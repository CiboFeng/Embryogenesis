"""Import Modules"""
import argparse
import numpy as np
from scipy.spatial.distance import squareform,pdist
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
cs={'zp':0.0,'zm':1.0,'esc':2.0}
c=cs[pdbfile.split('/')[-2].split('_')[-1]]
w=float(pdbfile.split('/')[-1].split('.pdb')[0].split('_')[1])
outname='rmsd'

""""""
t,r=pdbs(pdbfile)
nfr,nat,_=np.shape(r)
D0=squareform(pdist(r[0],'euclidean'))
d=np.zeros(nfr)
for i in range(nfr):
    d[i]=(np.mean((squareform(pdist(r[i],'euclidean'))-D0)**2))**0.5

"""Output Data"""
if not os.path.exists(f'{outname}.pyw'):
    try:
        os.mkdir(f'{outname}.pyw')
    except IOError:
        pass
f=open(f'{outname}.pyw/_','w')
f.write('Cell Type Index, Trajectory Index, Time, RMSD:\n')
f.close()
f=open(f'{outname}.pyw/{c}_{w}','w')
for i in range(nfr):
    f.write(f'{t[i]} {d[i]}\n')
f.close()