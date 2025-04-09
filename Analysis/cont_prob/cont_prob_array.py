"""Import Modules"""
import argparse
import numpy as np
import datetime
import multiprocessing
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pdbs import pdbs
from cont import cont

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-f",type=str,help="The .trj files")
args=parser.parse_args()
trjfile=args.f
cl=['zp','zm','esc']
sw=['zp_esc','zm_esc']
outname='cont_prob'
nat=600
nfr=56
ntr=[1214,974]
nsw=len(sw)
ncl=len(cl)

print(1,datetime.datetime.now())
"""Calculate Contact Probabilities for Cells"""
P0=np.zeros((ncl,nat,nat))
for i in range(ncl):
    _,r=pdbs(f'../trj/{cl[i]}/trj.pdb')
    P0[i]=cont(r/10,map='tanh',d0=0.6,beta=3).prob()
print(1,datetime.datetime.now())

"""Loop"""
P=np.zeros((nsw,nfr+2,nat,nat))
for i in range(nsw):
    """Define a Function for Reading .pdb Files Calculating Contact Probabilities"""
    def calc(j,q):
        P=np.zeros((nat,nat))
        _,r=pdbs(f'../trj/{sw[i]}/reshp_{j}/trj.pdb')
        P[:,:]=cont(r/10,map='tanh',d0=0.6,beta=3).prob()
        q.put((j,P))
        ### The j is for recording the order of queue.

    """Calculate Contact Probabilities in Parallel"""
    with multiprocessing.Manager() as manager:
        q=manager.Queue()
        procs=[]
        for j in range(nfr):
            p=multiprocessing.Process(target=calc,args=(j,q))
            p.start()
            procs.append(p)
        for p in procs:
            p.join()
        while not q.empty():
            j,p=q.get()
            P[i,j+1:j+2,:,:]=p
P[0,0]=P0[0]
P[0,-1]=P0[2]
P[1,0]=P0[1]
P[1,-1]=P0[2]

"""Output Data"""
np.save(f'cont_prob.npy',P)

