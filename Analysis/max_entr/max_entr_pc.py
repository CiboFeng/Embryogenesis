"""Import Modules"""
import argparse
import numpy as np
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from trj import trj
from cont import cont

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-s",type=str,help="The sysname name")
parser.add_argument("-w",type=str,help="The work directory")
parser.add_argument("-t",type=str,help="The .trj file")
parser.add_argument("-pc",type=str,help="The .npz file containing P and C")
args=parser.parse_args()
sysname=args.s
wkdir=args.w
trjfile=args.t
pcfile=args.pc
tpirfile=wkdir+f'/trn_pair_{sysname}.pyw'
d0=0.6
beta=3

"""Read Trained Pairs from the File"""
f=open(tpirfile)
lsf=f.readlines()
f.close()
ls_f=[]
for i in range(1,len(lsf)):
    l_f=lsf[i].strip('\n').split()
    ls_f.append(l_f)
tpir=[list(map(int,subls)) for subls in ls_f]

"""Read .trj File and Calculate"""
t,x,y,z=trj(f'{trjfile}')
p,C=cont(x,y,z,d0=d0,map='cos/sqr',tpir=tpir,beta=beta,wkdir=wkdir).intg()
P=cont(x,y,z,d0=d0,map='cos/sqr',tpir=None,beta=beta,wkdir=wkdir).prob()

"""Write to Files"""
np.savez(pcfile,p=p,C=C,P=P)
# def load_or_create():
#     try:
#         data=np.load(pcfile)
#         P0=data['P']
#         C0=data['C']
#         mark=data['Mark']
#         data.close()
#     except FileNotFoundError:
#         P0=np.zeros((1,nhipair))
#         C0=np.zeros((nhipair,nhipair))
#         mark=0
#         np.savez(pcfile,P=P0,C=C0,Mark=mark)
#     return P0,C0,mark
#
# def read_write():
#     lockfile=pcfile+".lock"
#     with open(lockfile,'w') as f:
#         fcntl.flock(f,fcntl.LOCK_EX)
#         P0,C0,mark=load_or_create()
#         P0+=P
#         C0+=C
#         mark+=1
#         np.savez(pcfile,P=P0,C=C0,Mark=mark)
#         fcntl.flock(f, fcntl.LOCK_UN)
#
# if __name__ == "__main__":
#     read_write()
#
# ######################################################
# def acquire_lock():
#     while True:
#         try:
#             os.makedirs(pcfile+'.lock')
#             return
#         except FileExistsError:
#             time.sleep(0.1)
#
# def release_lock():
#     os.rmdir(pcfile+'.lock')
#
# acquire_lock()
# try:
#     with np.load(pcfile) as data:
#         P0=data['P']
#         C0=data['C']
#         mark=data['Mark']
#         P0+=P
#         C0+=C
#         mark+=1
#     np.savez(pcfile,P=P0,C=C0,Mark=mark)
# finally:
#     release_lock()

