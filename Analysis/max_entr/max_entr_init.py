"""Import Modules"""
import argparse
import numpy as np
import os
import random
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-s",type=str,help="The sysname name")
parser.add_argument("-w",type=str,help="The work directory")
parser.add_argument("-a",type=str,help="The written alpha.pyw file")
parser.add_argument("-d",type=str,help="The directory of .dat files")
parser.add_argument("-n",type=str,help="The number of parallel trajectories")
args=parser.parse_args()
sysname=args.s
wkdir=args.w
oalffile=args.a
datdir=args.d
n=int(args.n)
pexpfile=wkdir+f'/exp_cont_prob_{sysname}.pyw'
alf0file=wkdir+f'/alf_{sysname}.pyw'
psimfile=wkdir+f'/sim_cont_prob_{sysname}.pyw'
mpirfile=wkdir+f'/mtn_pair_{sysname}.pyw'
tpirfile=wkdir+f'/trn_pair_{sysname}.pyw'
pmin1=0.1
pmin2=0.01
pmin3=0.001
rel_err=0.01
R0=0.4
dspace=10

"""Read Files, Determine the Maintained Pairs and Trained Pairs"""
p_exp_all=pyw(pexpfile,'Probability:')[0]
### Even though there is only one column of data, Pa is still a hierarchical list. So "[0]" is
### essential.
N=int((len(p_exp_all))**0.5)
P_exp=np.zeros((N,N))
for i in range(N):
    for j in range(N):
        P_exp[i,j]=p_exp_all[i*N+j]

if not os.path.exists(alf0file):
    tpir=[]
    for i in range(1,N,2):
        for j in range(1,i,2):
            if abs(i-j)>1 and P_exp[i,j]>=pmin1:
                tpir.append([i,j])
    for i in range(0,N,4):
        for j in range(0,i,4):
            if abs(i-j)>1 and P_exp[i,j]>=pmin2:
                tpir.append([i,j])
    ntpir=len(tpir)

else:
    p_sim_all=pyw(psimfile,'Probability:')[0]
    alf_0_all=pyw(alf0file,'Alpha')[0]
    P_sim=np.zeros((N,N))
    Alf_0=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            P_sim[i,j]=p_sim_all[i*N+j]
            Alf_0[i,j]=alf_0_all[i*N+j]
    f=open(mpirfile)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for i in range(1,len(lsf)):
        l_f=lsf[i].strip('\n').split()
        ls_f.append(l_f)
    mpir=[list(map(int,subls)) for subls in ls_f]
#####################################################################
    tpir=[]
    for i in range(0,N,2):
        for j in range(0,i,2):
            if [i,j] not in mpir and abs(i-j)>2 and P_exp[i,j]>=pmin3:
                if abs(P_sim[i,j]-P_exp[i,j])/P_exp[i,j]>rel_err: ### It is based on the previous "if" condition.
                    tpir.append([i,j])
    ntpir=len(tpir)
#####################################################################

pyw=open(tpirfile,'w')
pyw.write('The Trained Pairs (from 0):\n')
for i in range(ntpir):
    pyw.write(f'{tpir[i][0]} {tpir[i][1]}\n')

"""Write to Initial Alpha File"""
if not os.path.exists(alf0file):
    pyw=open(oalffile,'w')
    pyw.write('Alpha of Data-driven Term:\n')
    for i in range(N):
        for j in range(N):
            if [i,j] in tpir or [j,i] in tpir:
                pyw.write(f'{-0.01}\n')
            else:
                pyw.write(f'{0}\n')
else:
    pyw=open(oalffile,'w')
    pyw.write('Alpha of Data-driven Term:\n')
    for i in range(N):
        for j in range(N):
            if [i,j] in tpir or [j,i] in tpir:
                pyw.write(f'{-0.01}\n')
            elif [i,j] in mpir or [j,i] in mpir: ### The "if...if...else..." structure is wrong.
                pyw.write(f'{Alf_0[i,j]}\n')
            else:
                pyw.write(f'{0}\n')

"""Write to .dat Files Cyclically"""
for i in range(0,n+1):
    #"""Generate an Initial Conformation"""
    x=[0]
    y=[0]
    z=[0]
    for j in range(1,N):
        judge=True
        while judge:
            theta=random.uniform(0,np.pi)
            phi=random.uniform(0,2*np.pi)
            x1=x[j-1]+R0*np.sin(theta)*np.cos(phi)
            y1=y[j-1]+R0*np.sin(theta)*np.sin(phi)
            z1=z[j-1]+R0*np.cos(theta)
            dist=[]
            for k in range(len(x)):
                dist+=[((x1-x[k])**2+(y1-y[k])**2+(z1-z[k])**2)**0.5]
            dmin=min(dist)
            if (x1**2+y1**2+z1**2)**0.5<R0/2*(10*N)**(1/3)-R0/2 and dmin>0.07:
                judge=False
        x+=[x1]
        y+=[y1]
        z+=[z1]

    #"""Calculate the Spacial Range of the Box"""
    xlo=min(x)-dspace*(max(x)-min(x))
    xhi=max(x)+dspace*(max(x)-min(x))
    ylo=min(x)-dspace*(max(y)-min(y))
    yhi=max(x)+dspace*(max(y)-min(y))
    zlo=min(x)-dspace*(max(z)-min(z))
    zhi=max(x)+dspace*(max(z)-min(z))

    #"""Write to .dat File"""
    dat=open(datdir+f'/{sysname}_{i}.dat','w')
    dat.write('LAMMPS Data\n\n')
    dat.write(f'{N}\tatoms\n')
    dat.write(f'{N-1}\tbonds\n')
    dat.write(f'{N-2}\tangles\n\n')
    dat.write(f'{N}\tatom types\n')
    dat.write(f'1\tbond types\n')
    dat.write(f'1\tangle types\n\n')
    dat.write(f'{tot_len(xlo,16)}\t{tot_len(xhi,16)}\txlo\txhi\n')
    dat.write(f'{tot_len(ylo,16)}\t{tot_len(yhi,16)}\tylo\tyhi\n')
    dat.write(f'{tot_len(zlo,16)}\t{tot_len(zhi,16)}\tzlo\tzhi\n\n')
    dat.write('Atoms\n\n')
    for j in range(N):
        dat.write(f'{j+1}\t1\t{j+1}\t{tot_len(0,16)}\t{tot_len(x[j],16)}\t{tot_len(y[j],16)}\t{tot_len(z[j],16)}\n')
    dat.write('\n')
    dat.write('Bonds\n\n')
    for j in range(N-1):
        dat.write(f'{j+1}\t1\t{j+1}\t{j+2}\n')
    dat.write('\n')
    dat.write('Angles\n\n')
    for j in range(N-2):
        dat.write(f'{j+1}\t1\t{j+1}\t{j+2}\t{j+3}\n')
    dat.write('\n')
    dat.close()