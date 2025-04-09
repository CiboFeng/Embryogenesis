"""Import Modules"""
import numpy as np
import argparse
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('D:\Biophysics\MyPython\\functions')
from tot_len import tot_len
from pyw import pyw
import random
# from table_soft_core import table_soft_core
# import matplotlib.pyplot as plt
# from matplotlib.pyplot import MultipleLocator

"""Set Arguments"""
# parser=argparse.ArgumentParser()
# parser.add_argument("-w",type=str,help="The work directory")
# parser.add_argument("-a",type=str,help="The inputting alpha.pyw file")
# parser.add_argument("-c",type=str,help="The outputting .cfg file")
# parser.add_argument("-p",type=str,help="The outputting .par file")
# parser.add_argument("-t",type=str,help="The outputting .tab file")
# parser.add_argument("-thermo",type=int,help="Print the thermodynamic information every #thermo timesteps")
# parser.add_argument("-Th",type=float,help="The high temperature of annealing algorithm")
# parser.add_argument("-Tl",type=float,help="The low temperature of annealing algorithm")
# parser.add_argument("-na",type=int,help="The number of steps of simulated annealing")
# parser.add_argument("-nr",type=int,help="The number of steps of relaxation")
# parser.add_argument("-dump",type=int,help="Dump trajectory every #dump timesteps")
# parser.add_argument("-ne",type=int,help="The number of steps of equilibrium simulation")
# args=parser.parse_args()
# wkdir=args.w
# alffile=args.a
# cfgfile=args.c
# parfile=args.p
# tabfile=args.t
# nthermo=int(args.thermo)
# Thi=float(args.Th)
# Tlo=float(args.Tl)
# nstep_annl=int(args.na)
# nstep_relx=int(args.nr)
# ndump=int(args.dump)
# nstep_equl=int(args.ne)
sysname='1ubq'
mapfile='D:\Biophysics\MyPython\inputs\\b_cell_cont_map.pyw'
# mapfile=wkdir+'/tse_cont_map_1ubq_143.pyw'
alffile='alf_1ubq.pyw'
cfgfile='1ubq.cfg'
logfile='1ubq.log'
idatfile='1ubqi.dat'
parfile='1ubq.par'
trjfile='1ubq.trj'
odatfile='1ubqo.dat'
tabfile='1ubq.tab'
nthermo=100000
Thi=1008
Tlo=504
nstep_annl=2000000
nstep_relx=1000000
ndump=100000
nstep_equl=1000000000
pmin=0
R0=0.4
dspace=2
mass=1
epsilon=200
sigma=200
Ecut=4
r0=(2/(1+2**0.5))**(1/6)
rc=2
Kb=450
Ka=2
theta0=180
n=500
rmin=0.000001
rmax=5

"""Find the high-frequency contacted pairs"""
N=600
# Pa=pyw(mapfile,'Map')[0]
# ### Even though there is only one column of data, Pa is still a hierarchical list. So "[0]" is
# ### essential.
# N=int((len(Pa))**0.5)
# P=np.zeros((N,N),dtype=float)
# for i in range(N):
#     for j in range(N):
#         P[i,j]=Pa[i*N+j]
# hipair=[]
# for i in range(N):
#     for j in range(i):
#         if P[i,j]>=pmin and abs(i-j)>0:
#             hipair.append([i,j])

"""Write to Initial Alpha File"""
# alpha=open(alffile,'w') # If assigned to "pyw", the following pyw function will can't be used.
# alpha.write('Alpha of data-driven term:\n')
# for i in range(N**2):
#     if i%1000==0:
#         print(f'{i} of {N**2}: {i/N**2*100}%')
#     j,k=divmod(i,N)
#     if [j,k] in hipair or [k,j] in hipair:
#         alpha.write(f'{-0.01}\n')
#     else:
#         alpha.write(f'{0}\n')

"""Write to .gro Files Cyclically"""
for i in range(0,201):
    #"""Generate an Initial Conformation"""
    Rwall=R0/2*(10*N)**(1/3)-R0/2
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
            if (x1**2+y1**2+z1**2)**0.5<Rwall and dmin>0.07:
                judge = False
        x+=[x1]
        y+=[y1]
        z+=[z1]

    #"""Translation"""
    x=np.array(x)+dspace*2*Rwall+Rwall
    y=np.array(y)+dspace*2*Rwall+Rwall
    z=np.array(z)+dspace*2*Rwall+Rwall
    xhi=(dspace*2*Rwall+Rwall)*2
    yhi=(dspace*2*Rwall+Rwall)*2
    zhi=(dspace*2*Rwall+Rwall)*2

    #"""Write to .gro File"""
    dat=open(f'D:\Biophysics\MyPython\\applications\\land_swch\\max_entr\\gmx\\grofiles_600\\N{i}.gro','w')
    dat.write('Gro\n')
    dat.write(f'{N}\n')
    for j in range(N):
        if j+1<10:
            dat.write(f'    {j+1}GLY     CA    {j+1}   {tot_len(x[j],4)}   {tot_len(y[j],4)}   {tot_len(z[j],4)}\n')
        if j+1>=10 and j+1<100:
            dat.write(f'   {j+1}GLY     CA   {j+1}   {tot_len(x[j],4)}   {tot_len(y[j],4)}   {tot_len(z[j],4)}\n')
        if j+1>=100 and j+1<1000:
            dat.write(f'  {j+1}GLY     CA  {j+1}   {tot_len(x[j],4)}   {tot_len(y[j],4)}   {tot_len(z[j],4)}\n')
        if j+1>=1000 and j+1<10000:
            dat.write(f' {j+1}GLY     CA {j+1}   {tot_len(x[j],4)}   {tot_len(y[j],4)}   {tot_len(z[j],4)}\n')

    dat.write(f'{xhi} {yhi} {zhi}')

"""Read the Alpha File and Store alphas to Matrix"""
# alf=pyw(alffile,'Alpha')
# alpha=np.zeros((N,N),dtype=float)
# for i in range(N):
#     for j in range(N):
#         alpha[i,j]=alf[0][i*N+j]

"""Write to .top File"""
# par=open(parfile,'w')
#
# par.write('[ defaults ]\n')
# par.write('1  1  no\n\n')
# par.write('[ atomtypes ]\n')
# par.write('CA 1.000 0.000 A 0.00000E+00 0.10000E+01\n\n')
# par.write('[ moleculetype ]\n')
# par.write('Protein 1\n\n')
#
# par.write('[ atoms ]\n')
# for i in range(N):
#     par.write(f'{i+1}\tCA\t{i+1}\tGLY\tCA\t{i+1}\t{tot_len(0,4)}\t{tot_len(mass,4)}\n')
# par.write('\n')
#
# par.write('[ pairs ]\n')
# for i in range(N):
#     for j in range(i):
#         par.write(f'{j+1}\t{i+1}\t1\t{tot_len(0,6)}\t{tot_len(alpha[i,j],6)}\n') # pair_coeff i j ...: i<j
# par.write('\n')
#
# par.write('[ bonds ]\n')
# for i in range(N-1):
#     par.write(f'{i+1}\t{i+2}\t8\t0\t{tot_len(Kb,6)}\n')
#
# par.write('[ angles ]\n')
# for i in range(N-2):
#     par.write(f'{i+1}\t{i+2}\t{i+3}\t8\t0\t{tot_len(Ka,6)}\n')
#
# par.write('[ system ]\n')
# par.write('Protein\n\n')
# par.write('[ molecules ]\n')
# par.write('Protein 1')


