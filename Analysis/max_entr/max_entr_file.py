"""Import Modules"""
import argparse
import os
import numpy as np
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
parser.add_argument("-a",type=str,help="The inputting alpha.pyw file")
parser.add_argument("-c",type=str,help="The outputting .cfg file")
parser.add_argument("-p",type=str,help="The outputting .par file")
parser.add_argument("-thermo",type=int,help="Print the thermodynamic information every #thermo timesteps")
parser.add_argument("-Th",type=float,help="The high temperature of annealing algorithm")
parser.add_argument("-Tl",type=float,help="The low temperature of annealing algorithm")
parser.add_argument("-na",type=int,help="The number of steps of simulated annealing")
parser.add_argument("-nr",type=int,help="The number of steps of relaxation")
parser.add_argument("-dump",type=int,help="Dump trajectory every #dump timesteps")
parser.add_argument("-ne",type=int,help="The number of steps of equilibrium simulation")
args=parser.parse_args()
sysname=args.s
wkdir=args.w
alffile=args.a
cfgfile=args.c
parfile=args.p
nthermo=int(args.thermo)
Thi=float(args.Th)
Tlo=float(args.Tl)
nstep_annl=int(args.na)
nstep_relx=int(args.nr)
ndump=int(args.dump)
nstep_equl=int(args.ne)
mpirfile=wkdir+f'/mtn_pair_{sysname}.pyw'
tpirfile=wkdir+f'/trn_pair_{sysname}.pyw'
logfile=f'{sysname}.log'
idatfile=f'{sysname}i.dat'
trjfile=f'{sysname}.trj'
odatfile=f'{sysname}o.dat'
d0=0.6
beta=3
R0=0.4
mass=1
Kb=450
Ka=2

"""Read Maintained Pairs, Trained Pairs from the File"""
if os.path.exists(mpirfile):
    f=open(mpirfile)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for i in range(1,len(lsf)):
        l_f=lsf[i].strip('\n').split()
        ls_f.append(l_f)
    mpir=[list(map(int,subls)) for subls in ls_f]
else:
    mpir=[]
nmpir=len(mpir)

f=open(tpirfile)
lsf=f.readlines()
f.close()
ls_f=[]
for i in range(1,len(lsf)):
    l_f=lsf[i].strip('\n').split()
    ls_f.append(l_f)
tpir=[list(map(int,subls)) for subls in ls_f]

"""Read the Alpha File and Store alphas to Matrix"""
alf=pyw(alffile,'Alpha')[0]
N=int((len(alf))**0.5)
Alf=np.zeros((N,N),dtype=float)
for i in range(N):
    for j in range(N):
        Alf[i,j]=alf[i*N+j]

"""Write to .cfg file"""
cfg=open(cfgfile,'w')
cfg.write(f'log\t\t{logfile}\n')
cfg.write('units\t\tlj\n')
cfg.write('boundary\tm m m\n')
cfg.write('timestep\t0.01\n')
cfg.write('neighbor\t2.0 bin\n')
cfg.write('neigh_modify\tdelay 1 every 1\n')
cfg.write('comm_modify\tcutoff 7.0\n\n')
cfg.write('atom_style\tfull\n')
cfg.write('bond_style\tharmonic\n')
cfg.write('angle_style\tcosine\n')
cfg.write(f'pair_style\thybrid/overlay cosine/squared {d0+0.5*np.pi/beta} cosine/squared 0.44\n')
cfg.write('special_bonds\tlj 0.0 0.0 1.0\n\n')
cfg.write(f'read_data\t{idatfile}\n')
cfg.write(f'include\t\t{parfile}\n\n')
cfg.write(f'thermo\t\t{nthermo}\n')
cfg.write('thermo_style\tcustom step temp ke pe ebond eangle epair\n\n')
cfg.write(f'region\t\tnucleus sphere 0.0 0.0 0.0 {R0/2*(10*N)**(1/3)-R0/2+1.0} side in units box\n')
cfg.write('fix\t\twall all wall/region nucleus lj126 1.0 1.0 1.0\n\n')
cfg.write('minimize\t1.0e-6\t1.0e-8\t1000\t100000\n')
cfg.write(f'velocity\tall create {Thi} 5725 rot yes dist gaussian\n')
cfg.write('fix\t\tbalance all balance 1000 1.1 shift xyz 10 1.01\n\n')
cfg.write('fix\t\tnve all nve\n\n')
cfg.write(f'fix\t\tlangevin1 all langevin {Thi} {Tlo} 2.0 5725 zero yes\n')
cfg.write(f'run\t\t{int(nstep_annl)}\n')
cfg.write('unfix\t\tlangevin1\n\n')
cfg.write(f'fix\t\tlangevin2 all langevin {Tlo} {Tlo} 2.0 5725 zero yes\n')
cfg.write(f'run\t\t{nstep_relx}\n')
cfg.write('unfix\t\tlangevin2\n\n')
cfg.write(f'fix\t\tlangevin3 all langevin {Tlo} {Tlo} 2.0 5725 zero yes\n')
cfg.write(f'dump\t\tdump all custom {ndump} {trjfile} id type x y z\n')
cfg.write('dump_modify\tdump sort id\n')
cfg.write(f'run\t\t{nstep_equl}\n\n')
cfg.write(f'write_data\t{odatfile}')
cfg.close()

"""Write to .par File"""
par=open(parfile,'w')
par.write('#Mass\n\n') ### An empty line is required.
for i in range(N):
    par.write(f'mass\t{i+1}\t{tot_len(mass,16)}\n')
par.write('\n')
par.write('#Bond Coeff\n\n')
par.write(f'bond_coeff\t1\t{tot_len(Kb,16)}\t{tot_len(R0,16)}\n\n')
par.write('#Angle Coeff\n\n')
par.write(f'angle_coeff\t1\t{tot_len(Ka,16)}\n\n')
par.write('#Pair Coeff\n\n')
for i in range(N):
    for j in range(i+1):
        if [i,j] in mpir:
            par.write(f'pair_coeff\t{j+1}\t{i+1}\tcosine/squared\t1\t{tot_len(-Alf[i,j],16)}\t{tot_len(d0-0.5*np.pi/beta,16)}\t{tot_len(d0+0.5*np.pi/beta,16)}\n') ### I J..., I<J
        if [i,j] in tpir:
            par.write(f'pair_coeff\t{j+1}\t{i+1}\tcosine/squared\t1\t{tot_len(-Alf[i,j],16)}\t{tot_len(d0-0.5*np.pi/beta,16)}\t{tot_len(d0+0.5*np.pi/beta,16)}\n')
par.write(f'pair_coeff\t*\t*\tcosine/squared\t2\t{tot_len(-4,16)}\t{tot_len(0.36,16)}\t{tot_len(0.44,16)}\n')
### If the pair style is used multiple times in the pair_style command, then an additional numeric
### argument must also be specified which is a number from 1 to M where M is the number of times
### the sub-style was listed in the pair style command. The extra number indicates which instance
### of the sub-style these coefficients apply to.
par.close()

"""Write to .tab File"""
# tab=open(tabfile,'w')
# tab.write('# DATE: 2023-09-06 UNITS: real CONTRIBUTOR: Cibo_Feng\n\n')
# tab.write(f'SOFT_CORE_NOTANH\n')
# tab.write(f'N\t{n}\n\n')
# r,p,f=table_soft_core(n=n,rmin=rmin,rmax=rmax,r0=r0,epsilon=epsilon,sigma=sigma,Ecut=Ecut,alpha=0,rc=rc)
# for i in range(n):
#     tab.write(f'{i+1}\t{tot_len(r[i],16)}\t{tot_len(p[i],16)}\t{tot_len(f[i],16)}\n')
# tab.write('\n')
# for i in range(N):
#     for j in range(i-1):
#         if [i,j] in hipair:
#             tab.write(f'SOFT_CORE_{i+1}_{j+1}\n')
#             tab.write(f'N\t{n}\n\n')
#             r,U,F=table_soft_core(n=n,rmin=rmin,rmax=rmax,r0=r0,epsilon=epsilon,sigma=sigma,Ecut=Ecut,alpha=alpha[i,j],rc=rc)
#             for k in range(n):
#                 tab.write(f'{k+1}\t{tot_len(r[k],16)}\t{tot_len(U[k],16)}\t{tot_len(F[k],16)}\n')
#             tab.write('\n')

"""View the Force Fields"""
# rmin=sigma*((1+(600/epsilon)**0.5)/2)**(-1/6)
# rmax=R0*(1-np.exp(-2*600/Kb/R0**2))**0.5
# r=np.linspace(rmin,rmax,1000)
# Ub1=-0.5*Kb*R0**2*np.log(1-(r/R0)**2)
# Ub2=[]
# for i in range(len(r)):
#     if r[i]<2**(1/6)*sigma:
#         Ub2+=[4*epsilon*((sigma/r[i])**12-(sigma/r[i])**6)+epsilon]
#     else:
#         Ub2+=[0]
# Ub2=np.array(Ub2)
# Ub=Ub1+Ub2
#
# plt.figure(figsize=(14,7))
# wth=3
# size=20
# lenmaj=15
# lenmin=8
#
# ax = plt.subplot(221)
# plt.plot(r,Ub,'-',color='b',linewidth=wth)
# ax.minorticks_on()
# ax.tick_params(axis="both", which="major", direction="in", width=wth, length=lenmaj, labelsize=size)
# ax.tick_params(axis="both", which="minor", direction="in", width=wth, length=lenmin, labelsize=size)
# plt.xlim(-0.08,1.58)
# plt.ylim(-50,650)
# axis = plt.gca()
# axis.xaxis.set_major_locator(MultipleLocator(0.25))
# axis.xaxis.set_minor_locator(MultipleLocator(0.125))
# axis.yaxis.set_major_locator(MultipleLocator(150))
# axis.yaxis.set_minor_locator(MultipleLocator(75))
# axis.spines['bottom'].set_linewidth(wth)
# axis.spines['top'].set_linewidth(wth)
# axis.spines['left'].set_linewidth(wth)
# axis.spines['right'].set_linewidth(wth)
# plt.xlabel('$r$',fontsize=size)
# plt.ylabel('$U_\mathrm{bond}$',fontsize=size)
#
# theta=np.linspace(0.001*theta0,0.999*theta0,100)
# Ua=Ka*(1-np.cos(theta-theta0))
#
# ax = plt.subplot(222)
# plt.plot(theta,Ua,'-',color='b',linewidth=wth)
# ax.minorticks_on()
# ax.tick_params(axis="both",which="major",direction="in",width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis="both",which="minor",direction="in",width=wth,length=lenmin,labelsize=size)
# # plt.xlim(0,3.5)
# # plt.ylim(0,0.4)
# axis = plt.gca()
# axis.xaxis.set_major_locator(MultipleLocator(0.5))
# axis.xaxis.set_minor_locator(MultipleLocator(0.25))
# axis.yaxis.set_major_locator(MultipleLocator(1))
# axis.yaxis.set_minor_locator(MultipleLocator(0.5))
# axis.spines['bottom'].set_linewidth(wth)
# axis.spines['top'].set_linewidth(wth)
# axis.spines['left'].set_linewidth(wth)
# axis.spines['right'].set_linewidth(wth)
# plt.xlabel('$\\theta$',fontsize=size)
# plt.ylabel('$U_\mathrm{angle}$',fontsize=size)
#
# r=np.linspace(0.001*rc,0.999*rc,100)
# Up=1+np.cos(np.pi*r/rc)
#
# ax = plt.subplot(223)
# plt.plot(r,Up,'-',color='b',linewidth=wth)
# ax.minorticks_on()
# ax.tick_params(axis="both",which="major",direction="in",width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis="both",which="minor",direction="in",width=wth,length=lenmin,labelsize=size)
# # plt.xlim(0,5)
# # plt.ylim(0,2)
# axis = plt.gca()
# axis.xaxis.set_major_locator(MultipleLocator(1))
# axis.xaxis.set_minor_locator(MultipleLocator(0.5))
# axis.yaxis.set_major_locator(MultipleLocator(0.5))
# axis.yaxis.set_minor_locator(MultipleLocator(0.25))
# axis.spines['bottom'].set_linewidth(wth)
# axis.spines['top'].set_linewidth(wth)
# axis.spines['left'].set_linewidth(wth)
# axis.spines['right'].set_linewidth(wth)
# plt.xlabel('$r$',fontsize=size)
# plt.ylabel('$U_\mathrm{pair}$',fontsize=size)
#
# plt.tight_layout()
# plt.savefig(f'force_field_{sysname}.png',format='png')
# plt.savefig(f'force_field_{sysname}.pdf',format='pdf')
# plt.show()



