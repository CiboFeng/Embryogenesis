"""Import Modules"""
import numpy as np
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from cpmt_sgnl import cpmt_sgnl
from split_data import split_data
from insul_score import insul_score

"""Set Arguments"""
csvfile='../loci_distr/mm9.csv'
reffile='../cont_prob/ref/esc.dat'
tclfile1='gene.tcl'
tclfile2='cpmt.tcl'
tclfile3='tad.tcl'
nat=600
R=1.5
c1=9
c2=10

"""Read and Calculate the Gene Density"""
f=open(csvfile)
lsf=f.readlines()
Pg=[]
for i in range(len(lsf)):
    ls_f=lsf[i].strip('\n').split(',')
    if ls_f[1]=='chr15' or ls_f[1]=='"chr15"':
        pg=(float(ls_f[2])+float(ls_f[3]))/2/1e5
        if 250<pg<=850:
            Pg+=[pg-251]

bins=np.arange(0-1/2,nat,1)
hist,edge=np.histogram(Pg,bins=bins)
Pg=(edge[:-1]+edge[1:])/2
P_Pg=hist/np.sum(hist)

"""Color According to Genes"""
tcl=open(tclfile1,'w')
for i in range(nat):
    if P_Pg[i]>0:
        tcl.write('mol addrep 0\n')
        tcl.write(f'mol modstyle {i+1} 0 Tube {R} 12.000000\n')
        tcl.write(f'mol modselect {i+1} 0 index {i}\n')
        tcl.write(f'mol modcolor {i+1} 0 ColorID {c1}\n')
    else:
        tcl.write('mol addrep 0\n')
        tcl.write(f'mol modstyle {i+1} 0 Tube {R} 12.000000\n')
        tcl.write(f'mol modselect {i+1} 0 index {i}\n')
        tcl.write(f'mol modcolor {i+1} 0 ColorID {c2}\n')
tcl.write('mol showrep 0 0 0')
tcl.close()

# """Read and Calculate the Compartments"""
# P0=np.zeros((nat,nat))
# f=open(reffile,'r')
# lsf=f.readlines()
# for k in range(len(lsf)):
#     ls_f=lsf[k].strip('\n').split()
#     P0[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
# CS=cpmt_sgnl(P0)[0]
# if np.corrcoef(CS,P_Pg)[0,1]<0:
#     CS*=-1
#
# """Color According to Compartments"""
# x=np.arange(0,nat,1)
# x,CS=split_data(x,CS,0)
#
# tcl=open(tclfile2,'w')
# for i in range(len(CS)):
#     if np.mean(CS[i])>0:
#         tcl.write('mol addrep 0\n')
#         tcl.write(f'mol modstyle {i+1} 0 Tube {R} 12.000000\n')
#         tcl.write(f'mol modselect {i+1} 0 index {int(x[i][0])} to {int(x[i][-1])+1}\n')
#         tcl.write(f'mol modcolor {i+1} 0 ColorID {c1}\n')
#     else:
#         tcl.write('mol addrep 0\n')
#         tcl.write(f'mol modstyle {i+1} 0 Tube {R} 12.000000\n')
#         tcl.write(f'mol modselect {i+1} 0 index {int(x[i][0])} to {int(x[i][-1])+1}\n')
#         tcl.write(f'mol modcolor {i+1} 0 ColorID {c2}\n')
# tcl.write('mol showrep 0 0 0')
# tcl.close()
#
# """TAD"""
# pos=insul_score(P0,bin_size=10**5)[1]
# tcl=open(tclfile3,'w')
# for i in range(len(pos)-1):
#     tcl.write('mol addrep 0\n')
#     tcl.write(f'mol modstyle {i+1} 0 Tube {R} 12.000000\n')
#     tcl.write(f'mol modselect {i+1} 0 index {pos[i]} to {pos[i+1]}\n')
#     tcl.write(f'mol modcolor {i+1} 0 ColorID {c1}\n')
# tcl.write('mol showrep 0 0 0')
# tcl.close()
