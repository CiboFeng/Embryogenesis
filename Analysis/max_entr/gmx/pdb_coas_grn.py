"""Import Modules"""
import math as m
import numpy as np
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
from tot_len import tot_len

"""Set Arguments"""
ipdbfile='D:\Biophysics\MyPython\Inputs\\1ubq.pdb'
opdbfile='Polymer_1ubq_0.4.pdb'

"""Read"""
f=open(ipdbfile)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
x=[]
y=[]
z=[]
for i in range(len(ls_f)):
    if len(ls_f[i])>1:
        if ls_f[i][0]=='ATOM' and ls_f[i][2]=='CA':
            x+=[float(ls_f[i][6])]
            y+=[float(ls_f[i][7])]
            z+=[float(ls_f[i][8])]

"""Write to .pdb File"""
N=len(x)
pdb=open(opdbfile,'w')
pdb.write(f'REMARK Polymer Model, NumberofBeads   {N}, Separated by 4.00 (A)\n')
for j in range(N):
    serial=''
    for k in range(5-len(str(j+1))):
        serial+=' '
    serial+=str(j+1)
    sequence=''
    for k in range(4-len(str(j+1))):
        sequence+=' '
    sequence+=str(j+1)
    pdb.write('ATOM  '+serial+'  CA  GLY  '+sequence+'     '+f' {tot_len(x[j],5)}  {tot_len(y[j],5)}  {tot_len(z[j],5)}'+'  1.00  0.00\n')
pdb.write('END\n')
pdb.close()

"""Write to .gro Files"""
x=np.array(x)/10
y=np.array(y)/10
z=np.array(z)/10
xlo=np.min(x)
ylo=np.min(y)
zlo=np.min(z)
x-=xlo
y-=ylo
z-=zlo
xhi=np.max(x)
yhi=np.max(y)
zhi=np.max(z)
for i in range(1,101):
    dat=open(f'D:\Biophysics\MyPython\\applications\\chu_max_entr_prep\\grofiles_pdb_1ubq\\N{i}.gro','w')
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






