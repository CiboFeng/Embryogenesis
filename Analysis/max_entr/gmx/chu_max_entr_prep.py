import sys
sys.path.append('D://Biophysics//MyPython//functions')
from pyw import pyw
import numpy as np
from tot_len import tot_len
import random

# f=open('Contact_a_0.dat','w')
# for i in range(76):
#     for j in range(i+2,76):
#         f.write(f'{i+1}\t{j+1}\t-9.99999978E-03\n')

# file='prob_76.pyw'
# p=pyw(file,'probability')[0]
# N=int(len(p)**0.5)
# P=np.zeros((N,N))
# for i in range(N):
#     for j in range(N):
#         P[i,j]=p[i*N+j]
# f=open('PC_100000_iced_chr14_dense.matrix_P.dat','w')
# for i in range(N):
#     for j in range(N):
#         f.write(f'{i+1}\t{j+1}\t{P[i,j]}\n')
# f=open('Probability.dat','w')
# for i in range(N):
#     for j in range(i+2,N):
#         f.write(f'{i+1}\t{j+1}\t{P[i,j]}\n')


pdbfile=f'Polymer_600_0.4.pdb'
N=600
R0=4
Rwall=R0/2*(10*N)**(1/3)-R0/2
dspace=2
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

"""Translation"""
x=np.array(x)+dspace*2*Rwall+Rwall
y=np.array(y)+dspace*2*Rwall+Rwall
z=np.array(z)+dspace*2*Rwall+Rwall

xmin=min(x)
xmax=max(x)
ymin=min(y)
ymax=max(y)
zmin=min(z)
zmax=max(z)
pdb=open(pdbfile,'w')
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
