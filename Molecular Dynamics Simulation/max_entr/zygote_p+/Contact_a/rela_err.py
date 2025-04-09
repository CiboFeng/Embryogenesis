import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
file1='Probability_20.dat'
file2='../uij/PC_100000_iced_chr14_dense.matrix_P.dat'
N=600
P_sim=np.zeros((N,N))
f=open(file1)
lsf=f.readlines()
f.close()
for i in range(len(lsf)):
    ls_f=lsf[i].strip('\n').split()
    P_sim[int(float(ls_f[0]))-1,int(float(ls_f[1]))-1]=float(ls_f[2])
    P_sim[int(float(ls_f[1]))-1,int(float(ls_f[0]))-1]=float(ls_f[2])
P_exp=np.zeros((N,N))
f=open(file2)
lsf=f.readlines()
f.close()
for i in range(len(lsf)):
    ls_f=lsf[i].strip('\n').split()
    P_exp[int(float(ls_f[0]))-1,int(float(ls_f[1]))-1]=float(ls_f[2])
    P_exp[int(float(ls_f[1]))-1,int(float(ls_f[0]))-1]=float(ls_f[2])
dP=np.zeros((N,N))
for i in range(N):
    for j in range(N):
        if P_exp[i,j]!=0:
            dP[i,j]=abs(P_sim[i,j]-P_exp[i,j])/P_exp[i,j]
norm=colors.LogNorm()
plt.imshow(dP,norm=norm)
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()
