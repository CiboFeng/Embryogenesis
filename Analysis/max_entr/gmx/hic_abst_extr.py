import numpy as np

ifile='/hpc2hdd/home/chu-amat/cbfengphy/hic_fert/14th_sperm_pwk_rep123_100000_iced.matrix'
ofile='/hpc2hdd/home/chu-amat/cbfengphy/sperm/uij/PC_100000_iced_chr14_dense.matrix_P.dat'

f=open(ifile)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
ls_f=[list(map(float,sublist)) for sublist in ls_f]
N=int(ls_f[len(ls_f)-1][0])
P=np.zeros((N,N),dtype=float)
for i in range(len(ls_f)):
    #if len(ls_f[i])>0: #In order to make sure the list index in next line does not exceed the range.
    P[int(ls_f[i][0])-1,int(ls_f[i][1])-1]=ls_f[i][2]
    P[int(ls_f[i][1])-1,int(ls_f[i][0])-1]=ls_f[i][2]

idx_out=[]
for i in range(int(len(P)/2)):
    if P[i,i]!=0.0:
        idx_out+=[i]
        break
for i in range(len(P)-1,int(len(P)/2),-1):
    if P[i,i]!=0.0:
        idx_out+=[i]
        break
print(idx_out)

start=idx_out[0]
end=idx_out[1]+1
P=P[start:end,start:end]

f=open(ofile,'w')
N=len(P)
print(N)
for i in range(N):
    for j in range(N):
        f.write(f'{i+1} {j+1} {P[i,j]}\n')
f.close()
