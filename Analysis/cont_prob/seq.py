"""Import Modules"""
import numpy as np
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from dist_prob import dist_prob
from cpmt_sgnl import cpmt_sgnl
from insul_score import insul_score

"""Set Arguments"""
npyfile='cont_prob.npy'
outname='seq'
cl=['ZP','ZM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
t0=[0,0.01,
    0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
    0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
    2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
    20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
    200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
    2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
t0=[0,0.01,0.02,0.03,0.04]
t=[0,0.01,0.1,1.0,10.0,100.0,1000.0,10000.0]
t=[0,0.01,0.02,0.03,0.04]
idx=[0]
for i in range(len(t)):
    idx+=[t0.index(t[i])+1]
idx+=[-1]
idx_cl=[(0,0),(1,0),(0,-1)]

"""Read Data and Do Analysis"""
P_all=np.load(npyfile)
P=np.zeros((nsw,len(idx),nat,nat))
P_nor1=np.zeros((nsw,len(idx),nat,nat))
P_nor2=np.zeros((nsw,len(idx),nat,nat))
DP=np.zeros((nsw,len(idx),nat))
CS=np.zeros((nsw,len(idx),nat))
IS=np.zeros((nsw,len(idx),nat))
for i in range(nsw):
    P_temp=P_all[i]
    for j in range(len(idx)):
        P[i,j]=P_temp[idx[j]]
        DP[i,j]=dist_prob(P[i][j])
        CS[i,j]=cpmt_sgnl(P[i][j])[0]
        P_nor10=cpmt_sgnl(P[i][j])[1]
        P_nor10[P_nor10==0]=np.min(P_nor10[P_nor10>0])/2
        P_nor1[i,j]=np.log2(P_nor10)
        P_nor2[i,j]=np.ma.corrcoef(P_nor1[i][j])
        IS[i,j]=insul_score(P[i][j],bin_size=10**5)[0]

"""Output Data"""
pyw=open(f'{outname}.pyw','w')
for i in range(ncl):
    pyw.write(f'The Contact Probability of {cl[i]}:\n')
    for j in range(nat):
        for k in range(nat):
            pyw.write(f'{P[idx_cl[i]][j,k]}\n')
    pyw.write('\n')
    pyw.write(f'The Obs/Exp Matrix of {cl[i]}:\n')
    for j in range(nat):
        for k in range(nat):
            pyw.write(f'{P_nor1[idx_cl[i]][j,k]}\n')
    pyw.write('\n')
    pyw.write(f'The Correlation Coefficient Matrix of {cl[i]}:\n')
    for j in range(nat):
        for k in range(nat):
            pyw.write(f'{P_nor2[idx_cl[i]][j,k]}\n')
    pyw.write('\n')
    pyw.write(f'The Contact Probability v.s. Sequence Distance of {cl[i]}:\n')
    for j in range(nat):
        pyw.write(f'{DP[idx_cl[i]][j]}\n')
    pyw.write('\n')
    pyw.write(f'The Compartment Signal of {cl[i]}:\n')
    for j in range(nat):
        pyw.write(f'{CS[idx_cl[i]][j]}\n')
    pyw.write('\n')
    pyw.write(f'The Insulation Score of {cl[i]}:\n')
    for j in range(nat):
        pyw.write(f'{IS[idx_cl[i]][j]}\n')
    pyw.write('\n')

for i in range(nsw):
    for j in range(1,len(t)+1):
        pyw.write(f'The Contact Probability of {sw[i]} at t={t[j-1]}:\n')
        for k in range(nat):
            for l in range(nat):
                pyw.write(f'{P[i,j,k,l]}\n')
        pyw.write('\n')
        pyw.write(f'The Obs/Exp Matrix of {sw[i]} at t={t[j-1]}:\n')
        for k in range(nat):
            for l in range(nat):
                pyw.write(f'{P_nor1[i,j,k,l]}\n')
        pyw.write('\n')
        pyw.write(f'The Correlation Coefficient Matrix of {sw[i]} at t={t[j-1]}:\n')
        for k in range(nat):
            for l in range(nat):
                pyw.write(f'{P_nor2[i,j,k,l]}\n')
        pyw.write('\n')
        pyw.write(f'The Contact Probability v.s. Sequence Distance of {sw[i]} at t={t[j-1]}:\n')
        for k in range(nat):
            pyw.write(f'{DP[i,j,k]}\n')
        pyw.write('\n')
        pyw.write(f'The Compartment Signal of {sw[i]} at t={t[j-1]}:\n')
        for k in range(nat):
            pyw.write(f'{CS[i,j,k]}\n')
        pyw.write('\n')
        pyw.write(f'The Insulation Score of {sw[i]} at t={t[j-1]}:\n')
        for k in range(nat):
            pyw.write(f'{IS[i,j,k]}\n')
        pyw.write('\n')

pyw.close()


