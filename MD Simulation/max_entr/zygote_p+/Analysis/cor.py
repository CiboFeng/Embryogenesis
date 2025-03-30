import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

pyw=open('cor.pyw','w')
pyw.write(f'alf dp is cs:\n')
pyw.close()
a=0
while a>=0:
    cp_gd_hic_file=f'../ConstantT/a_{a}/120.28/Probability/Analysis/CP_GD_HiC.dat'
    cp_gd_sim_file=f'../ConstantT/a_{a}/120.28/Probability/Analysis/CP_GD_Sim.dat'
    cm_hic_file=f'../ConstantT/a_{a}/120.28/Probability/Analysis/CM_HiC.dat'
    cm_sim_file=f'../ConstantT/a_{a}/120.28/Probability/Analysis/CM_Sim.dat'
    insulationscore_hic_file=f'../ConstantT/a_{a}/120.28/Probability/Analysis/InsulationScore_HiC.dat'
    insulationscore_sim_file=f'../ConstantT/a_{a}/120.28/Probability/Analysis/InsulationScore_Sim.dat'

    if not os.path.exists(cp_gd_hic_file):
        break

    d_exp=[]
    dp_exp=[]
    f=open(cp_gd_hic_file)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    for i in range(len(ls_f)-1):
        d_exp+=[float(ls_f[i][0])]
        dp_exp+=[float(ls_f[i][1])]

    d_sim=[]
    dp_sim=[]
    f=open(cp_gd_sim_file)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    for i in range(len(ls_f)-1):
        d_sim+=[float(ls_f[i][0])]
        dp_sim+=[float(ls_f[i][1])]

    r2_dp=r2_score(dp_exp,dp_sim)

    idx1_exp=[]
    cs_exp=[]
    f=open(cm_hic_file)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    for i in range(len(ls_f)):
        idx1_exp+=[10*(i+1)]
        cs_exp+=[float(ls_f[i][0])]

    idx1_sim=[]
    cs_sim=[]
    f=open(cm_sim_file)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    for i in range(len(ls_f)):
        idx1_sim+=[10*(i+1)]
        cs_sim+=[float(ls_f[i][0])]
    cs_cor=np.corrcoef(cs_sim,cs_exp)[0,1]
    if cs_cor<0:
        cs_sim=-1*np.array(cs_sim)

    r2_cs=r2_score(cs_exp,cs_sim)

    idx2_exp=[]
    is_exp=[]
    f=open(insulationscore_hic_file)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    for i in range(len(ls_f)):
        idx2_exp+=[float(ls_f[i][0])]
        is_exp+=[float(ls_f[i][1])]

    idx2_sim=[]
    is_sim=[]
    f=open(insulationscore_sim_file)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    for i in range(len(ls_f)):
        idx2_sim+=[float(ls_f[i][0])]
        is_sim+=[float(ls_f[i][1])]

    r2_is=r2_score(is_exp,is_sim)

    pyw=open('cor.pyw','a')
    pyw.write(f'{a} {r2_dp} {r2_is} {r2_cs}\n')
    pyw.close()

    a+=1
