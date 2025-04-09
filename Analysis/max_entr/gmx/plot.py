import numpy as np
import matplotlib.pyplot as plt

cp_gd_hic_file='CP_GD_HiC.dat'
cp_gd_sim_file='CP_GD_Sim.dat'
cm_hic_file='CM_HiC.dat'
cm_sim_file='CM_Sim.dat'
insulationscore_hic_file='InsulationScore_HiC.dat'
insulationscore_sim_file='InsulationScore_Sim.dat'

d_exp=[]
cd_exp=[]
f=open(cp_gd_hic_file)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
for i in range(len(ls_f)-1):
    d_exp+=[float(ls_f[i][0])]
    cd_exp+=[float(ls_f[i][1])]

d_sim=[]
cd_sim=[]
f=open(cp_gd_sim_file)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
for i in range(len(ls_f)-1):
    d_sim+=[float(ls_f[i][0])]
    cd_sim+=[float(ls_f[i][1])]

idx1_exp=[]
c_exp=[]
f=open(cm_hic_file)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
for i in range(len(ls_f)):
    idx1_exp+=[10*(i+1)]
    c_exp+=[float(ls_f[i][0])]

idx1_sim=[]
c_sim=[]
f=open(cm_sim_file)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
for i in range(len(ls_f)):
    idx1_sim+=[10*(i+1)]
    c_sim+=[float(ls_f[i][0])]

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

# Compare=np.zeros((N,N))
# for i in range(N):
#     for j in range(N):
#         if i>j:
#             Compare[i,j]=P_exp[i,j]
#         else:
#             Compare[i,j]=P_sim[i,j]
# Compare[Compare==0]=np.min(Compare[Compare>0])/2
#
# plt.figure()
# norm=colors.SymLogNorm(0.001)
# plt.imshow(Alf_o,norm=norm)
# plt.gca().invert_yaxis()
# plt.xlabel('Index of particles')
# plt.ylabel('Index of particles')
# plt.colorbar()
# plt.savefig(opath+f'alf_{sysname}.png',format='png')
# plt.savefig(opath+f'alf_{sysname}.pdf',format='pdf')
#
# plt.figure()
# norm=colors.LogNorm()
# plt.imshow(Compare,norm=norm)
# plt.gca().invert_yaxis()
# plt.title('Upper left: Exp, Lower right: Sim')
# plt.xlabel('Index of particles')
# plt.ylabel('Index of particles')
# plt.colorbar()
# plt.savefig(opath+f'cont_prob_{sysname}.png',format='png')
# plt.savefig(opath+f'cont_prob_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(d_exp,cd_exp,label='Exp')
plt.plot(d_sim,cd_sim,label='Sim')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Sequence distance')
plt.ylabel('Contact probability')
plt.savefig(f'cont_dist.pdf',format='pdf')

plt.figure()
plt.plot(idx1_exp,c_exp,label='Exp')
plt.plot(idx1_sim,c_sim,label='Sim')
plt.legend(loc='best')
plt.xlabel('Index of particles')
plt.ylabel('Compartment signal')
plt.savefig(f'cmpt.pdf',format='pdf')

plt.figure()
plt.plot(idx2_exp,is_exp,label='Exp')
plt.plot(idx2_sim,is_sim,label='Sim')
plt.legend(loc='best')
plt.xlabel('Index of particles')
plt.ylabel('Insulation score')
plt.savefig(f'insul_score.pdf',format='pdf')
