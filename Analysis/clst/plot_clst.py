"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
from pyw import pyw
from dist_prob import dist_prob
from cpmt_sgnl import cpmt_sgnl
from insul_score import insul_score

"""Set Arguments"""
file1='clst1.pyw'
file2='dim_rdc.pyw'
file3='clst2.pyw'
file4=['cont_prob_zp.dat','cont_prob_zm.dat','cont_prob_esc.dat']
figname1=file1[:-5]
figname2='repro'
cls={'ZP':0.0,'ZM':1.0,'ESC':2.0}
ncl=3
sgm=4
nat=600
nhist=20

"""Read Data"""
q=pyw(file1,dir=True,mat=False)
Z=[]
for i in range(ncl):
    z=np.zeros((int(np.shape(q)[1]/4),4))
    for j in range(int(len(q[i])/4)):
        for k in range(4):
            z[j,k]=q[i][j*4+k][0]
    z[:,2]/=sgm
    Z+=[z]

q=pyw(file2,dir=True)
x=q[:,:,0]
y=q[:,:,1]

idx=[]
DD=np.zeros((ncl,nhist))
PP=np.zeros((ncl,nhist))
P_rdc=np.ones((ncl,nat,nat))
dp_rdc=np.zeros((ncl,nat))
cs_rdc=np.zeros((ncl,nat))
is_rdc=np.zeros((ncl,nat))
for i in range(ncl):
    idx+=pyw(file3,f'Index of {list(cls.keys())[i]}:',fmt='int')
    DD[i]=np.array(pyw(file3,f'{list(cls.keys())[i]},')[0])/sgm
    PP[i]=np.array(pyw(file3,f'{list(cls.keys())[i]},')[1])
    q=np.array(pyw(file3,f'Conformations of {list(cls.keys())[i]}:')[0])
    for j in range(nat):
        for k in range(j+2,nat):
            P_rdc[i,j,k]=q[j*nat+k-int(j*(j+5)/2+2)]
            P_rdc[i,k,j]=q[j*nat+k-int(j*(j+5)/2+2)]
    dp_rdc[i]=dist_prob(P_rdc[i])
    cs_rdc[i]=cpmt_sgnl(P_rdc[i])[0]
    is_rdc[i]=insul_score(P_rdc[i],bin_size=10**5)[0]

P_ful=np.ones((ncl,nat,nat))
dp_ful=np.zeros((ncl,nat))
cs_ful=np.zeros((ncl,nat))
is_ful=np.zeros((ncl,nat))
r2_p=np.zeros(ncl)
r2_dp=np.zeros(ncl)
r2_cs=np.zeros(ncl)
r2_is=np.zeros(ncl)
for i in range(ncl):
    f=open(file4[i],'r')
    lsf=f.readlines()
    f.close()
    for j in range(len(lsf)):
        ls_f=lsf[j].strip('\n').split()
        P_ful[i,int(float(ls_f[0])-1),int(float(ls_f[1])-1)]=float(ls_f[2])
        P_ful[i,int(float(ls_f[1])-1),int(float(ls_f[0])-1)]=float(ls_f[2])
    dp_ful[i]=dist_prob(P_ful[i])
    cs_ful[i]=cpmt_sgnl(P_ful[i])[0]
    is_ful[i]=insul_score(P_ful[i],bin_size=10**5)[0]
    r2_p[i]=r2_score(P_ful[i].reshape(-1),P_rdc[i].reshape(-1))
    r2_dp[i]=r2_score(dp_ful[i],dp_rdc[i])
    r2_cs[i]=r2_score(cs_ful[i],cs_rdc[i])
    r2_is[i]=r2_score(is_ful[i],is_rdc[i])

plt.figure()
norm=colors.LogNorm()
plt.imshow(np.abs(P_ful[0]-P_rdc[0]),norm=norm)
plt.show()

fig,ax=plt.subplots(3,3,figsize=(12,8))
for i in range(ncl):
    ax[i,0].plot(dp_ful[i])
    ax[i,0].plot(dp_rdc[i])
    ax[i,0].set_xscale('log')
    ax[i,0].set_yscale('log')
    ax[i,1].plot(cs_ful[i])
    ax[i,1].plot(cs_rdc[i])
    ax[i,2].plot(is_ful[i])
    ax[i,2].plot(is_rdc[i])
plt.savefig('sgnl.png',format='png')
plt.show()

"""Plot"""
# plt.figure(figsize=(14,7))
# dendrogram(Z[0])
# plt.xticks([])
# plt.savefig('valid.png',format='png')
# plt.show()

fig,ax=plt.subplots(3,3,figsize=(12,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002

for i in range(ncl):
    dendrogram(Z[i],ax=ax[i,0])
    ax[i,0].autoscale()
    # ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].spines['bottom'].set_linewidth(wth)
    ax[i,0].spines['top'].set_linewidth(wth)
    ax[i,0].spines['left'].set_linewidth(wth)
    ax[i,0].spines['right'].set_linewidth(wth)
    ax[i,0].set_xticks([])
    ax[i,0].set_ylabel(f'RMSD ({list(cls.keys())[i]}) /$\\sigma$',fontsize=size)

    ax[i,1].plot(x[i],y[i],'.',color='c',linewidth=wth,alpha=0.5)
    ax[i,1].plot(x[i,idx[i]],y[i,idx[i]],'.',color='k',linewidth=5*wth)
    ax[i,1].autoscale()
    ax[i,1].minorticks_on()
    ax[i,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].spines['bottom'].set_linewidth(wth)
    ax[i,1].spines['top'].set_linewidth(wth)
    ax[i,1].spines['left'].set_linewidth(wth)
    ax[i,1].spines['right'].set_linewidth(wth)
    ax[i,1].set_ylabel(f'PC2 ({list(cls.keys())[i]})',fontsize=size)

    ax[i,2].bar(DD[i],PP[i],color='skyblue',edgecolor='black',width=DD[i,1]-DD[i,0])
    ax[i,2].autoscale()
    ax[i,2].minorticks_on()
    ax[i,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,2].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].spines['bottom'].set_linewidth(wth)
    ax[i,2].spines['top'].set_linewidth(wth)
    ax[i,2].spines['left'].set_linewidth(wth)
    ax[i,2].spines['right'].set_linewidth(wth)
    ax[i,2].set_ylabel(f'Population ({list(cls.keys())[i]})',fontsize=size)

ax[2,0].set_xlabel('Cluster Index',fontsize=size)
ax[2,1].set_xlabel('PC1',fontsize=size)
ax[2,2].set_xlabel('RMSD /$\\sigma$',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname1}.png',format='png')
plt.savefig(f'{figname1}.pdf',format='pdf')
# plt.show()

fig,ax=plt.subplots(3,4,figsize=(12,8))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002

for i in range(ncl):
    ax[i,0].plot(P_ful[i].reshape(-1),P_rdc[i].reshape(-1),'.',color='c',linewidth=wth,alpha=0.5,label=f'r2={format(r2_p[i],".3f")}')
    ax[i,0].plot([np.min(P_ful[i]),np.max(P_ful[i])],[np.min(P_ful[i]),np.max(P_ful[i])],color='k',linewidth=wth)
    ax[i,0].autoscale()
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].spines['bottom'].set_linewidth(wth)
    ax[i,0].spines['top'].set_linewidth(wth)
    ax[i,0].spines['left'].set_linewidth(wth)
    ax[i,0].spines['right'].set_linewidth(wth)
    ax[i,0].set_ylabel(f'Rdc. CP ({list(cls.keys())[i]})',fontsize=size)
    ax[i,0].legend(loc='best')

    ax[i,1].plot(dp_ful[i],dp_rdc[i],'.',color='c',linewidth=wth,alpha=0.5,label=f'r2={format(r2_dp[i],".3f")}')
    ax[i,1].plot([np.min(dp_ful[i]),np.max(dp_ful[i])],[np.min(dp_ful[i]),np.max(dp_ful[i])],color='k',linewidth=wth)
    ax[i,1].autoscale()
    ax[i,1].minorticks_on()
    ax[i,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].spines['bottom'].set_linewidth(wth)
    ax[i,1].spines['top'].set_linewidth(wth)
    ax[i,1].spines['left'].set_linewidth(wth)
    ax[i,1].spines['right'].set_linewidth(wth)
    ax[i,1].set_ylabel(f'Rdc. DP ({list(cls.keys())[i]})',fontsize=size)
    ax[i,1].legend(loc='best')

    ax[i,2].plot(cs_ful[i],cs_rdc[i],'.',color='c',linewidth=wth,alpha=0.5,label=f'r2={format(r2_cs[i],".3f")}')
    ax[i,2].plot([np.min(cs_ful[i]),np.max(cs_ful[i])],[np.min(cs_ful[i]),np.max(cs_ful[i])],color='k',linewidth=wth)
    ax[i,2].autoscale()
    ax[i,2].minorticks_on()
    ax[i,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,2].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].spines['bottom'].set_linewidth(wth)
    ax[i,2].spines['top'].set_linewidth(wth)
    ax[i,2].spines['left'].set_linewidth(wth)
    ax[i,2].spines['right'].set_linewidth(wth)
    ax[i,2].set_ylabel(f'Rdc. CS ({list(cls.keys())[i]})',fontsize=size)
    ax[i,2].legend(loc='best')

    ax[i,3].plot(is_ful[i],is_rdc[i],'.',color='c',linewidth=wth,alpha=0.5,label=f'r2={format(r2_is[i],".3f")}')
    ax[i,3].plot([np.min(is_ful[i]),np.max(is_ful[i])],[np.min(is_ful[i]),np.max(is_ful[i])],color='k',linewidth=wth)
    ax[i,3].autoscale()
    ax[i,3].minorticks_on()
    ax[i,3].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,3].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,3].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,3].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,3].spines['bottom'].set_linewidth(wth)
    ax[i,3].spines['top'].set_linewidth(wth)
    ax[i,3].spines['left'].set_linewidth(wth)
    ax[i,3].spines['right'].set_linewidth(wth)
    ax[i,3].set_ylabel(f'Rdc. IS ({list(cls.keys())[i]})',fontsize=size)
    ax[i,3].legend(loc='best')

ax[2,0].set_xlabel('Ful. CP',fontsize=size)
ax[2,1].set_xlabel('Ful. DP',fontsize=size)
ax[2,2].set_xlabel('Ful. CS',fontsize=size)
ax[2,3].set_xlabel('Ful. IS',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname2}.png',format='png')
plt.savefig(f'{figname2}.pdf',format='pdf')
# plt.show()
