"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from insul_score import insul_score

"""Set Aruments"""
reffile='../cont_prob/ref/esc.dat'
file='dist_distr.npz'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
nhist=20
sgm=4.0
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
# L=[1,
#    2,3,4,5,6,7,8,9,10,
#    20,30,40,50,60,70,80,90,100,
#    200,300,400,500]
L=[1,2,4,8,16,32,64,128,256,512]
cl_idx=[(0,0),(1,0),(0,1)]

"""Read Data"""
data=np.load(file)
d0=data['d0']/sgm
P_d0=data['P_d0']*sgm
d0_mean=data['d0_mean']/sgm
d0_std=data['d0_std']/sgm

l=np.arange(0,nat/10,0.1)

"""Plot"""
fig,ax=plt.subplots(3,3,figsize=(20,10))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,0.6,1.0)],[(0.8,1.0,0.6),(0.8,1.0,0.6)]]
axlabel=['ZP','ZM','ESC']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=np.min(l[1:]),vmax=np.max(l[1:]))

for i in range(len(cl_idx)):
    for j in L:
        normc=norm(l[j])
        rgb=cmap(normc)
        ax[i,0].plot(d0[cl_idx[i]][j],P_d0[cl_idx[i]][j]/np.max(P_d0[cl_idx[i]][j]),'-',color=rgb,linewidth=wth)
        # ax[i,0].plot((d0[cl_idx[i]][j]-d0[cl_idx[i]][j][P_d0[cl_idx[i]][j]==np.max(P_d0[cl_idx[i]][j])])**2,np.log(P_d0[cl_idx[i]][j]),'-',color=rgb,linewidth=wth)
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,0].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i,0].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].spines['bottom'].set_linewidth(wth)
    ax[i,0].spines['top'].set_linewidth(wth)
    ax[i,0].spines['left'].set_linewidth(wth)
    ax[i,0].spines['right'].set_linewidth(wth)
    # ax[i,0].set_yscale('log')
    # ax[i,0].set_xlim(0,10)
    # ax[i,0].set_ylim(0,0.5)
    ax[i,0].set_ylabel('$P(d)~(\\sigma^{-1})$\nof %s'%axlabel[i],fontsize=size)

    ax[i,1].plot(l,d0_mean[cl_idx[i]],'-',color='k',linewidth=wth)
    ax[i,1].minorticks_on()
    ax[i,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,1].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i,1].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i,1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].spines['bottom'].set_linewidth(wth)
    ax[i,1].spines['top'].set_linewidth(wth)
    ax[i,1].spines['left'].set_linewidth(wth)
    ax[i,1].spines['right'].set_linewidth(wth)
    ax[i,1].set_xscale('log')
    ax[i,1].set_yscale('log')
    ax[i,1].set_ylabel('$d_\\mathrm{mean}~(\\sigma)$\nof %s'%axlabel[i],fontsize=size)

    ax[i,2].plot(l,d0_std[cl_idx[i]],'-',color='k',linewidth=wth)
    ax[i,2].minorticks_on()
    ax[i,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,2].xaxis.set_major_locator(MultipleLocator(xtick))
    # ax[i,2].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i,2].yaxis.set_major_locator(MultipleLocator(ytick))
    # ax[i,2].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].spines['bottom'].set_linewidth(wth)
    ax[i,2].spines['top'].set_linewidth(wth)
    ax[i,2].spines['left'].set_linewidth(wth)
    ax[i,2].spines['right'].set_linewidth(wth)
    ax[i,2].set_xscale('log')
    ax[i,2].set_yscale('log')
    ax[i,2].set_ylabel('$d_\\mathrm{std}~(\\sigma)$\nof %s'%axlabel[i],fontsize=size)
ax[2,0].set_xlabel('$d~(\\sigma)$',fontsize=size)
ax[2,1].set_xlabel('$L$ (Mb)',fontsize=size)
ax[2,2].set_xlabel('$L$ (Mb)',fontsize=size)

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$L$ (Mb)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0,0.85,1])
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()