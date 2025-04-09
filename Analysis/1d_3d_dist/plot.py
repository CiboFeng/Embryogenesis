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
file='1d_3d_dist.npz'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
sgm=4.0
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read Data"""
data=np.load(file)
d0_mean=data['d0_mean']/sgm
d0_std=data['d0_std']/sgm
d_mean=data['d_mean']/sgm
d_std=data['d_std']/sgm

l=np.arange(0,nat/10,0.1)

"""Plot"""
fig,ax=plt.subplots(2,3,figsize=(20,10))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))

for i in range(nsw):
    for j in range(1,nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,0].plot(l,d_mean[i,j,:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i,0].plot(l,d0_mean[i,j,:],'-',color='k',linewidth=2*wth)
        ax[i,0].plot(l,d0_mean[i,j,:],'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
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
    ax[i,0].legend(loc='best',fontsize=0.8*size)
    ax[i,0].set_ylabel('$d(l)~(\\sigma)$\nin %s'%sw[i],fontsize=size)

    for j in range(1,nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,1].plot(l,d_mean[i,j,:]-l**0.5,'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i,1].plot(l,d0_mean[i,j,:]-l**0.5,'-',color='k',linewidth=2*wth)
        ax[i,1].plot(l,d0_mean[i,j,:]-l**0.5,'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
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
    ax[i,1].legend(loc='best',fontsize=0.8*size)
    ax[i,1].set_ylabel('$d(l)-l^{1/2}/2~(\\sigma)$\nin %s'%sw[i],fontsize=size)

    for j in range(1,nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,2].plot(l[1:],d_mean[i,j,1:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i,2].plot(l[1:],d0_mean[i,j,1:],'-',color='k',linewidth=2*wth)
        ax[i,2].plot(l[1:],d0_mean[i,j,1:],'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
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
    ax[i,2].legend(loc='best',fontsize=0.8*size)
    ax[i,2].set_xscale('log')
    ax[i,2].set_yscale('log')
    ax[i,2].set_ylabel('$d(l)~(\\sigma)$\nin %s'%sw[i],fontsize=size)
ax[1,0].set_xlabel('Sequence Distance (Mb)',fontsize=size)
ax[1,1].set_xlabel('Sequence Distance (Mb)',fontsize=size)
ax[1,2].set_xlabel('Sequence Distance (Mb)',fontsize=size)

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$t~(\\tau)$',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0,0.85,1])
plt.savefig(f'{figname}1.png',format='png')
plt.savefig(f'{figname}1.pdf',format='pdf')
plt.show()

fig,ax=plt.subplots(2,1,figsize=(10,10))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(1,0.8,0),(1,0,0),(0,0,0)]
nodes=[0/9,4/9,9/9]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=0,vmax=nat/10)

for i in range(nsw):
    # for j in range(ntad-1,0,-1):
    #     normc=norm(j-1)
    #     rgb=cmap(normc)
    #     ax[i].plot(t[1:-9],r_mean[i,1:-9,j],'-',color=rgb,linewidth=wth)
        # ax[i].plot(r0_mean[i,k,:],'-',color='k',linewidth=2*wth)
        # ax[i].plot(r0_mean[i,k,:],'-',color=color0[k][i],linewidth=wth,label=label0[k][i])
    # for j in range(nat):
    for j in [0,2,4,8,16,32,64,128,256,512]:
        normc=norm(j/10)
        rgb=cmap(normc)
        ax[i].plot(t[1:-9],d_mean[i,1:-9,j]-d_mean[i,-10,j],'-',color=rgb,linewidth=wth)
    ax[i].minorticks_on()
    ax[i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i].spines['bottom'].set_linewidth(wth)
    ax[i].spines['top'].set_linewidth(wth)
    ax[i].spines['left'].set_linewidth(wth)
    ax[i].spines['right'].set_linewidth(wth)
    ax[i].set_xscale('log')
    # ax[i].set_yscale('log')
    # ax[i].legend(loc='best',fontsize=0.8*size)
    ax[i].set_ylabel('$d(l)-d(l)^\\mathrm{ESC}~(\\sigma)$\nin %s'%sw[i],fontsize=size)
ax[1].set_xlabel('Time $(\\tau)$',fontsize=size)

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Sequence Distance (Mb)',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0,0.85,1])
plt.savefig(f'{figname}2.png',format='png')
plt.savefig(f'{figname}2.pdf',format='pdf')
plt.show()

# fig,ax=plt.subplots(1,1,figsize=(9,5))
# wth=3
# size=30
# lenmaj=15
# lenmin=8
# xtick=1
# ytick=0.02
# color0=[[(1.0,0.8,0.6),(0.8,0.6,1.0)],[(0.8,1.0,0.6),(0.8,1.0,0.6)]]
# label0=[['ZP','ZM'],['ESC','ESC']]
# color=['b','r']
#
# for j in range(nsw):
#     for k in range(ncl):
#         if (j,k)!=(0,1):
#             ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[j,k],r0_mean[j,k]],'--',color=color0[k][j],linewidth=wth,label=label0[k][j])
#     ax.plot(t[1:-9],r_mean[j,1:-9],'-',color=color[j],linewidth=wth)
#     # ax.errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.set_xscale('log')
# # ax.xaxis.set_major_locator(MultipleLocator(xtick))
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# # ax.yaxis.set_major_locator(MultipleLocator(ytick))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# ax.legend(loc='lower right',fontsize=0.8*size)
# ax.set_xlabel(f'Time $(\\tau)$',fontsize=size)
# ax.set_ylabel('$\\langle R_{\\mathrm{g}} \\rangle ~(\\sigma)$',fontsize=size)
# ax.set_ylabel('$\\langle d_{\\mathrm{n}} \\rangle ~(\\sigma)$',fontsize=size)
#
# plt.tight_layout()
# plt.savefig(f'{figname}_mean.png', format='png')
# plt.savefig(f'{figname}_mean.pdf', format='pdf')
# plt.show()