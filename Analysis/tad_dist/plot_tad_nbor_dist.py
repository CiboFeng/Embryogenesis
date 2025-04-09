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
from tot_len import tot_len

"""Set Aruments"""
file='tad_nbor_dist.pyw'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nhist=50
sgm=4.0
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read Data"""
d0=np.zeros((nsw,ncl,nhist))
P_d0=np.zeros((nsw,ncl,nhist))
d0_mean=np.zeros((nsw,ncl))
d0_std=np.zeros((nsw,ncl))
for i in range(nsw):
    q=pyw(file,f'Probability Density (Reference Value, {sw[i]}):')
    for j in range(ncl):
        for k in range(nhist):
            d0[i,j,k]=q[1][j*nhist+k]/sgm
            P_d0[i,j,k]=q[2][j*nhist+k]*sgm
    q=pyw(file,f'TAD (Reference Value, {sw[i]}):')
    for j in range(ncl):
        d0_mean[i,j]=q[1][j]/sgm
        d0_std[i,j]=q[2][j]/sgm

d=np.zeros((nsw,nfr,nhist))
P_d=np.zeros((nsw,nfr,nhist))
d_mean=np.zeros((nsw,nfr))
d_std=np.zeros((nsw,nfr))
for i in range(nsw):
    q=pyw(file,f'Probability Density ({sw[i]}):')
    for j in range(nfr):
        for k in range(nhist):
            d[i,j,k]=q[1][j*nhist+k]/sgm
            P_d[i,j,k]=q[2][j*nhist+k]*sgm
    q=pyw(file,f'TAD ({sw[i]}):')
    for j in range(nfr):
        d_mean[i,j]=q[1][j]/sgm
        d_std[i,j]=q[2][j]/sgm

"""Plot"""
fig,ax=plt.subplots(2,1,figsize=(10,10))
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
        for k in range(1,nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i].plot(d[i,k,:],P_d[i,k,:],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[i].plot(d0[i,k,:],P_d0[i,k,:],'-',color='k',linewidth=2*wth)
            ax[i].plot(d0[i,k,:],P_d0[i,k,:],'-',color=color0[i][k],linewidth=wth,label=label0[i][k])
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
        ax[i].legend(loc='best',fontsize=0.8*size)
        ax[i].set_ylabel(f'Distribution\nin {sw[i]}',fontsize=size)
ax[1].set_xlabel('Neighbor Distance $(\\sigma)$',fontsize=size)

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
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,ax=plt.subplots(1,1,figsize=(9,5))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=['b','r']

ax2=ax.twinx()
for i in range(nsw):
    ax.fill_between(t[1:-9],d_mean[i,1:-9]-d_std[i,1:-9]/10,d_mean[i,1:-9]+d_std[i,1:-9]/10,color=color[i],alpha=0.1)
    for k in range(ncl):
        if (i,k)!=(0,1):
            ax.plot([min(t[1:-9]),max(t[1:-9])],[d0_mean[i,k],d0_mean[i,k]],'--',color=color0[i][k],linewidth=wth,label=label0[i][k])
    ax2.plot(t[1:-9],d_mean[i,1:-9],'-',color=color[i],linewidth=wth,label=sw[i])
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.set_xscale('log')
# ax.xaxis.set_major_locator(MultipleLocator(xtick))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_major_locator(MultipleLocator(ytick))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
ymax=max(np.max(d_mean[:,1:-9]),np.max(d0_mean))
ymin=min(np.min(d_mean[:,1:-9]),np.min(d0_mean))
ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
ax.legend(loc='lower right',fontsize=0.8*size)
ax.set_xlabel(f'Time $(\\tau)$',fontsize=size)
ax.set_ylabel('$\\langle d_{\\mathrm{n}} \\rangle ~(\\sigma)$',fontsize=size)
ax2.set_yticks([])
ax2.set_ylim(ax.get_ylim())
ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax2.legend(loc='upper right',fontsize=0.8*size)

plt.tight_layout()
plt.savefig(f'{figname}_mean.png', format='png')
plt.savefig(f'{figname}_mean.pdf', format='pdf')
plt.show()