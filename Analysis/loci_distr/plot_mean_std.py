"""Import Modules"""
import numpy as np
from scipy import interpolate
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='loci_distr.pyw'
figname='mean_std'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl=len(cl)
nfr=56
ts=-2
te=3
nt=10000
rmax=52.3424
sgm=4.0
idx1=[0,1,3,5]
idx2=[2,4,6]

"""Read Data"""
r0_mean=np.zeros((nsw,ncl,7))
r0_std=np.zeros((nsw,ncl,7))
for i in range(nsw):
    q=pyw(file,f'The Index of Cells, Mean, and Standard Deviation of Location of All Loci, of Location of Loci in '
               f'Compartment A, of Relative Location of Loci in Compartment A, Similar for Loci in Compartment B, and '
               f'Weighted by Gene Density ({sw[i]}):')
    for j in range(ncl):
        for k in idx1:
            r0_mean[i,j,k]=q[2*k+1][j]/rmax
            r0_std[i,j,k]=q[2*k+2][j]/rmax
        for k in idx2:
            r0_mean[i,j,k]=q[2*k+1][j]
            r0_std[i,j,k]=q[2*k+2][j]

r_mean=np.zeros((nsw,nfr,7))
r_std=np.zeros((nsw,nfr,7))
for i in range(nsw):
    q=pyw(file,f'Time, Mean, and Standard Deviation of Location of All Loci, of Location of Loci in Compartment A, of '
               f'Relative Location of Loci in Compartment A, Similar for Loci in Compartment B, and Weighted by Gene '
               f'Density ({sw[i]}):')
    for j in range(nfr):
        for k in idx1:
            r_mean[i,j,k]=q[2*k+1][j]/rmax
            r_std[i,j,k]=q[2*k+2][j]/rmax
        for k in idx2:
            r_mean[i,j,k]=q[2*k+1][j]
            r_std[i,j,k]=q[2*k+2][j]

"""Generate Times"""
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
t_=np.logspace(ts,te,nt)

"""Define a Function for Interpolating"""
def itplt(t0,x0,t):
   f=interpolate.interp1d(t0,x0,kind='cubic')
   x=f(t)
   return x

"""Plot"""
fig,ax=plt.subplots(1,len(idx1),figsize=(35,7))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['All','A','B','Gene']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))

for i in range(len(idx1)):
    for j in range(nsw):
        ax[i].scatter(itplt(t,r_mean[j,:,idx1[i]],t_),itplt(t,r_std[j,:,idx1[i]],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
        scatter=ax[i].scatter(itplt(t,r_mean[j,:,idx1[i]],t_),itplt(t,r_std[j,:,idx1[i]],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[i].scatter(r0_mean[j,k,idx1[i]],r0_std[j,k,idx1[i]],color=color0[j][k],marker='^',s=100*wth)
                ax[i].scatter([],[],color=color0[j][k],marker='^',s=100*wth,label=label0[j][k])
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
    ax[i].set_xlabel('$r^\\mathrm{mean}~(R_\\mathrm{m})$ of %s'%label[i],fontsize=size)
    ax[i].set_ylabel('$r^\\mathrm{std}~(R_\\mathrm{m})$ of %s'%label[i],fontsize=size)

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
# plt.show()

fig,ax=plt.subplots(1,len(idx2),figsize=(25,7))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['A','B','Gene']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))

for i in range(len(idx2)):
    for j in range(nsw):
        ax[i].scatter(itplt(t,r_mean[j,:,idx2[i]],t_),itplt(t,r_std[j,:,idx2[i]],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
        scatter=ax[i].scatter(itplt(t,r_mean[j,:,idx2[i]],t_),itplt(t,r_std[j,:,idx2[i]],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[i].scatter(r0_mean[j,k,idx2[i]],r0_std[j,k,idx2[i]],color=color0[j][k],marker='^',s=100*wth)
                ax[i].scatter([],[],color=color0[j][k],marker='^',s=100*wth,label=label0[j][k])
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
    ax[i].set_xlabel('$r_\\mathrm{rel}^\\mathrm{mean}$ of %s'%label[i],fontsize=size)
    ax[i].set_ylabel('$r_\\mathrm{rel}^\\mathrm{std}$ of %s'%label[i],fontsize=size)

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
plt.savefig(f'{figname}_rel.png',format='png')
plt.savefig(f'{figname}_rel.pdf',format='pdf')
# plt.show()