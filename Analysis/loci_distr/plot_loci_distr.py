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
import matplotlib.patheffects as PathEffects
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='loci_distr.pyw'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
nsw=len(sw)
ncl=len(cl)
nfr=56
nhist=50
rmax=52.3424
sgm=4.0
idx1=[0,1,3,5]
idx2=[2,4,6]

"""Read Data"""
r0=np.zeros((nsw,ncl,nhist,7))
P0=np.zeros((nsw,ncl,nhist,7))
r0_mean=np.zeros((nsw,ncl,7))
r0_std=np.zeros((nsw,ncl,7))
for i in range(nsw):
    q=pyw(file,f'The Index of Cells, Location of All Loci, Its Probability Density, Location of Loci in Compartment A, '
               f'Its Probability Density, Relative Location of Loci in Compartment A, Its Probability Density, Similar '
               f'for Loci in Compartment B, and Weighted by Gene Density ({sw[i]}):')
    for j in range(ncl):
        for k in range(nhist):
            for l in idx1:
                r0[i,j,k,l]=q[2*l+1][j*nhist+k]/rmax
                P0[i,j,k,l]=q[2*l+2][j*nhist+k]*rmax
            for l in idx2:
                r0[i,j,k,l]=q[2*l+1][j*nhist+k]
                P0[i,j,k,l]=q[2*l+2][j*nhist+k]
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
# rho0=P0/(4*np.pi*r0**2)
rho0=P0

t=np.zeros(nfr)
r=np.zeros((nsw,nfr,nhist,7))
P=np.zeros((nsw,nfr,nhist,7))
r_mean=np.zeros((nsw,nfr,7))
r_std=np.zeros((nsw,nfr,7))
for i in range(nsw):
    q=pyw(file,f'Time, Location of All Loci, Its Probability Density, Location of Loci in Compartment A, Its Probability '
               f'Density, Relative Location of Loci in Compartment A, Its Probability Density, Similar for Loci in '
               f'Compartment B, and Weighted by Gene Density ({sw[i]}):')
    for j in range(nfr):
        t[j]=q[0][j*nhist]
        for k in range(nhist):
            for l in idx1:
                r[i,j,k,l]=q[2*l+1][j*nhist+k]/rmax
                P[i,j,k,l]=q[2*l+2][j*nhist+k]*rmax
            for l in idx2:
                r[i,j,k,l]=q[2*l+1][j*nhist+k]
                P[i,j,k,l]=q[2*l+2][j*nhist+k]
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
# rho=P/(4*np.pi*r**2)
rho=P

"""Plot"""
fig,ax=plt.subplots(2,4,figsize=(22,10))
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

for i in range(nsw):
    for j in range(len(idx1)):
        for k in range(nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i,j].plot(r[i,k,:,idx1[j]],rho[i,k,:,idx1[j]],color=rgb,linestyle='-',linewidth=wth)
        for k in range(ncl):
            ax[i,j].plot(r0[i,k,:,idx1[j]],rho0[i,k,:,idx1[j]],color='k',linestyle='-',linewidth=2*wth)
            ax[i,j].plot(r0[i,k,:,idx1[j]],rho0[i,k,:,idx1[j]],color=color0[i][k],linestyle='-',linewidth=wth,label=label0[i][k])
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax[i,j].xaxis.set_major_locator(MultipleLocator(xtick))
        ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax[i,j].yaxis.set_major_locator(MultipleLocator(ytick))
        ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        ax[i,0].legend(loc='upper right',fontsize=0.8*size)
ax[1,0].set_xlabel('$r~(R_\\mathrm{m})$ of %s'%label[0],fontsize=size)
ax[1,1].set_xlabel('$r~(R_\\mathrm{m})$ of %s'%label[1],fontsize=size)
ax[1,2].set_xlabel('$r~(R_\\mathrm{m})$ of %s'%label[2],fontsize=size)
ax[1,3].set_xlabel('$r~(R_\\mathrm{m})$ of %s'%label[3],fontsize=size)
ax[0,0].set_ylabel(f'$\\rho$ of {sw[0]}',fontsize=size)
ax[1,0].set_ylabel(f'$\\rho$ of {sw[1]}',fontsize=size)

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

fig,ax=plt.subplots(2,3,figsize=(22,10))
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

for i in range(nsw):
    for j in range(len(idx2)):
        for k in range(nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i,j].plot(r[i,k,:,idx2[j]],rho[i,k,:,idx2[j]],linestyle='-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[i,j].plot(r0[i,k,:,idx2[j]],rho0[i,k,:,idx2[j]],color='k',linestyle='-',linewidth=2*wth)
            ax[i,j].plot(r0[i,k,:,idx2[j]],rho0[i,k,:,idx2[j]],color=color0[i][k],linestyle='-',linewidth=wth,label=label0[i][k])
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax[i,j].xaxis.set_major_locator(MultipleLocator(xtick))
        ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax[i,j].yaxis.set_major_locator(MultipleLocator(ytick))
        ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        ax[i,0].legend(loc='upper right',fontsize=0.8*size,)
ax[1,0].set_xlabel('$r_\\mathrm{rel}$ of %s'%label[0],fontsize=size)
ax[1,1].set_xlabel('$r_\\mathrm{rel}$ of %s'%label[1],fontsize=size)
ax[1,2].set_xlabel('$r_\\mathrm{rel}$ of %s'%label[2],fontsize=size)
ax[0,0].set_ylabel(f'$\\rho$ of {sw[0]}',fontsize=size)
ax[1,0].set_ylabel(f'$\\rho$ of {sw[1]}',fontsize=size)

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
plt.show()

fig,ax=plt.subplots(4,2,figsize=(20,16))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=['r','g','b','m']
label=['All','A','B','Gene']

for i in range(nsw):
    ax2=ax[i,0].twinx()
    for j in range(len(idx1)):
        # ax[i].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i,0].plot([min(t[1:-9]),t0-gap0],[r0_mean[i,k,idx1[j]],r0_mean[i,k,idx1[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
            if k==1:
                ax[i,0].plot([t0+gap0,max(t[1:-9])],[r0_mean[i,k,idx1[j]],r0_mean[i,k,idx1[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
    for k in range(ncl):
        ax[i,0].plot([],[],color=color0[i][k],linestyle='--',linewidth=wth,label=label0[i][k])
    for j in range(len(idx1)):
        ax2.plot(t[1:-9],r_mean[i,1:-9,idx1[j]],c=color[j],linestyle='-',linewidth=wth,label=f'{label[j]}')
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,0].set_xscale('log')
    # ax[i,0].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i,0].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].spines['bottom'].set_linewidth(wth)
    ax[i,0].spines['top'].set_linewidth(wth)
    ax[i,0].spines['left'].set_linewidth(wth)
    ax[i,0].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_mean[i,1:-9,idx1]),np.max(r0_mean[i,:,idx1]))
    ymin=min(np.min(r_mean[i,1:-9,idx1]),np.min(r0_mean[i,:,idx1]))
    ax[i,0].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    ax[i,0].set_ylabel('$r^\\mathrm{mean}~(R_\\mathrm{m})$\nof %s'%sw[i],fontsize=size)
    if i==0:
        ax[i,0].legend(loc='upper right',fontsize=0.8*size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i,0].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    if i==0:
        ax2.legend(loc='lower left',fontsize=0.8*size)
ax[1,0].set_xlabel(f'$t~(\\tau)$',fontsize=size)

color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=['r','g','b','m']
label=['All','A','B','Gene']

for i in range(nsw):
    ax2=ax[i+2,0].twinx()
    for j in range(len(idx1)):
        # ax[i+2,0].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i+2,0].plot([min(t[1:-9]),t0-gap0],[r0_std[i,k,idx1[j]],r0_std[i,k,idx1[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
            if k==1:
                ax[i+2,0].plot([t0+gap0,max(t[1:-9])],[r0_std[i,k,idx1[j]],r0_std[i,k,idx1[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
    for k in range(ncl):
        ax[i+2,0].plot([],[],color=color0[i][k],linestyle='--',linewidth=wth,label=label0[i][k])
    for j in range(len(idx1)):
        ax2.plot(t[1:-9],r_std[i,1:-9,idx1[j]],c=color[j],linestyle='-',linewidth=wth,label=f'{label[j]}')
    ax[i+2,0].minorticks_on()
    ax[i+2,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i+2,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i+2,0].set_xscale('log')
    # ax[i+2,0].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i+2,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i+2,0].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i+2,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i+2,0].spines['bottom'].set_linewidth(wth)
    ax[i+2,0].spines['top'].set_linewidth(wth)
    ax[i+2,0].spines['left'].set_linewidth(wth)
    ax[i+2,0].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_std[i,1:-9,idx1]),np.max(r0_std[i,:,idx1]))
    ymin=min(np.min(r_std[i,1:-9,idx1]),np.min(r0_std[i,:,idx1]))
    ax[i+2,0].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    ax[i+2,0].set_ylabel('$r^\\mathrm{std}~(R_\\mathrm{m})$\nof %s'%sw[i],fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i+2,0].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[3,0].set_xlabel(f'$t~(\\tau)$',fontsize=size)

color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=['g','b','m']
label=['A','B','Gene']

for i in range(nsw):
    ax2=ax[i,1].twinx()
    for j in range(len(idx2)):
        # ax[i,1].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i,1].plot([min(t[1:-9]),t0-gap0],[r0_mean[i,k,idx2[j]],r0_mean[i,k,idx2[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
            if k==1:
                ax[i,1].plot([t0+gap0,max(t[1:-9])],[r0_mean[i,k,idx2[j]],r0_mean[i,k,idx2[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
    for k in range(ncl):
        ax[i,1].plot([],[],color=color0[i][k],linestyle='--',linewidth=wth,label=label0[i][k])
    for j in range(len(idx2)):
        ax2.plot(t[1:-9],r_mean[i,1:-9,idx2[j]],c=color[j],linestyle='-',linewidth=wth,label=f'{label[j]}')
    ax[i,1].minorticks_on()
    ax[i,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,1].set_xscale('log')
    # ax[i,1].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i,1].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i,1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].spines['bottom'].set_linewidth(wth)
    ax[i,1].spines['top'].set_linewidth(wth)
    ax[i,1].spines['left'].set_linewidth(wth)
    ax[i,1].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_mean[i,1:-9,idx2]),np.max(r0_mean[i,:,idx2]))
    ymin=min(np.min(r_mean[i,1:-9,idx2]),np.min(r0_mean[i,:,idx2]))
    ax[i,1].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    ax[i,1].set_ylabel('Mean $r_\\mathrm{rel}^\\mathrm{mean}$\nof %s'%sw[i],fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i,1].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[1,1].set_xlabel(f'$t~(\\tau)$',fontsize=size)

color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=['g','b','m']
label=['A','B','Gene']

for i in range(nsw):
    ax2=ax[i+2,1].twinx()
    for j in range(len(idx2)):
        # ax[i+2,1].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i+2,1].plot([min(t[1:-9]),t0-gap0],[r0_std[i,k,idx2[j]],r0_std[i,k,idx2[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
            if k==1:
                ax[i+2,1].plot([t0+gap0,max(t[1:-9])],[r0_std[i,k,idx2[j]],r0_std[i,k,idx2[j]]],color=color0[i][k],linestyle='--',linewidth=wth)
    for k in range(ncl):
        ax[i+2,1].plot([],[],color=color0[i][k],linestyle='--',linewidth=wth,label=label0[i][k])
    for j in range(len(idx2)):
        ax2.plot(t[1:-9],r_std[i,1:-9,idx2[j]],c=color[j],linestyle='-',linewidth=wth,label=f'{label[j]}')
    ax[i+2,1].minorticks_on()
    ax[i+2,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i+2,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i+2,1].set_xscale('log')
    # ax[i+2,1].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i+2,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i+2,1].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i+2,1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i+2,1].spines['bottom'].set_linewidth(wth)
    ax[i+2,1].spines['top'].set_linewidth(wth)
    ax[i+2,1].spines['left'].set_linewidth(wth)
    ax[i+2,1].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_std[i,1:-9,idx2]),np.max(r0_std[i,:,idx2]))
    ymin=min(np.min(r_std[i,1:-9,idx2]),np.min(r0_std[i,:,idx2]))
    ax[i+2,1].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    ax[i+2,1].set_ylabel('$r_\\mathrm{rel}^\\mathrm{std}$\nof %s'%sw[i],fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i+2,1].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax[3,1].set_xlabel(f'$t~(\\tau)$',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_mean_std.png', format='png')
plt.savefig(f'{figname}_mean_std.pdf', format='pdf')
plt.show()
