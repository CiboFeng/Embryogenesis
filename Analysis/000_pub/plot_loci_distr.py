"""Import Modules"""
import numpy as np
import seaborn as sns
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='../loci_distr/loci_distr.pyw'
figname='loci_distr'
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
ts=-2
te=3
nt=10000

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
fig,ax=plt.subplots(3,4,figsize=(25,15))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=[[0.5,0.5,0.5,0.5],[0.5,0.5,0.5,0.5],[0.05,0.05,0.05,0.05]]
ytick=[[1,1,1,1],[1,1,1,1],[0.01,0.01,0.01,0.01]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['All','A','B','Gene']
axlabel=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))
ann=['(A)','(D)','(G)','(J)','(B)','(E)','(H)','(K)','(C)','(F)','(I)','(L)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

for i in range(nsw):
    lines0=[]
    for j in range(len(idx1)):
        for k in range(nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i,j].plot(r[i,k,:,idx1[j]],rho[i,k,:,idx1[j]],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[i,j].plot(r0[i,k,:,idx1[j]],rho0[i,k,:,idx1[j]],color=(0.5,0.5,0.5),linestyle='-',linewidth=2*wth)
            lines0+=ax[i,j].plot(r0[i,k,:,idx1[j]],rho0[i,k,:,idx1[j]],color=color0[i][k],linestyle='-',linewidth=wth,label=label0[i][k])
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax[i,j].xaxis.set_major_locator(MultipleLocator(xtick[i][j]))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(xtick[i][j]/2))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(ytick[i][j]))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(ytick[i][j]/2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        legend=ax[i,0].legend(loc='upper right',handlelength=0.0,handletextpad=0.0,fontsize=0.8*size)
        for text,line in zip(legend.get_texts(),lines0):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
            text.set_color(line.get_color())
        ax[i,j].set_xlabel('$r/R_\\mathrm{m}$ of %s'%label[j],fontsize=size)
        ax[i,j].set_ylabel('$f(r/R_\\mathrm{m})$\nin %s'%axlabel[i],fontsize=size)

for i in range(len(idx1)):
    scatters0=[]
    for j in range(nsw):
        ax[2,i].scatter(itplt(t,r_mean[j,:,idx1[i]],t_),itplt(t,r_std[j,:,idx1[i]],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
        scatter=ax[2,i].scatter(itplt(t,r_mean[j,:,idx1[i]],t_),itplt(t,r_std[j,:,idx1[i]],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[2,i].scatter(r0_mean[j,k,idx1[i]],r0_std[j,k,idx1[i]],color=(0.5,0.5,0.5),marker='^',s=200*wth)
                scatters0+=[ax[2,i].scatter(r0_mean[j,k,idx1[i]],r0_std[j,k,idx1[i]],color=color0[j][k],marker='^',s=100*wth,label=label0[j][k])]
    ax[2,i].minorticks_on()
    ax[2,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[2,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[2,i].xaxis.set_major_locator(MultipleLocator(xtick[2][i]))
    ax[2,i].xaxis.set_minor_locator(MultipleLocator(xtick[2][i]/2))
    ax[2,i].yaxis.set_major_locator(MultipleLocator(ytick[2][i]))
    ax[2,i].yaxis.set_minor_locator(MultipleLocator(ytick[2][i]/2))
    ax[2,i].spines['bottom'].set_linewidth(wth)
    ax[2,i].spines['top'].set_linewidth(wth)
    ax[2,i].spines['left'].set_linewidth(wth)
    ax[2,i].spines['right'].set_linewidth(wth)
    ax[2,i].set_xlabel('$\\langle{}r/R_\\mathrm{m}\\rangle{}$ of %s'%label[i],fontsize=size)
    ax[2,i].set_ylabel('$\\langle{}\\langle{}r/R_\\mathrm{m}\\rangle{}\\rangle{}$ of %s'%label[i],fontsize=size)

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

fig,ax=plt.subplots(3,3,figsize=(25,15))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=[[1,1,1],[1,1,1],[0.02,0.02,0.02]]
ytick=[[1,1,1],[1,1,1],[0.05,0.02,0.05]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['A','B','Gene']
axlabel=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))
ann=['(A)','(D)','(G)','(B)','(E)','(H)','(C)','(F)','(I)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

for i in range(nsw):
    lines0=[]
    for j in range(len(idx2)):
        for k in range(nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i,j].plot(r[i,k,:,idx2[j]],rho[i,k,:,idx2[j]],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[i,j].plot(r0[i,k,:,idx2[j]],rho0[i,k,:,idx2[j]],color=(0.5,0.5,0.5),linestyle='-',linewidth=2*wth)
            lines0+=ax[i,j].plot(r0[i,k,:,idx2[j]],rho0[i,k,:,idx2[j]],color=color0[i][k],linestyle='-',linewidth=wth,label=label0[i][k])
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax[i,j].xaxis.set_major_locator(MultipleLocator(xtick[i][j]))
        ax[i,j].xaxis.set_minor_locator(MultipleLocator(xtick[i][j]/2))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(ytick[i][j]))
        ax[i,j].yaxis.set_minor_locator(MultipleLocator(ytick[i][j]/2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        legend=ax[i,0].legend(loc='upper right',handlelength=0.0,handletextpad=0.0,fontsize=0.8*size,)
        for text,line in zip(legend.get_texts(),lines0):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
            text.set_color(line.get_color())
        ax[i,j].set_xlabel('$r/R_\\mathrm{g}$ of %s'%label[j],fontsize=size)
        ax[i,j].set_ylabel('$f(r/R_\\mathrm{g})$\nin %s'%axlabel[i],fontsize=size)

for i in range(len(idx2)):
    scatters0=[]
    for j in range(nsw):
        ax[2,i].scatter(itplt(t,r_mean[j,:,idx2[i]],t_),itplt(t,r_std[j,:,idx2[i]],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
        scatter=ax[2,i].scatter(itplt(t,r_mean[j,:,idx2[i]],t_),itplt(t,r_std[j,:,idx2[i]],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[2,i].scatter(r0_mean[j,k,idx2[i]],r0_std[j,k,idx2[i]],color=(0.5,0.5,0.5),marker='^',s=200*wth)
                scatters0+=[ax[2,i].scatter(r0_mean[j,k,idx2[i]],r0_std[j,k,idx2[i]],color=color0[j][k],marker='^',s=100*wth,label=label0[j][k])]
    ax[2,i].minorticks_on()
    ax[2,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[2,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[2,i].xaxis.set_major_locator(MultipleLocator(xtick[2][i]))
    ax[2,i].xaxis.set_minor_locator(MultipleLocator(xtick[2][i]/2))
    ax[2,i].yaxis.set_major_locator(MultipleLocator(ytick[2][i]))
    ax[2,i].yaxis.set_minor_locator(MultipleLocator(ytick[2][i]/2))
    ax[2,i].spines['bottom'].set_linewidth(wth)
    ax[2,i].spines['top'].set_linewidth(wth)
    ax[2,i].spines['left'].set_linewidth(wth)
    ax[2,i].spines['right'].set_linewidth(wth)
    ax[2,i].set_xlabel('$\\langle{}r/R_\\mathrm{g}\\rangle{}$ of %s'%label[i],fontsize=size)
    ax[2,i].set_ylabel('$\\langle{}\\langle{}r/R_\\mathrm{g}\\rangle{}\\rangle{}$ of %s'%label[i],fontsize=size)

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

fig,ax=plt.subplots(2,4,figsize=(25,10))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=[]
ytick=[[0.1,0.01,0.1,0.05],[0.1,0.01,0.1,0.05]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
ann=['(A)','(C)','(E)','(G)','(B)','(D)','(F)','(H)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

for i in range(nsw):
    color=[(1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.9,1.0),(1.0,0.0,1.0)]
    label=['All','A','B','Gene']
    axlabel=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
    ax2=ax[i,0].twinx()
    lines0=[]
    lines=[]
    for j in range(len(idx1)):
        # ax[i,0].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i,0].plot([min(t[1:-9]),t0-gap0],[r0_mean[i,k,idx1[j]],r0_mean[i,k,idx1[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,0].plot([min(t[1:-9]),t0-gap0],[r0_mean[i,k,idx1[j]],r0_mean[i,k,idx1[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
            if k==1:
                ax[i,0].plot([t0+gap0,max(t[1:-9])],[r0_mean[i,k,idx1[j]],r0_mean[i,k,idx1[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,0].plot([t0+gap0,max(t[1:-9])],[r0_mean[i,k,idx1[j]],r0_mean[i,k,idx1[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
    for k in range(ncl):
        lines0+=ax[i,0].plot([],[],color=color0[i][k],linestyle='--',linewidth=2*wth,label=label0[i][k])
    for j in range(len(idx1)):
        ax2.plot(t[1:-9],r_mean[i,1:-9,idx1[j]],'-',c='k',linewidth=3*wth)
        lines+=ax2.plot(t[1:-9],r_mean[i,1:-9,idx1[j]],'-',c=color[j],linewidth=2*wth,label=f'{label[j]}')
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,0].set_xscale('log')
    # ax[i,0].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].yaxis.set_major_locator(MultipleLocator(ytick[i][0]))
    ax[i,0].yaxis.set_minor_locator(MultipleLocator(ytick[i][0]/2))
    ax[i,0].spines['bottom'].set_linewidth(wth)
    ax[i,0].spines['top'].set_linewidth(wth)
    ax[i,0].spines['left'].set_linewidth(wth)
    ax[i,0].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_mean[i,1:-9,idx1]),np.max(r0_mean[i,:,idx1]))
    ymin=min(np.min(r_mean[i,1:-9,idx1]),np.min(r0_mean[i,:,idx1]))
    ax[i,0].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    if i==1:
        ax[i,0].set_xlabel('$t~(\\tau)$',fontsize=size)
    ax[i,0].set_ylabel('$\\langle{}r/R_\\mathrm{m}\\rangle{}$ in %s'%axlabel[i],fontsize=size)
    if i==0:
        legend=ax[i,0].legend(bbox_to_anchor=(0.0,0.3),loc='lower left',handlelength=0.0,handletextpad=0.0,columnspacing=0.5,fontsize=0.8*size,ncol=len(lines0))
        ### Taking the lower left cornor of the legend box as the anchor point to align to bbox_to_anchor.
        for text,line in zip(legend.get_texts(),lines0):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
            text.set_color(line.get_color())
    if i==1:
        legend=ax[i,0].legend(loc='lower left',handlelength=0.0,handletextpad=0.0,columnspacing=0.5,fontsize=0.8*size,ncol=len(lines0))
        for text,line in zip(legend.get_texts(),lines0):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
            text.set_color(line.get_color())
    ax2.set_yticks([])
    ax2.set_ylim(ax[i,0].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    if i==0:
        legend=ax2.legend(loc='lower left',handlelength=0.0,handletextpad=0.0,columnspacing=0.5,fontsize=0.8*size,ncol=len(lines))
        for text,line in zip(legend.get_texts(),lines):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground='k')])
            text.set_color(line.get_color())

    ax2=ax[i,1].twinx()
    lines0=[]
    lines=[]
    for j in range(len(idx1)):
        # ax[i,1].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i,1].plot([min(t[1:-9]),t0-gap0],[r0_std[i,k,idx1[j]],r0_std[i,k,idx1[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,1].plot([min(t[1:-9]),t0-gap0],[r0_std[i,k,idx1[j]],r0_std[i,k,idx1[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
            if k==1:
                ax[i,1].plot([t0+gap0,max(t[1:-9])],[r0_std[i,k,idx1[j]],r0_std[i,k,idx1[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,1].plot([t0+gap0,max(t[1:-9])],[r0_std[i,k,idx1[j]],r0_std[i,k,idx1[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
    for k in range(ncl):
        lines0+=ax[i,1].plot([],[],color=color0[i][k],linestyle='--',linewidth=2*wth,label=label0[i][k])
    for j in range(len(idx1)):
        ax2.plot(t[1:-9],r_std[i,1:-9,idx1[j]],'-',c='k',linewidth=3*wth)
        lines+=ax2.plot(t[1:-9],r_std[i,1:-9,idx1[j]],'-',c=color[j],linewidth=2*wth,label=f'{label[j]}')
    ax[i,1].minorticks_on()
    ax[i,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,1].set_xscale('log')
    # ax[i,1].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].yaxis.set_major_locator(MultipleLocator(ytick[i][1]))
    ax[i,1].yaxis.set_minor_locator(MultipleLocator(ytick[i][1]/2))
    ax[i,1].spines['bottom'].set_linewidth(wth)
    ax[i,1].spines['top'].set_linewidth(wth)
    ax[i,1].spines['left'].set_linewidth(wth)
    ax[i,1].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_std[i,1:-9,idx1]),np.max(r0_std[i,:,idx1]))
    ymin=min(np.min(r_std[i,1:-9,idx1]),np.min(r0_std[i,:,idx1]))
    ax[i,1].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    if i==1:
        ax[i,1].set_xlabel('$t~(\\tau)$',fontsize=size)
    ax[i,1].set_ylabel('$\\langle{}\\langle{}r/R_\\mathrm{m}\\rangle{}\\rangle{}$ in %s'%axlabel[i],fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i,1].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)

    color=[(0.0,1.0,0.0),(0.0,0.9,1.0),(1.0,0.0,1.0)]
    label=['A','B','Gene']
    axlabel=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
    ax2=ax[i,2].twinx()
    lines0=[]
    lines=[]
    for j in range(len(idx2)):
        # ax[i,2].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i,2].plot([min(t[1:-9]),t0-gap0],[r0_mean[i,k,idx2[j]],r0_mean[i,k,idx2[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,2].plot([min(t[1:-9]),t0-gap0],[r0_mean[i,k,idx2[j]],r0_mean[i,k,idx2[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
            if k==1:
                ax[i,2].plot([t0+gap0,max(t[1:-9])],[r0_mean[i,k,idx2[j]],r0_mean[i,k,idx2[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,2].plot([t0+gap0,max(t[1:-9])],[r0_mean[i,k,idx2[j]],r0_mean[i,k,idx2[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
    for k in range(ncl):
        lines0+=ax[i,2].plot([],[],color=color0[i][k],linestyle='--',linewidth=2*wth,label=label0[i][k])
    for j in range(len(idx2)):
        ax2.plot(t[1:-9],r_mean[i,1:-9,idx2[j]],'-',c='k',linewidth=3*wth)
        lines+=ax2.plot(t[1:-9],r_mean[i,1:-9,idx2[j]],'-',c=color[j],linewidth=2*wth,label=f'{label[j]}')
    ax[i,2].minorticks_on()
    ax[i,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,2].set_xscale('log')
    # ax[i,2].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,2].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].yaxis.set_major_locator(MultipleLocator(ytick[i][2]))
    ax[i,2].yaxis.set_minor_locator(MultipleLocator(ytick[i][2]/2))
    ax[i,2].spines['bottom'].set_linewidth(wth)
    ax[i,2].spines['top'].set_linewidth(wth)
    ax[i,2].spines['left'].set_linewidth(wth)
    ax[i,2].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_mean[i,1:-9,idx2]),np.max(r0_mean[i,:,idx2]))
    ymin=min(np.min(r_mean[i,1:-9,idx2]),np.min(r0_mean[i,:,idx2]))
    ax[i,2].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    if i==1:
        ax[i,2].set_xlabel('$t~(\\tau)$',fontsize=size)
    ax[i,2].set_ylabel('$\\langle{}r/R_\\mathrm{g}\\rangle{}$ in %s'%axlabel[i],fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i,2].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)

    ax2=ax[i,3].twinx()
    lines0=[]
    lines=[]
    for j in range(len(idx2)):
        # ax[i,3].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
        for k in range(ncl):
            t0=np.exp((np.log(min(t[1:-9]))+np.log(max(t[1:-9])))/2)
            gap0=np.exp(0.1*(np.log(min(t[1:-9]))-np.log(max(t[1:-9]))))
            if k==0:
                ax[i,3].plot([min(t[1:-9]),t0-gap0],[r0_std[i,k,idx2[j]],r0_std[i,k,idx2[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,3].plot([min(t[1:-9]),t0-gap0],[r0_std[i,k,idx2[j]],r0_std[i,k,idx2[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
            if k==1:
                ax[i,3].plot([t0+gap0,max(t[1:-9])],[r0_std[i,k,idx2[j]],r0_std[i,k,idx2[j]]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                ax[i,3].plot([t0+gap0,max(t[1:-9])],[r0_std[i,k,idx2[j]],r0_std[i,k,idx2[j]]],color=color0[i][k],linestyle='-',linewidth=2*wth)
    for k in range(ncl):
        lines0+=ax[i,3].plot([],[],color=color0[i][k],linestyle='--',linewidth=2*wth,label=label0[i][k])
    for j in range(len(idx2)):
        ax2.plot(t[1:-9],r_std[i,1:-9,idx2[j]],'-',c='k',linewidth=3*wth)
        lines+=ax2.plot(t[1:-9],r_std[i,1:-9,idx2[j]],'-',c=color[j],linewidth=2*wth,label=f'{label[j]}')
    ax[i,3].minorticks_on()
    ax[i,3].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,3].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[i,3].set_xscale('log')
    # ax[i,3].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,3].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,3].yaxis.set_major_locator(MultipleLocator(ytick[i][3]))
    ax[i,3].yaxis.set_minor_locator(MultipleLocator(ytick[i][3]/2))
    ax[i,3].spines['bottom'].set_linewidth(wth)
    ax[i,3].spines['top'].set_linewidth(wth)
    ax[i,3].spines['left'].set_linewidth(wth)
    ax[i,3].spines['right'].set_linewidth(wth)
    ymax=max(np.max(r_std[i,1:-9,idx2]),np.max(r0_std[i,:,idx2]))
    ymin=min(np.min(r_std[i,1:-9,idx2]),np.min(r0_std[i,:,idx2]))
    ax[i,3].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    if i==1:
        ax[i,3].set_xlabel('$t~(\\tau)$',fontsize=size)
    ax[i,3].set_ylabel('$\\langle{}\\langle{}r/R_\\mathrm{g}\\rangle{}\\rangle{}$ in %s'%axlabel[i],fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[i,3].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_mean_std.png', format='png')
plt.savefig(f'{figname}_mean_std.pdf', format='pdf')
# plt.show()

