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
file='cpmt_size.pyw'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nhist=20
sgm=4.0
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read Data"""
r0=np.zeros((nsw,ncl,2,nhist))
P_r0=np.zeros((nsw,ncl,2,nhist))
r0_mean=np.zeros((nsw,ncl,2))
r0_std=np.zeros((nsw,ncl,2))
for i in range(nsw):
    q=pyw(file,f'Probability Density (Reference Value, {sw[i]}):')
    for j in range(ncl):
        for k in range(nhist):
            r0[i,j,0,k]=q[1][j*nhist+k]/sgm
            P_r0[i,j,0,k]=q[2][j*nhist+k]*sgm
            r0[i,j,1,k]=q[3][j*nhist+k]/sgm
            P_r0[i,j,1,k]=q[4][j*nhist+k]*sgm
    q=pyw(file,f'and B (Reference Value, {sw[i]}):')
    for j in range(ncl):
        r0_mean[i,j,0]=q[1][j]/sgm
        r0_std[i,j,0]=q[2][j]/sgm
        r0_mean[i,j,1]=q[3][j]/sgm
        r0_std[i,j,1]=q[4][j]/sgm

r=np.zeros((nsw,nfr,2,nhist))
P_r=np.zeros((nsw,nfr,2,nhist))
r_mean=np.zeros((nsw,nfr,2))
r_std=np.zeros((nsw,nfr,2))
for i in range(nsw):
    q=pyw(file,f'Probability Density ({sw[i]}):')
    for j in range(nfr):
        for k in range(nhist):
            r[i,j,0,k]=q[1][j*nhist+k]/sgm
            P_r[i,j,0,k]=q[2][j*nhist+k]*sgm
            r[i,j,1,k]=q[3][j*nhist+k]/sgm
            P_r[i,j,1,k]=q[4][j*nhist+k]*sgm
    q=pyw(file,f'and B ({sw[i]}):')
    for j in range(nfr):
        r_mean[i,j,0]=q[1][j]/sgm
        r_std[i,j,0]=q[2][j]/sgm
        r_mean[i,j,1]=q[3][j]/sgm
        r_std[i,j,1]=q[4][j]/sgm

"""Plot"""
fig,ax=plt.subplots(2,2,figsize=(20,10))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
axlabel=['A','B']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))

for i in range(nsw):
    for j in range(2):
        for k in range(1,nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i,j].plot(r[i,k,j,:],P_r[i,k,j,:],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[i,j].plot(r0[i,k,j,:],P_r0[i,k,j,:],'-',color='k',linewidth=2*wth)
            ax[i,j].plot(r0[i,k,j,:],P_r0[i,k,j,:],'-',color=color0[i][k],linewidth=wth,label=label0[i][k])
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
        ax[i,j].legend(loc='best',fontsize=0.8*size)
        if j==0:
            ax[i,j].set_ylabel(f'Distribution\nin {sw[i]}',fontsize=size)
        if i==1:
            ax[i,j].set_xlabel(f'Radius of Gyration of {axlabel[j]} $(\\sigma)$',fontsize=size)

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
color=[[(0,0.5,0),(0.5,0,0)],[(0,1,0),(1,0,0)]]

ax2=ax.twinx()
for i in range(nsw):
    for j in range(2):
        # ax.fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
        # for k in range(ncl):
        #     if (i,k)!=(0,1):
        #         ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
        ax2.plot(t[1:-9],r_mean[i,1:-9,j],'-',color=color[i][j],linewidth=wth,label=axlabel[j]+' in '+sw[i])
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
# ymax=max(np.max(r_mean[:,1:-9,:]),np.max(r0_mean))
# ymin=min(np.min(r_mean[:,1:-9,:]),np.min(r0_mean))
# ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
# ax.legend(loc='lower right',fontsize=0.8*size)
ax.set_ylim(ax2.get_ylim())
ax.set_xlabel(f'Time $(\\tau)$',fontsize=size)
ax.set_ylabel('$\\langle R_{\\mathrm{g}} \\rangle ~(\\sigma)$',fontsize=size)
ax2.set_yticks([])
# ax2.set_ylim(ax.get_ylim())
ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax2.legend(loc='upper right',fontsize=0.6*size)

plt.tight_layout()
plt.savefig(f'{figname}_mean.png', format='png')
plt.savefig(f'{figname}_mean.pdf', format='pdf')
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
color=[[(0,0.5,0),(0.5,0,0)],[(0,1,0),(1,0,0)]]

ax2=ax.twinx()
for i in range(nsw):
    for j in range(2):
        # ax.fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
        # for k in range(ncl):
        #     if (i,k)!=(0,1):
        #         ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
        ax2.plot(t[1:-9],(r_mean[i,1:-9,j]-r_mean[i,-10,j])/(r_mean[i,1,j]-r_mean[i,-10,j]),'-',color=color[i][j],linewidth=wth,label=axlabel[j]+' in '+sw[i])
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
# ymax=max(np.max(r_mean[:,1:-9,:]),np.max(r0_mean))
# ymin=min(np.min(r_mean[:,1:-9,:]),np.min(r0_mean))
# ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
# ax.legend(loc='lower right',fontsize=0.8*size)
ax.set_ylim(ax2.get_ylim())
ax.set_xlabel(f'Time $(\\tau)$',fontsize=size)
ax.set_ylabel('$\\frac{\\langle{}R_\\mathrm{g}\\rangle-\\langle{}R_\\mathrm{g}^\\mathrm{end}\\rangle}'
              '{\\langle{}R_\\mathrm{g}^\\mathrm{start}\\rangle-\\langle{}R_\\mathrm{g}^\\mathrm{end}\\rangle}$',fontsize=size)
ax2.set_yticks([])
# ax2.set_ylim(ax.get_ylim())
ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax2.legend(loc='upper right',fontsize=0.6*size)

plt.tight_layout()
plt.savefig(f'{figname}_mean1.png', format='png')
plt.savefig(f'{figname}_mean1.pdf', format='pdf')
plt.show()