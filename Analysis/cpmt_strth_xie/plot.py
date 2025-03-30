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

"""Set Arguments"""
file='cpmt_strth.npz'
figname=file[:-4]
cl=['Z','8C','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
idx_cl=[0,3,4]

"""Read Data and Process"""
data=np.load(file)
CSs_ref=data['CSs_ref']
CSs=data['CSs']
x=np.arange(0,3,1)
ISs_tot_ref=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        ISs_tot_ref[i,j]=(np.abs(CSs_ref[i,idx_cl[j],0])+np.abs(CSs_ref[i,idx_cl[j],1]))/np.abs(CSs_ref[i,idx_cl[j],2])
ISs_tot=np.zeros((nsw,nfr))
for i in range(nsw):
    for j in range(nfr):
        ISs_tot[i,j]=(np.abs(CSs[i,j,0])+np.abs(CSs[i,j,1]))/np.abs(CSs[i,j,2])

"""Plot"""
fig,ax=plt.subplots(2,1,figsize=(10,8))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[(0.8,0.8,0.8),(0.4,0.4,0.4),(0.0,0.0,0.0)]
linestyle0=[['-.','--'],['-.','--'],[':',':']]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))

for i in range(nsw):
    for j in range(nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i].plot(x,CSs[i,j,:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i].plot(x,CSs_ref[i,idx_cl[j],:],'-.',color=color0[j],linestyle=linestyle0[j][i],linewidth=wth,label=f'{cl[j]}')
    ax[i].minorticks_on()
    ax[i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i].tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax[i].tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    xticks_positions=[0,1,2]
    xticks_labels=['AA','BB','AB']
    ax[i].set_xticks(xticks_positions,xticks_labels)
    # ax[i].xaxis.set_major_locator(MultipleLocator(xtick))
    # ax[i].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i].spines['bottom'].set_linewidth(wth)
    ax[i].spines['top'].set_linewidth(wth)
    ax[i].spines['left'].set_linewidth(wth)
    ax[i].spines['right'].set_linewidth(wth)
    ax[i].set_ylabel(f'Cpmt. Strth.\nof {sw[i]}',fontsize=size)
# ax[1].set_xlabel('Dist. to Bd. (Mb)',fontsize=size)

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('Time $(\\tau)$',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)
plt.tight_layout(rect=[0,0,0.85,1])

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()

fig,ax=plt.subplots(2,1,figsize=(9,8))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[(0.8,0.8,0.8),(0.4,0.4,0.4),(0.0,0.0,0.0)]
linestyle0=[['-.','--'],['-.','--'],[':',':']]

for i in range(nsw):
    ax[i].plot(t[:-9],ISs_tot[i,:-9],'-',color='r',linewidth=wth)
    for j in range(ncl):
        ax[i].plot([min(t[:-9]),max(t[:-9])],[ISs_tot_ref[i,j],ISs_tot_ref[i,j]],color=color0[j],linestyle=linestyle0[j][i],linewidth=wth,label=f'{cl[j]}')
    ax[i].minorticks_on()
    ax[i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i].xaxis.set_major_locator(MultipleLocator(xtick))
    # ax[i].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[i].yaxis.set_major_locator(MultipleLocator(ytick))
    ax[i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i].spines['bottom'].set_linewidth(wth)
    ax[i].spines['top'].set_linewidth(wth)
    ax[i].spines['left'].set_linewidth(wth)
    ax[i].spines['right'].set_linewidth(wth)
    ax[i].set_xscale('log')
    ax[i].set_ylabel('Cpmt. Strth.\nof %s'%sw[i],fontsize=size)
    ax[i].legend(loc='lower right',fontsize=0.8*size)
ax[1].set_xlabel('Time ($\\tau$)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_tot.png',format='png')
plt.savefig(f'{figname}_tot.pdf',format='pdf')
plt.show()
