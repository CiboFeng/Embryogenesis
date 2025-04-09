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
from matplotlib import rcParams
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file1='../dist_prob_evol/dist_prob.npz'
file2='../insul_strth/insul_strth.npz'
file3='../cpmt_strth_dekker/cpmt_strth.npz'
file4='../1d_3d_dist/1d_3d_dist.npz'
figname='str'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
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
idx_cl=[0,4]
mb=0.1
scp=5

"""Compartmentalization Strength"""
data=np.load(file3)
CSs_ref=data['CSs_ref']
CSs=data['CSs']
CSs_tot_ref=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        CSs_tot_ref[i,j]=(np.abs(CSs_ref[i,idx_cl[j],0])+np.abs(CSs_ref[i,idx_cl[j],1]))+np.abs(CSs_ref[i,idx_cl[j],2])
CSs_tot=np.zeros((nsw,nfr))
for i in range(nsw):
    for j in range(nfr):
        CSs_tot[i,j]=(np.abs(CSs[i,j,0])+np.abs(CSs[i,j,1]))+np.abs(CSs[i,j,2])

"""Insulation Strength"""
data=np.load(file2)
ISs=data['ISs']
ISs_ref=data['ISs_ref']
d2=np.arange(-scp,scp+1,1)*mb

ISs_tot_ref=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        ISs_tot_ref[i,j]=np.max(ISs_ref[i,idx_cl[j]])-np.min(ISs_ref[i,idx_cl[j]])

ISs_tot=np.zeros((nsw,nfr))
for i in range(nsw):
    for j in range(nfr):
        ISs_tot[i,j]=np.max(ISs[i,j])-np.min(ISs[i,j])

"""Distance Probability"""
data=np.load(file1)
DP_ref=data['DP_ref']
DP=data['DP']
d1=np.arange(0,nat,1)*mb

DP_nor_ref=DP_ref*d1[np.newaxis,np.newaxis,:]

DP_nor=DP*d1[np.newaxis,np.newaxis,:]

slpl_ref=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        slpl_ref[i,j]=np.polyfit(np.log(d1[5:70]),np.log(DP_ref[i,idx_cl[j],5:70]),1)[0]

slpl=np.zeros((nsw,nfr))
for i in range(nsw):
    for j in range(nfr):
        slpl[i,j]=np.polyfit(np.log(d1[5:70]),np.log(DP[i,j,5:70]),1)[0]

"""Spatial versus Sequence Distance"""
s=np.arange(0,nat/10,0.1)

data=np.load(file4)
d0_mean=data['d0_mean']/sgm
d0_std=data['d0_std']/sgm
slp0_mean=np.zeros((nsw,ncl))
slp0_std=np.zeros((nsw,ncl))
for i in range(nsw):
    for j in range(ncl):
        # par=curve_fit(f,s[100:],d0_mean[i,j,100:])[0][0]
        # err=1-r2_score(f(s[100:],par),d0_mean[i,j,100:])
        # slp0_mean[i,j]=par
        # slp0_std[i,j]=err
        slp0_mean[i,j]=np.mean(d0_mean[i,j,100:]/s[np.newaxis,np.newaxis,100:]**0.2,axis=-1)
        slp0_std[i,j]=np.std(d0_mean[i,j,100:]/s[np.newaxis,np.newaxis,100:]**0.2,axis=-1)

d_mean=data['d_mean']/sgm
d_std=data['d_std']/sgm
slp_mean=np.zeros((nsw,nfr))
slp_std=np.zeros((nsw,nfr))
for i in range(nsw):
    for j in range(nfr):
        # par=curve_fit(f,s[100:],d_mean[i,j,100:])[0][0]
        # err=1-r2_score(f(s[100:],par),d_mean[i,j,100:])
        # slp_mean[i,j]=par
        # slp_std[i,j]=err
        slp_mean[i,j]=np.mean(d_mean[i,j,100:]/s[np.newaxis,np.newaxis,100:]**0.2,axis=-1)
        slp_std[i,j]=np.std(d_mean[i,j,100:]/s[np.newaxis,np.newaxis,100:]**0.2,axis=-1)

"""Merge"""
Q_ref=np.zeros((4,nsw,ncl))
Q_ref[0]=CSs_tot_ref
Q_ref[1]=ISs_tot_ref
Q_ref[2]=slpl_ref
Q_ref[3]=slp0_mean
Q=np.zeros((4,nsw,nfr))
Q[0]=CSs_tot
Q[1]=ISs_tot
Q[2]=slpl
Q[3]=slp_mean

E_ref=np.zeros((4,nsw,ncl))
E_ref[3]=slp0_std
E=np.zeros((4,nsw,nfr))
E[3]=slp_std

"""Plot"""
fig,ax=plt.subplots(3,4,figsize=(30,15))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=0.0
ytick=[[0.5,0.3,0.0,10],[0.5,0.3,0.0,10],[0.4,0.2,0.3,1]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))
ann=['(A)','(D)','(G)','(J)','(B)','(E)','(H)','(K)','(C)','(F)','(I)','(L)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

for i in range(nsw):
    lines0=[]
    for j in range(nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        for k in range(3):
            ax[i,0].plot([k-0.4,k+0.4],[CSs[i,j,k],CSs[i,j,k]],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        for k in range(3):
            ax[i,0].plot([k-0.4,k+0.4],[CSs_ref[i,idx_cl[j],k],CSs_ref[i,idx_cl[j],k]],color=(0.5,0.5,0.5),linestyle='-',linewidth=2*wth)
            ax[i,0].plot([k-0.4,k+0.4],[CSs_ref[i,idx_cl[j],k],CSs_ref[i,idx_cl[j],k]],color=color0[i][j],linestyle='-',linewidth=wth)
        lines0+=ax[i,0].plot([],[],'-.',color=color0[i][j],linestyle='-',linewidth=wth,label=f'{label0[i][j]}')
    ax[i,0].minorticks_on()
    ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,0].tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax[i,0].tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    xticks_positions=[0,1,2]
    xticks_labels=['AA','BB','AB']
    ax[i,0].set_xticks(xticks_positions,xticks_labels)
    # ax[i,0].xaxis.set_major_locator(MultipleLocator(xtick))
    # ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].yaxis.set_major_locator(MultipleLocator(ytick[i][0]))
    ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,0].spines['bottom'].set_linewidth(wth)
    ax[i,0].spines['top'].set_linewidth(wth)
    ax[i,0].spines['left'].set_linewidth(wth)
    ax[i,0].spines['right'].set_linewidth(wth)
    ax[i,0].set_ylim(-0.8,0.43)
    ax[i,0].set_ylabel('$\\overline{S_{\\mathrm{c}}}$ in %s'%label[i],fontsize=size)
    ax[i,0].set_xlabel('Type',fontsize=size)
    legend=ax[i,0].legend(loc='upper right',handlelength=0,handletextpad=0,fontsize=0.8*size,)
    for text,line in zip(legend.get_texts(),lines0):
        text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
        text.set_color(line.get_color())

    lines0=[]
    for j in range(nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,1].plot(d2,ISs[i,j,:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i,1].plot(d2,ISs_ref[i,idx_cl[j],:],color=(0.5,0.5,0.5),linestyle='-',linewidth=2*wth)
        lines0+=ax[i,1].plot(d2,ISs_ref[i,idx_cl[j],:],color=color0[i][j],linestyle='-',linewidth=wth,label=f'{label0[i][j]}')
    ax[i,1].minorticks_on()
    ax[i,1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,1].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].yaxis.set_major_locator(MultipleLocator(ytick[i][1]))
    ax[i,1].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,1].spines['bottom'].set_linewidth(wth)
    ax[i,1].spines['top'].set_linewidth(wth)
    ax[i,1].spines['left'].set_linewidth(wth)
    ax[i,1].spines['right'].set_linewidth(wth)
    ax[i,1].set_ylim(-0.42,0.25)
    ax[i,1].set_xlabel('$\\xi_\\mathrm{b}$ (Mb)',fontsize=size)
    ax[i,1].set_ylabel('$\\overline{S_{\\mathrm{i}}}$ in %s'%label[i],fontsize=size)

    lines0=[]
    for j in range(nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,2].plot(d1[1:],DP_nor[i,j,1:],'-',color=rgb,linewidth=wth/4)
    for j in range(ncl):
        ax[i,2].plot(d1[1:],DP_nor_ref[i,idx_cl[j],1:],color=(0.5,0.5,0.5),linestyle='-',linewidth=2*wth)
        ax[i,2].plot(d1[1:],DP_nor_ref[i,idx_cl[j],1:],color=color0[i][j],linestyle='-',linewidth=wth)
        lines0+=ax[i,2].plot([],[],'-.',color=color0[i][j],linestyle='--',linewidth=wth,label=f'{label0[i][j]}')
    ax[i,2].minorticks_on()
    ax[i,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
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
    ax[i,2].set_ylim(0.005,0.2)
    ax[i,2].set_xlabel('$\\xi$ (Mb)',fontsize=size)
    ax[i,2].set_ylabel('$p/\\xi^{-1}$ in %s'%label[i],fontsize=size)

    lines0=[]
    for j in range(1,nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,3].plot(s,d_mean[i,j,:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i,3].plot(s,d0_mean[i,j,:],'-',color=(0.5,0.5,0.5),linewidth=2*wth)
        lines0+=ax[i,3].plot(s,d0_mean[i,j,:],'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
    ax[i,3].minorticks_on()
    ax[i,3].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,3].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,3].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,3].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,3].yaxis.set_major_locator(MultipleLocator(ytick[i][3]))
    ax[i,3].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,3].spines['bottom'].set_linewidth(wth)
    ax[i,3].spines['top'].set_linewidth(wth)
    ax[i,3].spines['left'].set_linewidth(wth)
    ax[i,3].spines['right'].set_linewidth(wth)
    ax[i,3].set_xlabel('$\\xi$ (Mb)',fontsize=size)
    ax[i,3].set_ylabel('$d~(\\sigma)$\nin %s'%label[i],fontsize=size)

color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0.0,1.0,1.0),(1.0,0.3,0.3)]
ylabel=['$A_{\\mathrm{c}}$','$A_{\\mathrm{i}}$','$\\gamma$','$\\omega~(\\sigma)$']

for i in range(4):
    ax2=ax[2,i].twinx()
    lines0=[]
    lines=[]
    for j in range(nsw):
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[2,i].plot([min(t[1:-9]),max(t[1:-9])],[Q_ref[i,j,k],Q_ref[i,j,k]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                lines0+=ax[2,i].plot([min(t[1:-9]),max(t[1:-9])],[Q_ref[i,j,k],Q_ref[i,j,k]],color=color0[j][k],linestyle='-',linewidth=2*wth,label=f'{label0[j][k]}')
        ax2.plot(t[1:-9],Q[i,j,1:-9],'-',color='k',linewidth=3*wth)
        lines+=ax2.plot(t[1:-9],Q[i,j,1:-9],'-',color=color[j],linewidth=2*wth,label=f'{label[j]}')
    ax[2,i].minorticks_on()
    ax[2,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[2,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[2,i].xaxis.set_major_locator(MultipleLocator(xtick))
    # ax[2,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[2,i].yaxis.set_major_locator(MultipleLocator(ytick[2][i]))
    ax[2,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[2,i].spines['bottom'].set_linewidth(wth)
    ax[2,i].spines['top'].set_linewidth(wth)
    ax[2,i].spines['left'].set_linewidth(wth)
    ax[2,i].spines['right'].set_linewidth(wth)
    ymax=max(np.max(Q[i,:,1:-9]),np.max(Q_ref[i]))
    ymin=min(np.min(Q[i,:,1:-9]),np.min(Q_ref[i]))
    ax[2,i].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    ax[2,i].set_xscale('log')
    ax[2,i].set_ylabel(ylabel[i],fontsize=size)
    ax[2,i].set_xlabel(r'$t~(\tau)$',fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[2,i].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    if i==0:
        legend=ax2.legend(bbox_to_anchor=(1.0,0.06),loc='lower right',handlelength=0.0,handletextpad=0.0,fontsize=0.8*size)
        for text,line in zip(legend.get_texts(),lines):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground='k')])
            text.set_color(line.get_color())
        # legend.set_bbox_to_anchor((0.8,0.75))

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label(r'$t~(\tau)$',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)
plt.tight_layout(rect=[0,0,0.85,1])

plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()