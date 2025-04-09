"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
from insul_score import insul_score

"""Set Aruments"""
reffile='../cont_prob/ref/esc.dat'
file1='../tad_dist/tad_dist.pyw'
file2='../1d_3d_dist/1d_3d_dist.npz'
figname='scl_law'
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

"""Define TAD"""
P_ref=np.zeros((nat,nat))
f=open(reffile,'r')
lsf=f.readlines()
for k in range(len(lsf)):
    ls_f=lsf[k].strip('\n').split()
    P_ref[int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])

pos=insul_score(P_ref,bin_size=10**5)[1]
k=1
while k<len(pos):
    if pos[k]-pos[k-1]<=1:
        pos.pop(k)
        k-=1
    k+=1
ntad=len(pos)-1

"""Read Data"""
data=np.load(file2)
d0_mean=data['d0_mean']/sgm
d0_std=data['d0_std']/sgm
d_mean=data['d_mean']/sgm
d_std=data['d_std']/sgm

l=np.arange(0,nat/10,0.1)

D0_mean=np.zeros((nsw,ncl,ntad))
D0_std=np.zeros((nsw,ncl,ntad))
for i in range(nsw):
    q=pyw(file1,f'TAD (Reference Value, {sw[i]}):')
    for j in range(ncl):
        for k in range(ntad):
            D0_mean[i,j,k]=q[2][j*ntad+k]/sgm
            D0_std[i,j,k]=q[3][j*ntad+k]/sgm

D_mean=np.zeros((nsw,nfr,ntad))
D_std=np.zeros((nsw,nfr,ntad))
for i in range(nsw):
    q=pyw(file1,f'TAD ({sw[i]}):')
    for j in range(nfr):
        for k in range(ntad):
            D_mean[i,j,k]=q[2][j*ntad+k]/sgm
            D_std[i,j,k]=q[3][j*ntad+k]/sgm

L=np.arange(0,ntad,1)

"""Plot"""
fig=plt.figure(figsize=(25,10))
gs=gridspec.GridSpec(2,15,width_ratios=([1,0.05,0.1,1.2]*4)[:-1],wspace=0,hspace=0.3)
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
ann=[['(A)','(A)','(A)','(A)','(C)','(C)','(C)','(C)','(E)','(E)','(E)','(E)','(G)','(G)','(G)'],
     ['(B)','(B)','(B)','(B)','(D)','(D)','(D)','(D)','(F)','(F)','(F)','(F)','(H)','(H)','(H)']]

cmap=[]
norm=[]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
norm+=[colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))]
color=[(1,0.8,0),(1,0,0),(0,0,0)]
nodes=[0/9,4/9,9/9]
cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
norm+=[colors.Normalize(vmin=1,vmax=nat/10)]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
norm+=[colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))]
color=[(1,0.8,0),(1,0,0),(0,0,0)]
nodes=[0/9,4/9,9/9]
cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
norm+=[colors.Normalize(vmin=1,vmax=ntad)]

for i in range(nsw):
    ax=fig.add_subplot(gs[i,0])
    ax.annotate(ann[i][0],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    lines0=[]
    for j in range(1,nfr-9):
        normc=norm[2](t[j])
        rgb=cmap[2](normc)
        ax.plot(l[1:],d_mean[i,j,1:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax.plot(l[1:],d0_mean[i,j,1:],'-',color=(0.5,0.5,0.5),linewidth=2*wth)
        lines0+=ax.plot(l[1:],d0_mean[i,j,1:],'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_major_locator(MultipleLocator(xtick))
    # ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(ytick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if i==1:
        ax.set_xlabel('$\\xi$ (Mb)',fontsize=size)
    ax.set_ylabel('$d~(\\sigma)$\nin %s'%label[i],fontsize=size)
    legend=ax.legend(loc='lower right',handlelength=0.0,handletextpad=0.0,fontsize=0.8*size)
    for text,line in zip(legend.get_texts(),lines0):
        text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
        text.set_color(line.get_color())

    ax=fig.add_subplot(gs[i,4])
    ax.annotate(ann[i][4],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    for j in [0,2,4,8,16,32,64,128,256,512]:
        normc=norm[3](j/10)
        rgb=cmap[3](normc)
        ax.plot(t[1:-9],d_mean[i,1:-9,j]-d0_mean[i,-1,j],'-',color=rgb,linewidth=wth)
        # ax.plot(r0_mean[i,k,:],'-',color='k',linewidth=2*wth)
        # ax.plot(r0_mean[i,k,:],'-',color=color0[k][i],linewidth=wth,label=label0[k][i])
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_major_locator(MultipleLocator(xtick))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(ytick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.legend(loc='best',fontsize=0.8*size)
    if i==1:
        ax.set_xlabel('$t~(\\tau)$',fontsize=size)
    ax.set_ylabel('$d-d^\\mathrm{ESC}~(\\sigma)$\nin %s'%label[i],fontsize=size)

    ax=fig.add_subplot(gs[i,8])
    ax.annotate(ann[i][8],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    lines0=[]
    for j in range(1,nfr-9):
        normc=norm[0](t[j])
        rgb=cmap[0](normc)
        ax.plot(L[1:],D_mean[i,j,1:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax.plot(L[1:],D0_mean[i,j,1:],'-',color=(0.5,0.5,0.5),linewidth=2*wth)
        lines0+=ax.plot(L[1:],D0_mean[i,j,1:],'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_major_locator(MultipleLocator(xtick))
    # ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(ytick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.set_yscale('log')
    if i==1:
        ax.set_xlabel('$\\xi^{\\mathrm{TAD}}$',fontsize=size)
    ax.set_ylabel('$d^\\mathrm{TAD}~(\\sigma)$\nin %s'%label[i],fontsize=size)
    legend=ax.legend(loc='lower right',handlelength=0.0,handletextpad=0.0,fontsize=0.8*size)
    for text,line in zip(legend.get_texts(),lines0):
        text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
        text.set_color(line.get_color())

    ax=fig.add_subplot(gs[i,12])
    ax.annotate(ann[i][12],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    for j in [0,2,4,8,16,32,64]:
        normc=norm[1](j)
        rgb=cmap[1](normc)
        ax.plot(t[1:-9],D_mean[i,1:-9,j]-D0_mean[i,-1,j],'-',color=rgb,linewidth=wth)
        # ax.plot(r0_mean[i,k,:],'-',color='k',linewidth=2*wth)
        # ax.plot(r0_mean[i,k,:],'-',color=color0[k][i],linewidth=wth,label=label0[k][i])
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_major_locator(MultipleLocator(xtick))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(ytick))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    # ax.legend(loc='best',fontsize=0.8*size)
    if i==1:
        ax.set_xlabel('$t~(\\tau)$',fontsize=size)
    ax.set_ylabel('$d^\\mathrm{TAD}-d^\\mathrm{TAD,ESC}~(\\sigma)$\nin %s'%label[i],fontsize=size)

cbar_ax=fig.add_subplot(gs[:,2])
sm=plt.cm.ScalarMappable(cmap=cmap[2],norm=norm[2])
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.ax.set_title('$t~(\\tau)$',fontsize=size,pad=20)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_subplot(gs[:,6])
sm=plt.cm.ScalarMappable(cmap=cmap[3],norm=norm[3])
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.ax.set_title('$\\xi$ (Mb)',fontsize=size,pad=20)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_subplot(gs[:,10])
sm=plt.cm.ScalarMappable(cmap=cmap[0],norm=norm[0])
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.ax.set_title('$t~(\\tau)$',fontsize=size,pad=20)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

cbar_ax=fig.add_subplot(gs[:,14])
sm=plt.cm.ScalarMappable(cmap=cmap[1],norm=norm[1])
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.ax.set_title('$\\xi^{\\mathrm{TAD}}$',fontsize=size,pad=20)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0,1,1])
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()