"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
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
from insul_score import insul_score

"""Set Aruments"""
reffile='../cont_prob/ref/esc.dat'
file1='../tad_size/tad_size.pyw'
file2='../tad_dist/tad_nbor_dist.pyw'
file3='../tad_dist/tad_dist.pyw'
figname='tad'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nat=600
nhist=50
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
S=np.arange(0,ntad,1)

"""Define the Fitting Function"""
def F(x,a):
    return a*np.log(x)

def f(x,a):
    return a*x**0.2

"""Read Data"""
R0=np.zeros((nsw,ncl,nhist,2))
P_R0=np.zeros((nsw,ncl,nhist,2))
R0_mean=np.zeros((nsw,ncl,2))
R0_std=np.zeros((nsw,ncl,2))
D0_mean=np.zeros((nsw,ncl,ntad))
D0_std=np.zeros((nsw,ncl,ntad))
Slp0_mean=np.zeros((nsw,ncl))
Slp0_std=np.zeros((nsw,ncl))
for i in range(nsw):
    q=pyw(file1,f'Probability Density (Reference Value, {sw[i]}):')
    for j in range(ncl):
        for k in range(nhist):
            R0[i,j,k,0]=q[1][j*nhist+k]/sgm
            P_R0[i,j,k,0]=q[2][j*nhist+k]*sgm
    q=pyw(file1,f'TAD (Reference Value, {sw[i]}):')
    for j in range(ncl):
        R0_mean[i,j,0]=q[1][j]/sgm
        R0_std[i,j,0]=q[2][j]/sgm
    q=pyw(file2,f'Probability Density (Reference Value, {sw[i]}):')
    for j in range(ncl):
        for k in range(nhist):
            R0[i,j,k,1]=q[1][j*nhist+k]/sgm
            P_R0[i,j,k,1]=q[2][j*nhist+k]*sgm
    q=pyw(file2,f'TAD (Reference Value, {sw[i]}):')
    for j in range(ncl):
        R0_mean[i,j,1]=q[1][j]/sgm
        R0_std[i,j,1]=q[2][j]/sgm
    q=pyw(file3,f'TAD (Reference Value, {sw[i]}):')
    for j in range(ncl):
        for k in range(ntad):
            D0_mean[i,j,k]=q[2][j*ntad+k]/sgm
            D0_std[i,j,k]=q[3][j*ntad+k]/sgm

for i in range(nsw):
    for j in range(ncl):
        # # par=np.polyfit(np.log(s[10:]),np.log(d0_mean[i,j,10:]),1)
        # # err=1-r2_score(par[0]*np.log(s[10:])+par[1],np.log(d0_mean[i,j,10:]))
        # par=curve_fit(F,S[10:],D0_mean[i,j,10:])[0][0]
        # err=1-r2_score(F(S[10:],par),D0_mean[i,j,10:])
        # # par=np.mean(d0_mean[i,j,10:]-d0_mean[i,-1,10:])
        # # err=0
        # Slp0_mean[i,j]=par
        # Slp0_std[i,j]=err
        Slp0_mean[i,j]=np.mean(D0_mean[i,j,10:]/np.log(S[np.newaxis,np.newaxis,10:]),axis=-1)
        Slp0_std[i,j]=np.std(D0_mean[i,j,10:]/np.log(S[np.newaxis,np.newaxis,10:]),axis=-1)

R=np.zeros((nsw,nfr,nhist,2))
P_R=np.zeros((nsw,nfr,nhist,2))
R_mean=np.zeros((nsw,nfr,2))
R_std=np.zeros((nsw,nfr,2))
D_mean=np.zeros((nsw,nfr,ntad))
D_std=np.zeros((nsw,nfr,ntad))
Slp_mean=np.zeros((nsw,nfr))
Slp_std=np.zeros((nsw,nfr))
for i in range(nsw):
    q=pyw(file1,f'Probability Density ({sw[i]}):')
    for j in range(nfr):
        for k in range(nhist):
            R[i,j,k,0]=q[1][j*nhist+k]/sgm
            P_R[i,j,k,0]=q[2][j*nhist+k]*sgm
    q=pyw(file1,f'TAD ({sw[i]}):')
    for j in range(nfr):
        R_mean[i,j,0]=q[1][j]/sgm
        R_std[i,j,0]=q[2][j]/sgm
    q=pyw(file2,f'Probability Density ({sw[i]}):')
    for j in range(nfr):
        for k in range(nhist):
            R[i,j,k,1]=q[1][j*nhist+k]/sgm
            P_R[i,j,k,1]=q[2][j*nhist+k]*sgm
    q=pyw(file2,f'TAD ({sw[i]}):')
    for j in range(nfr):
        R_mean[i,j,1]=q[1][j]/sgm
        R_std[i,j,1]=q[2][j]/sgm
    q=pyw(file3,f'TAD ({sw[i]}):')
    for j in range(nfr):
        for k in range(ntad):
            D_mean[i,j,k]=q[2][j*ntad+k]/sgm
            D_std[i,j,k]=q[3][j*ntad+k]/sgm

for i in range(nsw):
    for j in range(nfr):
        # # par=np.polyfit(np.log(s[10:]),np.log(d_mean[i,j,10:]),1)
        # # err=1-r2_score(par[0]*np.log(s[10:])+par[1],np.log(d_mean[i,j,10:]))
        # par=curve_fit(F,S[10:],D_mean[i,j,10:])[0][0]
        # err=1-r2_score(F(S[10:],par),D_mean[i,j,10:])
        # # par=np.mean(d_mean[i,j,10:]-d_mean[i,-1,10:])
        # # err=0
        # Slp_mean[i,j]=par
        # Slp_std[i,j]=err
        Slp_mean[i,j]=np.mean(D_mean[i,j,10:]/np.log(S[np.newaxis,np.newaxis,10:]),axis=-1)
        Slp_std[i,j]=np.std(D_mean[i,j,10:]/np.log(S[np.newaxis,np.newaxis,10:]),axis=-1)




# par=np.zeros((nsw,nfr,2))
# for i in range(nsw):
#     for j in range(nfr):
#         par[i,j]=np.polyfit(np.log(s[10:]),np.log(d_mean[i,j,10:]),1)
#         # err=1-r2_score(par[0]*np.log(s[5:])+par[1],np.log(d_mean[i,j,5:]))
# fig,ax=plt.subplots(2,1,figsize=(8,8))
# wth=3
# size=30
# lenmaj=15
# lenmin=8
# xtick=1
# ytick=0.02
# color0=[[(1.0,0.8,0.6),(0.8,0.6,1.0)],[(0.8,1.0,0.6),(0.8,1.0,0.6)]]
# label0=[['ZP','ZM'],['ESC','ESC']]
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))
# label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
#
# for i in range(nsw):
#     # for j in range(1,nfr-9):
#     for j in [0,10,20,30,40,50]:
#         normc=norm(t[j])
#         rgb=cmap(normc)
#         ax[i].plot(np.log(S[10:]),D_mean[i,j,10:],'-',color=rgb,linewidth=wth)
#         ax[i].plot(np.log(S[10:]),Slp_mean[i,j]*np.log(S[10:]),'-',color='k',linewidth=wth)
#         ax[i].plot(np.log(S[10:]),Slp_mean[i,j]*np.log(S[10:]),'-',color=rgb,linewidth=wth/2)
#         # ax[i].plot(S[10:],D_mean[i,j,10:],'-',color=rgb,linewidth=wth)
#         # ax[i].plot(S[10:],Slp_mean[i,j]*np.log(S[10:]),'-',color='k',linewidth=wth)
#         # ax[i].plot(S[10:],Slp_mean[i,j]*np.log(S[10:]),'-',color=rgb,linewidth=wth/2)
#         # ax[i].plot(np.log(s[100:]),np.log(d_mean[i,j,100:]),'-',color=rgb,linewidth=wth)
#         # ax[i].plot(np.log(s[100:]),0.2*np.log(s[100:])+np.log(slp_mean[i,j]),'-',color='k',linewidth=wth)
#         # ax[i].plot(np.log(s[100:]),0.2*np.log(s[100:])+np.log(slp_mean[i,j]),'-',color=rgb,linewidth=wth/2)
#         # ax[i].plot(s[100:],d_mean[i,j,100:],'-',color=rgb,linewidth=wth)
#         # ax[i].plot(s[100:],slp_mean[i,j]*s[100:]**0.2,'-',color='k',linewidth=wth)
#         # ax[i].plot(s[100:],slp_mean[i,j]*s[100:]**0.2,'-',color=rgb,linewidth=wth/2)
#     ax[i].minorticks_on()
#     ax[i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax[i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     # ax[i].xaxis.set_major_locator(MultipleLocator(xtick))
#     ax[i].xaxis.set_minor_locator(AutoMinorLocator(2))
#     # ax[i].yaxis.set_major_locator(MultipleLocator(ytick))
#     ax[i].yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i].spines['bottom'].set_linewidth(wth)
#     ax[i].spines['top'].set_linewidth(wth)
#     ax[i].spines['left'].set_linewidth(wth)
#     ax[i].spines['right'].set_linewidth(wth)
# plt.show()



"""Plot"""
fig,ax=plt.subplots(3,3,figsize=(25,15))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=[[1.0,0.3,10,],[1.0,0.3,10],[0.1,0.1,0.5]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
ann=['(A)','(D)','(G)','(B)','(E)','(H)','(C)','(F)','(I)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

for i in range(nsw):
    for j in range(2):
        lines0=[]
        for k in range(1,nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[i,j].plot(R[i,k,:,j],P_R[i,k,:,j],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[i,j].plot(R0[i,k,:,j],P_R0[i,k,:,j],'-',color=(0.5,0.5,0.5),linewidth=2*wth)
            lines0+=ax[i,j].plot(R0[i,k,:,j],P_R0[i,k,:,j],'-',color=color0[i][k],linewidth=wth,label=label0[i][k])
        ax[i,j].minorticks_on()
        ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax[i,j].xaxis.set_major_locator(MultipleLocator(xtick))
        ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].yaxis.set_major_locator(MultipleLocator(ytick[i][j]))
        ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[i,j].spines['bottom'].set_linewidth(wth)
        ax[i,j].spines['top'].set_linewidth(wth)
        ax[i,j].spines['left'].set_linewidth(wth)
        ax[i,j].spines['right'].set_linewidth(wth)
        legend=ax[i,0].legend(loc='upper right',handlelength=0.0,handletextpad=0.0,fontsize=0.8*size)
        for text,line in zip(legend.get_texts(),lines0):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
            text.set_color(line.get_color())
    ax[i,0].set_xlabel('$R^{\\mathrm{TAD}}_{\\mathrm{g}}~(\\sigma)$',fontsize=size)
    ax[i,1].set_xlabel('$d^{\\mathrm{TAD}}_{\\mathrm{n}}~(\\sigma)$',fontsize=size)
    ax[i,0].set_ylabel('$f(R^{\\mathrm{TAD}}_{\\mathrm{g}})~(\\sigma^{-1})$\nin %s'%label[i],fontsize=size)
    ax[i,1].set_ylabel('$f(d^{\\mathrm{TAD}}_{\\mathrm{n}})~(\\sigma^{-1})$\nin %s'%label[i],fontsize=size)

color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0.0,1.0,1.0),(1.0,0.3,0.3)]
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']

for i in range(2):
    ax2=ax[2,i].twinx()
    lines0=[]
    lines=[]
    for j in range(nsw):
        ax[2,i].fill_between(t[1:-9],R_mean[j,1:-9,i]-R_std[j,1:-9,i]/10,R_mean[j,1:-9,i]+R_std[j,1:-9,i]/10,color=color[j],alpha=0.1)
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[2,i].plot([min(t[1:-9]),max(t[1:-9])],[R0_mean[j,k,i],R0_mean[j,k,i]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
                lines0+=ax[2,i].plot([min(t[1:-9]),max(t[1:-9])],[R0_mean[j,k,i],R0_mean[j,k,i]],color=color0[j][k],linestyle='-',linewidth=2*wth,label=label0[j][k])
        ax2.plot(t[1:-9],R_mean[j,1:-9,i],color='k',linestyle='-',linewidth=3*wth)
        lines+=ax2.plot(t[1:-9],R_mean[j,1:-9,i],color=color[j],linestyle='-',linewidth=2*wth,label=label[j])
        # ax[i].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
    ax[2,i].minorticks_on()
    ax[2,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[2,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[2,i].set_xscale('log')
    # ax[2,i].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[2,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[2,i].yaxis.set_major_locator(MultipleLocator(ytick[2][i]))
    ax[2,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[2,i].spines['bottom'].set_linewidth(wth)
    ax[2,i].spines['top'].set_linewidth(wth)
    ax[2,i].spines['left'].set_linewidth(wth)
    ax[2,i].spines['right'].set_linewidth(wth)
    ymax=max(np.max(R_mean[:,1:-9,i]),np.max(R0_mean[:,:,i]))
    ymin=min(np.min(R_mean[:,1:-9,i]),np.min(R0_mean[:,:,i]))
    ax[2,i].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    ax[2,i].set_xlabel(f'$t~(\\tau)$',fontsize=size)
    ax2.set_yticks([])
    ax2.set_ylim(ax[2,i].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    if i==0:
        legend=ax2.legend(loc='center right',handlelength=0.0,handletextpad=0.0,columnspacing=0.5,fontsize=0.8*size,ncol=1)
        for text,line in zip(legend.get_texts(),lines):
            text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground='k')])
            text.set_color(line.get_color())
ax[2,0].set_ylabel('$\\langle{}R^{\\mathrm{TAD}}_{\\mathrm{g}}\\rangle{}~(\\sigma)$',fontsize=size)
ax[2,1].set_ylabel('$\\langle{}d^{\\mathrm{TAD}}_{\\mathrm{n}}\\rangle{}~(\\sigma)$',fontsize=size)

for i in range(nsw):
    lines0=[]
    for j in range(1,nfr-9):
        normc=norm(t[j])
        rgb=cmap(normc)
        ax[i,2].plot(D_mean[i,j,:],'-',color=rgb,linewidth=wth)
    for j in range(ncl):
        ax[i,2].plot(D0_mean[i,j,:],'-',color=(0.5,0.5,0.5),linewidth=2*wth)
        lines0+=ax[i,2].plot(D0_mean[i,j,:],'-',color=color0[i][j],linewidth=wth,label=label0[i][j])
    ax[i,2].minorticks_on()
    ax[i,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[i,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax[i,2].xaxis.set_major_locator(MultipleLocator(xtick))
    ax[i,2].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].yaxis.set_major_locator(MultipleLocator(ytick[i][2]))
    ax[i,2].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[i,2].spines['bottom'].set_linewidth(wth)
    ax[i,2].spines['top'].set_linewidth(wth)
    ax[i,2].spines['left'].set_linewidth(wth)
    ax[i,2].spines['right'].set_linewidth(wth)
    ax[i,2].set_xlabel('$\\xi^{\\mathrm{TAD}}$',fontsize=size)
    ax[i,2].set_ylabel('$d^{\\mathrm{TAD}}~(\\sigma)$\nin %s'%label[i],fontsize=size)

ax2=ax[2,2].twinx()
lines0=[]
lines=[]
for j in range(nsw):
    ax[2,2].fill_between(t[1:-9],Slp_mean[j,1:-9]-Slp_std[j,1:-9],Slp_mean[j,1:-9]+Slp_std[j,1:-9],color=color[j],alpha=0.1)
    for k in range(ncl):
        if (j,k)!=(0,1):
            ax[2,2].plot([min(t[1:-9]),max(t[1:-9])],[Slp0_mean[j,k],Slp0_mean[j,k]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
            lines0+=ax[2,2].plot([min(t[1:-9]),max(t[1:-9])],[Slp0_mean[j,k],Slp0_mean[j,k]],color=color0[j][k],linestyle='-',linewidth=2*wth,label=label0[j][k])
    ax2.plot(t[1:-9],Slp_mean[j,1:-9],color='k',linestyle='-',linewidth=3*wth)
    lines+=ax2.plot(t[1:-9],Slp_mean[j,1:-9],color=color[j],linestyle='-',linewidth=2*wth,label=label[j])
    # ax[i].errorbar(t,r_mean[i,:,j],yerr=r_std[i,:,j],fmt='none',ecolor=colors[i],capsize=3)
ax[2,2].minorticks_on()
ax[2,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[2,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax[2,2].set_xscale('log')
# ax[2,2].xaxis.set_major_locator(MultipleLocator(xtick))
ax[2,2].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[2,2].yaxis.set_major_locator(MultipleLocator(ytick[2][2]))
ax[2,2].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[2,2].spines['bottom'].set_linewidth(wth)
ax[2,2].spines['top'].set_linewidth(wth)
ax[2,2].spines['left'].set_linewidth(wth)
ax[2,2].spines['right'].set_linewidth(wth)
ymax=max(np.max(Slp_mean[:,1:-9]),np.max(Slp0_mean))
ymin=min(np.min(Slp_mean[:,1:-9]),np.min(Slp0_mean))
ax[2,2].set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
ax[2,2].set_xlabel(f'$t~(\\tau)$',fontsize=size)
ax[2,2].set_ylabel('$\\mathit{\\Omega}~(\\sigma)$',fontsize=size)
ax2.set_yticks([])
ax2.set_ylim(ax[2,2].get_ylim())
ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)

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
