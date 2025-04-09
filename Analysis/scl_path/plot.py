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
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import matplotlib.patheffects as PathEffects
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
npzfile1='scl_path_dist.npz'
npzfile2='scl_path_prob_mean.npz'
figname='scl_path'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
ntr=[1214,974]
scl0=20
scl=[]
k=0
while scl0*(2**(k+1)-1)<=nat:
    scl.append([scl0*(2**k-1),scl0*(2**(k+1)-1)])
    k+=1
scl.append([scl0*(2**k-1),nat])
nsc=len(scl)
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
ts=-2
te=3
nt=1000

"""Read Data"""
data=np.load(npzfile1)
Q0=100*data['Q0']
Q=[]
for i in range(nsw):
    Q+=[100*data[f'Q{i+1}']]

Q_mean=np.zeros((nsw,nfr,nsc))
for i in range(nsw):
   Q_mean[i]=np.mean(Q[i],axis=1)

data=np.load(npzfile2)
P0=100*data['Q0']
P=[]
for i in range(nsw):
    P+=[100*data[f'Q{i+1}']]

P_mean=np.zeros((nsw,nfr,nsc))
for i in range(nsw):
   P_mean[i]=np.mean(P[i],axis=1)

"""Define a Function for Interpolating"""
def itplt(t,x,t_):
   f=interpolate.interp1d(t,x,kind='cubic')
   x_=f(t_)
   return x_

"""Generate Times"""
t_=np.logspace(ts,te,nt)
# t_=t

"""Plot"""
fig,ax=plt.subplots(nsc,nsc,figsize=(30,23))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=3
ytick=3
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
axlabel=['0-2Mb','2-6Mb','6-14Mb','14-30Mb','30-60Mb']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=mcolors.LogNorm(vmin=10**ts,vmax=10**te)

for i in range(nsc):
   for j in range(nsc):
      if i>j:
         for k in range(nsw):
             ax[i,j].scatter(itplt(t,Q_mean[k,:,j],t_),itplt(t,Q_mean[k,:,i],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
             scatter=ax[i,j].scatter(itplt(t,Q_mean[k,:,j],t_),itplt(t,Q_mean[k,:,i],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
             for l in range(ncl):
                if (k,l)!=(0,1):
                   ax[i,j].scatter(Q0[k,l,j],Q0[k,l,i],color='k',marker='^',s=150*wth)
                   ax[i,j].scatter(Q0[k,l,j],Q0[k,l,i],color=color0[k][l],marker='^',s=100*wth,label=label0[k][l])
      elif i<j:
         for k in range(nsw):
            ax[i,j].scatter(itplt(t,P_mean[k,:,j],t_),itplt(t,P_mean[k,:,i],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
            scatter=ax[i,j].scatter(itplt(t,P_mean[k,:,j],t_),itplt(t,P_mean[k,:,i],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
            for l in range(ncl):
               if (k,l)!=(0,1):
                  ax[i,j].scatter(P0[k,l,j],P0[k,l,i],color='k',marker='^',s=150*wth)
                  ax[i,j].scatter(P0[k,l,j],P0[k,l,i],color=color0[k][l],marker='^',s=100*wth,label=label0[k][l])
      else:
         for k in range(nsw):
            ax[i,j].scatter(itplt(t,Q_mean[k,:,j],t_),itplt(t,P_mean[k,:,i],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
            scatter=ax[i,j].scatter(itplt(t,Q_mean[k,:,j],t_),itplt(t,P_mean[k,:,i],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
            for l in range(ncl):
               if (k,l)!=(0,1):
                  ax[i,j].scatter(Q0[k,l,j],P0[k,l,i],color='k',marker='^',s=150*wth)
                  ax[i,j].scatter(Q0[k,l,j],P0[k,l,i],color=color0[k][l],marker='^',s=100*wth,label=label0[k][l])
      ax[i,j].minorticks_on()
      ax[i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
      ax[i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
      # ax[i,j].xlim(115,225)
      # ax[i,j].ylim(2,32)
      ax[i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
      ax[i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
      ax[i,j].spines['bottom'].set_linewidth(wth)
      ax[i,j].spines['top'].set_linewidth(wth)
      ax[i,j].spines['left'].set_linewidth(wth)
      ax[i,j].spines['right'].set_linewidth(wth)
      if i>j:
         ax[i,j].set_xlabel('$Q_{\\mathrm{L}}\\times{}10^{-2}$ (%s)'%axlabel[j],fontsize=size)
         ax[i,j].set_ylabel('$Q_{\\mathrm{L}}\\times{}10^{-2}$ (%s)'%axlabel[i],fontsize=size)
      elif i<j:
         ax[i,j].set_xlabel('$P_{\\mathrm{L}}\\times{}10^{-2}$ (%s)'%axlabel[j],fontsize=size)
         ax[i,j].set_ylabel('$P_{\\mathrm{L}}\\times{}10^{-2}$ (%s)'%axlabel[i],fontsize=size)
      else:
         ax[i,j].set_xlabel('$Q_{\\mathrm{L}}\\times{}10^{-2}$ (%s)'%axlabel[j],fontsize=size)
         ax[i,j].set_ylabel('$P_{\\mathrm{L}}\\times{}10^{-2}$ (%s)'%axlabel[i],fontsize=size)
      if (i,j)==(0,0):
         ax[i,j].legend(loc='best',fontsize=0.8*size)

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
