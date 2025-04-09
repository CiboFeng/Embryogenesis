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
file='umap_obs_exp.pyw'
figname=file[:-4]
cl=['ZP_','ZM_','ESC_','ZP','ZM','E2CP','E2CM','L2CP','L2CM','8CP','8CM','ESC']
# cl=['ZP','ZM','ESC']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl=len(cl)
ts=-2
te=3
nt=1000

"""Read Data"""
pc10=[]
pc20=[]
for i in range(ncl):
    pc10.append(pyw(file,f'Component of {cl[i]}:')[0])
    pc20.append(pyw(file,f'Component of {cl[i]}:')[1])
pc1=[]
pc2=[]
for i in range(nsw):
    pc1.append(pyw(file,f'Component of {sw[i]}:')[0])
    pc2.append(pyw(file,f'Component of {sw[i]}:')[1])

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
def itplt(t,x,t_):
   f=interpolate.interp1d(t,x,kind='cubic')
   x_=f(t_)
   return x_

"""Plot"""
plt.figure(figsize=(8,7))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=3
ytick=3
color0=[(0.8,0.8,0.8),(0.8,0.8,0.8),(0.0,0.0,0.0),(0.8,0.8,0.8),(0.8,0.8,0.8),(0.6,0.6,0.6),(0.6,0.6,0.6),(0.4,0.4,0.4),(0.4,0.4,0.4),(0.2,0.2,0.2),(0.2,0.2,0.2),(0.0,0.0,0.0)]
marker0=['p','*','h','s','^','s','^','s','^','s','^','o']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=mcolors.LogNorm(vmin=10**ts,vmax=10**te)

ax=plt.subplot(111)
for i in range(nsw):
    plt.scatter(itplt(t,pc2[i],t_),itplt(t,pc1[i],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
    scatter=plt.scatter(itplt(t,pc2[i],t_),itplt(t,pc1[i],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
for i in range(ncl):
    plt.scatter(pc20[i],pc10[i],color=color0[i],marker=marker0[i],s=40*wth,label=cl[i])
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# plt.xlim(115,225)
# plt.ylim(2,32)
axis=plt.gca()
axis.xaxis.set_minor_locator(AutoMinorLocator(2))
axis.yaxis.set_minor_locator(AutoMinorLocator(2))
axis.spines['bottom'].set_linewidth(wth)
axis.spines['top'].set_linewidth(wth)
axis.spines['left'].set_linewidth(wth)
axis.spines['right'].set_linewidth(wth)
plt.xlabel(f'PC2 of Obs. Exp',fontsize=size)
plt.ylabel(f'PC1 of Obs. Exp',fontsize=size)
plt.legend(loc='best',fontsize=size/2)
divider=make_axes_locatable(ax)
cax=divider.append_axes('right',size='5%',pad=0.05)
cbar=plt.colorbar(scatter,cax=cax)
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()