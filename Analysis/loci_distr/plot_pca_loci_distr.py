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
file='pca_loci_distr.pyw'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl=len(cl)
nfr=56
ts=-2
te=3
nt=10000

"""Read Data"""
perct_=pyw(file,'Percentage')
perct=np.zeros((2,3))
for i in range(2):
    for j in range(3):
        perct[i,j]=perct_[j*2+i][0]
pc0=np.zeros((nsw,ncl,2,3))
pc=np.zeros((nsw,nfr,2,3))
for i in range(nsw):
    pc0[i,:,0,0]=pyw(file,f'({sw[i]}, Reference):')[0]
    pc0[i,:,1,0]=pyw(file,f'({sw[i]}, Reference):')[1]
    pc0[i,:,0,1]=pyw(file,f'({sw[i]}, Reference):')[2]
    pc0[i,:,1,1]=pyw(file,f'({sw[i]}, Reference):')[3]
    pc0[i,:,0,2]=pyw(file,f'({sw[i]}, Reference):')[4]
    pc0[i,:,1,2]=pyw(file,f'({sw[i]}, Reference):')[5]
    pc[i,:,0,0]=pyw(file,f'({sw[i]}):')[0]
    pc[i,:,1,0]=pyw(file,f'({sw[i]}):')[1]
    pc[i,:,0,1]=pyw(file,f'({sw[i]}):')[2]
    pc[i,:,1,1]=pyw(file,f'({sw[i]}):')[3]
    pc[i,:,0,2]=pyw(file,f'({sw[i]}):')[4]
    pc[i,:,1,2]=pyw(file,f'({sw[i]}):')[5]

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
   f=interpolate.interp1d(t0,x0,kind='linear')
   x=f(t)
   return x

"""Plot"""
wth=3
size=30
lenmaj=15
lenmin=8
xtick=3
ytick=3
color0=[(0.6,0.6,0.6),(0.0,0.0,0.0)]
marker0=[['s','o'],['^','o']]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=mcolors.LogNorm(vmin=10**ts,vmax=10**te)
subp=[131,132,133]
label=['A','B','Gene']

plt.figure(figsize=(18,6))
for i in range(3):
    ax=plt.subplot(subp[i])
    for j in range(nsw):
        plt.scatter(itplt(t,pc[j,:,1,i],t_),itplt(t,pc[j,:,0,i],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
        scatter=plt.scatter(itplt(t,pc[j,:,1,i],t_),itplt(t,pc[j,:,0,i],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
        # scatter=plt.scatter(pc[j,:,1,i],pc[j,:,0,i],c=t,cmap=cmap,norm=norm)
        for k in range(ncl):
            plt.scatter(pc0[j,k,1,i],pc0[j,k,0,i],color=color0[k],marker=marker0[j][k],s=40*wth)
            if (j,k)!=(0,1):
                plt.scatter([],[],color=color0[k],marker=marker0[j][k],s=40*wth,label=label0[j][k])
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
    ax.set_xlabel(f'PC2 of {label[i]} ({tot_len(perct[1,i]*100,3)}%)',fontsize=size)
    ax.set_ylabel(f'PC1 of {label[i]} ({tot_len(perct[0,i]*100,3)}%)',fontsize=size)
    if i==0:
        plt.legend(loc='best',fontsize=size)
    if i==2:
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