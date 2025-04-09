"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
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
file='r2_cont_prob.npy'
# file2='r2_cont_prob_cross.npy'
figname=file[:-4]
cl=['Z','8C','ESC']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl=len(cl)
nfr=56
idx_cl=[0,3,4]

"""Read Data"""
r2=np.load(file)
# r2_cross=np.zeros((nsw,nsw-1,nfr))
# r2_cross[0,0]=np.diagonal(np.load(file2)[0,1])
# r2_cross[1,0]=np.diagonal(np.load(file2)[1,0])
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Plot"""
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
subp=[211,212]
markerstyle=[]
color0=['r','g','b']
label0=[['ZP','8CP','ESC'],['ZM','8CM','ESC']]

plt.figure(figsize=(10,10))
for i in range(nsw):
    ax1=plt.subplot(subp[i])
    ax2=ax1.twinx()
    for j in [0,2]:
        ax1.plot(t[1:-9],r2[i,idx_cl[j],i,1:-9],color=color0[j],linestyle='-',linewidth=wth)
    for j in range(ncl):
        if j in [0,2]:
            ax2.plot([],[],'-',color=color0[j],linestyle='-',linewidth=wth,label=label0[i][j])
        if j in [1]:
            ax2.plot(t[1:-9],r2[i,idx_cl[j],i,1:-9],'-',color=color0[j],linestyle='-',linewidth=wth,label=label0[i][j])
    # plt.plot(t[:-9],r2_cross[i,0,:-9],'-',color='k',linewidth=wth,label=sw[-i-1])
    plt.xscale('log')
    ax1.minorticks_on()
    ax1.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax1.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # plt.xlim(115,225)
    # plt.ylim(0.92,1)
    axis=plt.gca()
    # ax1.xaxis.set_major_locator(MultipleLocator(xtick))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax1.yaxis.set_major_locator(MultipleLocator(ytick))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
    axis.spines['bottom'].set_linewidth(wth)
    axis.spines['top'].set_linewidth(wth)
    axis.spines['left'].set_linewidth(wth)
    axis.spines['right'].set_linewidth(wth)
    axis.spines['right'].set_color(color0[1])
    ax1.set_xlabel(f'Time in {sw[i]}',fontsize=size)
    ax1.set_ylabel(f'$R^2$ of Cont. Prob.',fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,colors=color0[1])
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,colors=color0[1])
    # ax2.xaxis.set_major_locator(MultipleLocator(xtick))
    # ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax2.yaxis.set_major_locator(MultipleLocator(ytick))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.set_ylabel(f'$R^2$ of Cont. Prob.',fontsize=size)
    ax2.yaxis.label.set_color(color0[1])
    ax2.legend(loc='upper right',fontsize=0.8*size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()