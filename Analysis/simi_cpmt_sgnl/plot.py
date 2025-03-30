"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='r2_cpmt_sgnl.npy'
figname=file[:-4]
cl=['ZP','ZM','E2CP','E2CM','L2CP','L2CM','8CP','8CM','ESC']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl=len(cl)
nfr=56

"""Read Data"""
r2=np.load(file)
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
color=['r','r','g','g','b','b','m','m','c']
linestyle=[['-','--','-','--','-','--','-','--','-'],
           ['--','-','--','-','--','-','--','-','-']]

plt.figure(figsize=(10,10))
for i in range(nsw):
    ax=plt.subplot(subp[i])
    for j in range(ncl):
        plt.plot(t[:-9],r2[j,i,:-9],'-',color=color[j],linewidth=wth,label=cl[j],linestyle=linestyle[i][j])
    # plt.plot(t,r2_cross[i,0],'-',color='k',linewidth=wth,label=sw[-i-1])
    plt.xscale('log')
    ax.minorticks_on()
    ax.tick_params(axis="both",which="major",direction="in",width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis="both",which="minor",direction="in",width=wth,length=lenmin,labelsize=size)
    # plt.xlim(115,225)
    # plt.ylim(0.92,1)
    axis=plt.gca()
    # axis.xaxis.set_major_locator(MultipleLocator(xtick))
    # axis.xaxis.set_minor_locator(AutoMinorLocator(2))
    # axis.yaxis.set_major_locator(MultipleLocator(ytick))
    axis.yaxis.set_minor_locator(AutoMinorLocator(2))
    axis.spines['bottom'].set_linewidth(wth)
    axis.spines['top'].set_linewidth(wth)
    axis.spines['left'].set_linewidth(wth)
    axis.spines['right'].set_linewidth(wth)
    plt.xlabel(f'Time in {sw[i]}',fontsize=size)
    plt.ylabel(f'$R^2$ of Cpmt. Sgnl.',fontsize=size)
    plt.legend(loc='lower right',fontsize=size/2)
# ax=plt.subplot(326)
# plt.legend(loc='center',fontsize=size/4)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()