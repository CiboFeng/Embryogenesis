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
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='cpmt_strth.npz'
figname=file[:-4]
cl=['Z','ESC']
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
idx_cl=[0,4]

"""Read Data and Process"""
data=np.load(file)
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

"""Plot"""
fig,ax=plt.subplots(1,1,figsize=(9,5))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,0.6,1.0)],[(0.8,1.0,0.6),(0.8,1.0,0.6)]]
label0=[['ZP','ZM'],['ESC','ESC']]
color=['b','r']

for i in range(nsw):
    ax.plot(t[1:-9],CSs[i,1:-9,0],color=color[i],linestyle='-',linewidth=wth,label=f'{sw[i]}')
    ax.plot(t[1:-9],CSs[i,1:-9,1],color=color[i],linestyle=':',linewidth=wth,label=f'{sw[i]}')
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
# ymax=max(np.max(CSs_tot[:,1:-9]),np.max(CSs_tot_ref))
# ymin=min(np.min(CSs_tot[:,1:-9]),np.min(CSs_tot_ref))
# ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
ax.set_xscale('log')
ax.set_ylabel('Cpmt. Strth.',fontsize=size)
ax.legend(loc='lower right',fontsize=0.8*size)
ax.set_xlabel('Time ($\\tau$)',fontsize=size)

plt.tight_layout()
# plt.savefig(f'{figname}_tot.png',format='png')
# plt.savefig(f'{figname}_tot.pdf',format='pdf')
plt.show()
