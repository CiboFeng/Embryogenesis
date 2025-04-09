"""Import Modules"""
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
from pyw import pyw
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

"""Set Arguments"""
file='rmsd.pyw'
figname=file[:-4]
ws=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
nw=len(ws)
sgm=4

"""Read Data"""
q=pyw(file,dir=True)
t=q[:,:,:,0]
d=q[:,:,:,1]/sgm
dm=np.mean(d,axis=1)

"""Plot"""
fig,ax=plt.subplots(3,1,figsize=(8,10))
wth=2
size=20
lenmaj=15
lenmin=8
xtick=0.1
ytick=20
# ytick=0.002
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=min(ws),vmax=max(ws))

for i in range(nw):
    normc=norm(ws[i])
    rgb=cmap(normc)
    ax[0].plot(t[0,i],d[0,i],color=rgb,linewidth=wth,alpha=0.5)
ax[0].plot(t[0,0],dm[0],color='k',linewidth=wth)
ax[0].autoscale()
ax[0].minorticks_on()
ax[0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[0].spines['bottom'].set_linewidth(wth)
ax[0].spines['top'].set_linewidth(wth)
ax[0].spines['left'].set_linewidth(wth)
ax[0].spines['right'].set_linewidth(wth)
# ax[0].set_xlabel('Time (ps)',fontsize=size)
ax[0].set_ylabel('RMSD ZP ($\\sigma$)',fontsize=size)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
# ax[0].yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))
ticks=[2,3,4,5,6]
ax[0].yaxis.set_major_locator(ticker.FixedLocator(ticks))
ax[0].yaxis.set_major_formatter(ticker.FixedFormatter([str(tick) for tick in ticks]))
# ax[0].yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))

for i in range(nw):
    normc=norm(ws[i])
    rgb=cmap(normc)
    ax[1].plot(t[1,i],d[1,i],color=rgb,linewidth=wth,alpha=0.5)
ax[1].plot(t[1,0],dm[1],color='k',linewidth=wth)
ax[1].autoscale()
ax[1].minorticks_on()
ax[1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[1].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[1].spines['bottom'].set_linewidth(wth)
ax[1].spines['top'].set_linewidth(wth)
ax[1].spines['left'].set_linewidth(wth)
ax[1].spines['right'].set_linewidth(wth)
# ax[1].set_xlabel('Time (ps)',fontsize=size)
ax[1].set_ylabel('RMSD ZM ($\\sigma$)',fontsize=size)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
# ax[1].yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))
ticks=[2,3,4,5,6]
ax[1].yaxis.set_major_locator(ticker.FixedLocator(ticks))
ax[1].yaxis.set_major_formatter(ticker.FixedFormatter([str(tick) for tick in ticks]))
# ax[1].yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))

for i in range(nw):
    normc=norm(ws[i])
    rgb=cmap(normc)
    ax[2].plot(t[2,i],d[2,i],color=rgb,linewidth=wth,alpha=0.5)
ax[2].plot(t[2,0],dm[2],color='k',linewidth=wth)
ax[2].autoscale()
ax[2].minorticks_on()
ax[2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax[2].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[2].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[2].spines['bottom'].set_linewidth(wth)
ax[2].spines['top'].set_linewidth(wth)
ax[2].spines['left'].set_linewidth(wth)
ax[2].spines['right'].set_linewidth(wth)
ax[2].set_xlabel('Time (ps)',fontsize=size)
ax[2].set_ylabel('RMSD ESC ($\\sigma$)',fontsize=size)
ax[2].set_xscale('log')
ax[2].set_yscale('log')
# ax[2].yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))
ticks=[2,3,4,5,6]
ax[2].yaxis.set_major_locator(ticker.FixedLocator(ticks))
ax[2].yaxis.set_major_formatter(ticker.FixedFormatter([str(tick) for tick in ticks]))
# ax[2].yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()