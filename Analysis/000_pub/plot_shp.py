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
import matplotlib.patheffects as PathEffects
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='../shp/shp.pyw'
figname='shp'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
nsw=len(sw)
ncl=len(cl)
ntr=[1214,974]
nfr=56
sgm=4

"""Read Data"""
mean0=np.zeros((nsw,ncl,6))
for i in range(nsw):
    q=pyw(file,f'The Index of Cells, the Length of Principal Axis 1, 2, 3, Radius of Gyration, Aspherical Quantity, and '
               f'S ({sw[i]}):')
    mean0[i,:,0]=np.array(q[1])/sgm
    mean0[i,:,1]=np.array(q[2])/sgm
    mean0[i,:,2]=np.array(q[3])/sgm
    mean0[i,:,3]=np.array(q[4])/sgm
    mean0[i,:,4]=np.array(q[5])
    mean0[i,:,5]=np.array(q[6])

mean=np.zeros((nsw,nfr,6))
std=np.zeros((nsw,nfr,6))
for i in range(nsw):
    PA=np.zeros((nfr,ntr[i],3))
    Rg=np.zeros((nfr,ntr[i]))
    Dlt=np.zeros((nfr,ntr[i]))
    S=np.zeros((nfr,ntr[i]))
    q=pyw(file,f'The Index of Frames, The Index of Trajectories, the Length of Principal Axis 1, 2, 3, Radius of '
               f'Gyration, Aspherical Quantity, and S ({sw[i]}):')
    for j in range(len(q[0])):
        PA[int(q[0][j]),int(q[1][j]),0]=q[2][j]/sgm
        PA[int(q[0][j]),int(q[1][j]),1]=q[3][j]/sgm
        PA[int(q[0][j]),int(q[1][j]),2]=q[4][j]/sgm
        Rg[int(q[0][j]),int(q[1][j])]=q[5][j]/sgm
        Dlt[int(q[0][j]),int(q[1][j])]=q[6][j]
        S[int(q[0][j]),int(q[1][j])]=q[7][j]
    mean[i,:,:3]=np.mean(PA,axis=1)
    std[i,:,:3]=np.std(PA,axis=1)
    mean[i,:,3]=np.mean(Rg,axis=1)
    std[i,:,3]=np.std(Rg,axis=1)
    mean[i,:,4]=np.mean(Dlt,axis=1)
    std[i,:,4]=np.std(Dlt,axis=1)
    mean[i,:,5]=np.mean(S,axis=1)
    std[i,:,5]=np.std(S,axis=1)

"""Generate Times"""
t0=[0,0.01,
    0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
    0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
    2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
    20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
    200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
    2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
t=np.logspace(-5,4,10000)

"""Define a Function for Interpolating"""
def itplt(t0,x0,t):
   f=interpolate.interp1d(t0,x0,kind='cubic')
   x=f(t)
   return x

"""Plot"""
fig,ax=plt.subplots(2,3,figsize=(25,13))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=0.1
ytick=[[],[]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0.0,1.0,1.0),(1.0,0.3,0.3)]
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
axlabel=['R_{\\mathrm{PA1}}~(\\sigma)',
       'R_{\\mathrm{PA2}}~(\\sigma)',
       'R_{\\mathrm{PA3}}~(\\sigma)',
       'R_{\\mathrm{g}}~(\\sigma)',
       '\\mathit{\\Delta}',
       'S']
nax=len(axlabel)
ann=['(A)','(B)','(C)','(D)','(E)','(F)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

for i in range(nax):
    j,k=divmod(i,3)
    for l in range(nsw):
        ax[j,k].fill_between(t0[1:-9],mean[l,1:-9,i]-std[l,1:-9,i],mean[l,1:-9,i]+std[l,1:-9,i],color=color[l],alpha=0.1)
        for m in range(ncl):
            ax[j,k].plot([min(t0[1:-9]),max(t0[1:-9])],[mean0[l,m,i],mean0[l,m,i]],color=(0.5,0.5,0.5),linestyle='-',linewidth=3*wth)
            ax[j,k].plot([min(t0[1:-9]),max(t0[1:-9])],[mean0[l,m,i],mean0[l,m,i]],color=color0[l][m],linestyle='-',linewidth=2*wth)
    for l in range(nsw):
        ax[j,k].plot(t0[1:-9],mean[l,1:-9,i],color='k',linewidth=3*wth)
        ax[j,k].plot(t0[1:-9],mean[l,1:-9,i],color=color[l],linewidth=2*wth)
    ax[j,k].autoscale()
    ax[j,k].minorticks_on()
    ax[j,k].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[j,k].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[j,k].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[j,k].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[j,k].spines['bottom'].set_linewidth(wth)
    ax[j,k].spines['top'].set_linewidth(wth)
    ax[j,k].spines['left'].set_linewidth(wth)
    ax[j,k].spines['right'].set_linewidth(wth)
    ax[j,k].set_ylabel('$%s$' % axlabel[i],fontsize=size)
    ax[j,k].set_xscale('log')
lines0=[]
lines=[]
for l,m in [(0,0),(1,0),(1,1)]:
    lines0+=ax[0,1].plot([],[],color=color0[m][l],linestyle='-',linewidth=wth,label=f'{label0[l][m]}')
for l in range(nsw):
    lines+=ax[0,0].plot([],[],color=color[l],linewidth=wth,label=f'{label[l]}')
legend=ax[0,1].legend(loc='upper right',handlelength=0,handletextpad=0,fontsize=0.8*size)
for text,line in zip(legend.get_texts(),lines0):
    text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
    text.set_color(line.get_color())
legend=ax[0,0].legend(loc='upper right',handlelength=0,handletextpad=0,fontsize=0.8*size)
for text,line in zip(legend.get_texts(),lines):
    text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground='k')])
    text.set_color(line.get_color())
ax[1,0].set_xlabel('$t~(\\tau)$',fontsize=size)
ax[1,1].set_xlabel('$t~(\\tau)$',fontsize=size)
ax[1,2].set_xlabel('$t~(\\tau)$',fontsize=size)

plt.tight_layout()
# plt.subplots_adjust(left=0.2,right=0.85,top=0.95,bottom=0.2)
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()