"""Import Modules"""
import numpy as np
import seaborn as sns
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import LogLocator,FixedLocator,AutoMinorLocator
from matplotlib.pyplot import MultipleLocator,FormatStrFormatter
import matplotlib.patheffects as PathEffects
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
npyfile=f'scl_half_life_prob_pc_.npy'
pywfile=f'scl_half_life_prob_pc_.pyw'
figname=npyfile[:-4]
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
j=0
while scl0*(2**(j+1)-1)<=nat/2:
    scl.append([scl0*(2**j-1),scl0*(2**(j+1)-1)])
    j+=1
# scl=[[0,5],[5,10],[10,20],[20,50],[50,100],[100,300]]
# scl=[[0,5],[5,15],[15,35],[35,75],[75,155]]
nsc=len(scl)
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read Data"""
Q=np.load(npyfile)
t1mean=np.zeros((nsw,nsc))
t2mean=np.zeros((nsw,nsc))
for i in range(nsw):
   q=pyw(pywfile,f'({sw[i]}):')
   t1mean[i]=q[0]
   t2mean[i]=q[1]

x=np.arange(0,nsc,1)

"""Define a Function for Interpolating"""
def itplt(t,x,t_):
   f=interpolate.interp1d(t,x,kind='linear')
   x_=f(t_)
   return x_

"""Define a Function for Axis Scale"""
def symlog(x,a):
   cond=[x<-10**a,(x>=-10**a)&(x<10**a),x>=10**a]
   func=[lambda x:-np.log10(-x)+a-1,
         lambda x:x/10**a,
         lambda x:np.log10(x)-a+1]
   y=np.piecewise(x,cond,func)
   return y

def symexp(x,a):
   cond=[x<-10**a,(x>=-10**a)&(x<10**a),x>=10**a]
   func=[lambda x:-10**(-x+a-1),
         lambda x:10**a*x,
         lambda x:10**(x+a-1)]
   y=np.piecewise(x,cond,func)
   return y

""""""
Qmean=np.mean(symlog(Q,0),axis=-1)
Qstd=np.std(symlog(Q,0),axis=-1)

"""Plot"""
plt.figure(figsize=(16,12))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=3
ytick=3
subp=[221,222,223,224]

for i in range(nsc):
   ax=plt.subplot(subp[i])
   for j in range(nat):
      ax.plot(t[1:-9],symlog(Q[0,1:-9,i,j],0),'-',c='c',linewidth=wth/5,zorder=0)
   ax.fill_between(t[1:-9],Qmean[0,1:-9,i]-Qstd[0,1:-9,i],Qmean[0,1:-9,i]+Qstd[0,1:-9,i],color=(1.0,0.8,0.8),zorder=1)
   ax.plot(t[1:-9],Qmean[0,1:-9,i],'-',c='r',linewidth=wth,zorder=2)
   ax.plot([t[1],t[-10]],[0.5,0.5],'--',c='k',linewidth=wth,zorder=3)
   plt.xscale('log')
   ax.minorticks_on()
   ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
   ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
   ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
   # plt.ylim(-3.2,3.2)
   axis=plt.gca()
   # axis.xaxis.set_major_locator(MultipleLocator(ticks))
   axis.xaxis.set_minor_locator(AutoMinorLocator(2))
   # axis.yaxis.set_major_locator(MultipleLocator(ytick))
   # axis.yaxis.set_minor_locator(AutoMinorLocator(2))
   ymax=np.ceil(symlog(np.max(Q[0][1:-9,:,i]),0))
   ymin=np.ceil(symlog(np.min(Q[0][1:-9,:,i]),0))-1
   y_positions=np.arange(ymin,ymax+1,1)
   y_labels=['$-10^{%s}$'%int(np.log10(-symexp(y_positions[i],0))) for i in range(-int(ymin))]
   y_labels+=['0']
   y_labels+=['$10^{%s}$'%int(np.log10(symexp(y_positions[i],0))) for i in range(-int(ymin)+1,int(ymax-ymin)+1)]
   y_min,y_max=plt.ylim()
   filt_positions=[val for val in y_positions if (y_min<=val<=y_max)]
   filt_labels=[y_labels[i] for i in range(len(y_labels)) if (y_min<=y_positions[i]<=y_max)]
   y_positions=filt_positions
   y_labels=filt_labels
   plt.yticks(y_positions,y_labels)
   axis.spines['bottom'].set_linewidth(wth)
   axis.spines['top'].set_linewidth(wth)
   axis.spines['left'].set_linewidth(wth)
   axis.spines['right'].set_linewidth(wth)
   plt.xlabel(f'$t~(\\tau)$',fontsize=size)
   plt.ylabel('$\\frac{Q-Q_\\mathrm{s}}{Q_\\mathrm{e}-Q_\\mathrm{s}}$'+f'({int(scl[i][0]/10)}-{int(scl[i][1]/10)} Mb)',fontsize=size)

plt.tight_layout()
plt.savefig(f'{figname}_trj.png',format='png')
plt.savefig(f'{figname}_trj.pdf',format='pdf')
plt.show()

# plt.figure(figsize=(16,12))
# wth=3
# size=30
# lenmaj=15
# lenmin=8
# xtick=3
# ytick=3
# subp=[221,222,223,224]
#
# for i in range(nsc):
#    ax=plt.subplot(subp[i])
#    for j in range(nat):
#       ax.plot(t,symlog(Q[0,:,i,j],0),'-',c='c',linewidth=wth/5,zorder=0)
#    ax.fill_between(t,Qmean[0,:,i]-Qstd[0,:,i],Qmean[0,:,i]+Qstd[0,:,i],color=(1.0,0.8,0.8),zorder=1)
#    ax.plot(t,Qmean[0,:,i],'-',c='r',linewidth=wth,zorder=2)
#    ax.plot([t[0],t[-1]],[0.5,0.5],'--',c='k',linewidth=wth,zorder=3)
#    plt.xscale('log')
#    ax.minorticks_on()
#    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#    ax.tick_params(axis='x',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=0,labelsize=size)
#    # plt.ylim(-3.2,3.2)
#    axis=plt.gca()
#    # axis.xaxis.set_major_locator(MultipleLocator(ticks))
#    axis.xaxis.set_minor_locator(AutoMinorLocator(2))
#    # axis.yaxis.set_major_locator(MultipleLocator(ytick))
#    # axis.yaxis.set_minor_locator(AutoMinorLocator(2))
#    ymax=np.ceil(symlog(np.max(Q[0][:,:,i]),0))
#    ymin=np.ceil(symlog(np.min(Q[0][:,:,i]),0))-1
#    y_positions=np.arange(ymin,ymax+1,1)
#    y_labels=['$-10^{%s}$'%int(np.log10(-symexp(y_positions[i],0))) for i in range(-int(ymin))]
#    y_labels+=['0']
#    y_labels+=['$10^{%s}$'%int(np.log10(symexp(y_positions[i],0))) for i in range(-int(ymin)+1,int(ymax-ymin)+1)]
#    y_min,y_max=plt.ylim()
#    filt_positions=[val for val in y_positions if (y_min<=val<=y_max)]
#    filt_labels=[y_labels[i] for i in range(len(y_labels)) if (y_min<=y_positions[i]<=y_max)]
#    y_positions=filt_positions
#    y_labels=filt_labels
#    plt.yticks(y_positions,y_labels)
#    axis.spines['bottom'].set_linewidth(wth)
#    axis.spines['top'].set_linewidth(wth)
#    axis.spines['left'].set_linewidth(wth)
#    axis.spines['right'].set_linewidth(wth)
#    plt.xlabel(f'$t~(\\tau)$',fontsize=size)
#    plt.ylabel('$\\frac{Q-Q_\\mathrm{s}}{Q_\\mathrm{e}-Q_\\mathrm{s}}$'+f'({int(scl[i][0]/10)}-{int(scl[i][1]/10)} Mb)',fontsize=size)
#
# plt.tight_layout()
# plt.savefig(f'{figname}_trj.png',format='png')
# plt.savefig(f'{figname}_trj.pdf',format='pdf')
# plt.show()

plt.figure(figsize=(10,8))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=3
ytick=3
color=[(0.0,1.0,1.0),(1.0,0.3,0.3)]
color_=['g','b']

ax1=plt.subplot(111)
ax2=ax1.twinx()
lines=[]
for i in range(nsw):
   ax1.bar(x+(2*i-1)*0.2,t1mean[i],width=0.4,color=color[i],edgecolor=color_[0],linewidth=wth)
   ax2.plot(x+(2*i-1)*0.2,t2mean[i]/t1mean[i]**2,'.',color=color_[1],markersize=6*wth)
   ax2.plot(x+(2*i-1)*0.2,t2mean[i]/t1mean[i]**2,'-',color=color_[1],linewidth=2*wth)
   ax2.plot(x+(2*i-1)*0.2,t2mean[i]/t1mean[i]**2,'.',color=color[i],markersize=4*wth)
   lines+=ax2.plot(x+(2*i-1)*0.2,t2mean[i]/t1mean[i]**2,'-',color=color[i],linewidth=wth,label=sw[i])
ax1.minorticks_on()
ax1.tick_params(axis='x',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax1.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=color_[0],labelcolor=color_[0])
ax1.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax1.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=color_[0],labelcolor=color_[0])
# ax1.xaxis.set_major_locator(MultipleLocator(xtick))
ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax1.yaxis.set_major_locator(MultipleLocator(ytick))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
axis=plt.gca()
axis.spines['bottom'].set_linewidth(wth)
axis.spines['top'].set_linewidth(wth)
axis.spines['left'].set_linewidth(wth)
axis.spines['right'].set_linewidth(wth)
axis.spines['left'].set_color(color_[0])
axis.spines['right'].set_color(color_[1])
xtick=x
xlabel=['0-2Mb','2-6Mb','6-14Mb','14-30Mb']
# xlabel=['0-0.5Mb','0.5-1Mb','1-2Mb','2-5Mb','5-10Mb','10-30Mb']
# xlabel=['0-0.5Mb','0.5-1.5Mb','1.5-3.5Mb','3.5-7.5Mb','7.5-15.5Mb']
ax1.set_xticks(xtick)
ax1.set_xticklabels(xlabel,rotation=45,ha='right')
ax1.set_xlabel(f'Distance Range',fontsize=size)
ax1.set_ylabel('$\\langle{}T_{1/2}\\rangle{}~(\\tau)$',fontsize=size)
ax1.set_yscale('log')
ax1.yaxis.label.set_color(color_[0])
ax2.minorticks_on()
ax2.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,color=color_[1],labelcolor=color_[1])
ax2.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,color=color_[1],labelcolor=color_[1])
# ax2.xaxis.set_major_locator(MultipleLocator(xtick))
# ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax2.yaxis.set_major_locator(MultipleLocator(ytick))
ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
ax2.set_ylabel('$\\langle{}T_{1/2}^2\\rangle{}/\\langle{}T_{1/2}\\rangle{}^2$',fontsize=size)
ax2.yaxis.label.set_color(color_[1])
legend=ax2.legend(loc='upper center',handlelength=0.0,handletextpad=0.0,columnspacing=0.5,fontsize=0.8*size,ncol=1)
for text,line in zip(legend.get_texts(),lines):
   # text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground='k')])
   text.set_color(line.get_color())
   text.set_fontsize(size)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()
