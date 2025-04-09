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
from matplotlib.lines import Line2D
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file='../shp/shp.pyw'
figname='shp_2d'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
nsw=len(sw)
ncl=len(cl)
ntr=[1214,974]
nfr=56
label=['R_{\\mathrm{PA1}}~(\\sigma)',
       'R_{\\mathrm{PA2}}~(\\sigma)',
       'R_{\\mathrm{PA3}}~(\\sigma)',
       'R_{\\mathrm{g}}~(\\sigma)',
       '\\Delta',
       'S']
nl=len(label)
sgm=4
ts=-2
te=3
nt=1000

"""Read Data"""
PA0=np.zeros((nsw,ncl,3))
Rg0=np.zeros((nsw,ncl))
Dlt0=np.zeros((nsw,ncl))
S0=np.zeros((nsw,ncl))
for i in range(nsw):
    q=pyw(file,f'The Index of Cells, the Length of Principal Axis 1, 2, 3, Radius of Gyration, Aspherical Quantity, and '
               f'S ({sw[i]}):')
    PA0[i,:,0]=np.array(q[1])/sgm
    PA0[i,:,1]=np.array(q[2])/sgm
    PA0[i,:,2]=np.array(q[3])/sgm
    Rg0[i]=np.array(q[4])/sgm
    Dlt0[i]=np.array(q[5])
    S0[i]=np.array(q[6])

PA=[]
Rg=[]
Dlt=[]
S=[]
PA_mean=np.zeros((nsw,nfr,3))
Rg_mean=np.zeros((nsw,nfr))
Dlt_mean=np.zeros((nsw,nfr))
S_mean=np.zeros((nsw,nfr))
for i in range(nsw):
    PA+=[np.zeros((nfr,ntr[i],3))]
    Rg+=[np.zeros((nfr,ntr[i]))]
    Dlt+=[np.zeros((nfr,ntr[i]))]
    S+=[np.zeros((nfr,ntr[i]))]
    q=pyw(file,f'The Index of Frames, The Index of Trajectories, the Length of Principal Axis 1, 2, 3, Radius of '
               f'Gyration, Aspherical Quantity, and S ({sw[i]}):')
    for j in range(len(q[0])):
        PA[i][int(q[0][j]),int(q[1][j]),0]=q[2][j]/sgm
        PA[i][int(q[0][j]),int(q[1][j]),1]=q[3][j]/sgm
        PA[i][int(q[0][j]),int(q[1][j]),2]=q[4][j]/sgm
        Rg[i][int(q[0][j]),int(q[1][j])]=q[5][j]/sgm
        Dlt[i][int(q[0][j]),int(q[1][j])]=q[6][j]
        S[i][int(q[0][j]),int(q[1][j])]=q[7][j]
    PA_mean[i]=np.mean(PA[i],axis=1)
    Rg_mean[i]=np.mean(Rg[i],axis=1)
    Dlt_mean[i]=np.mean(Dlt[i],axis=1)
    S_mean[i]=np.mean(S[i],axis=1)

"""Generate Times"""
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
t_=np.logspace(ts,te,nt)
t__=np.logspace(ts,te,int(nt/3))

"""Define a Function for Interpolating"""
def itplt(t,x,t_):
   f=interpolate.interp1d(t,x,kind='cubic')
   x_=f(t_)
   return x_

def itplt_chnn(x):
    f=interpolate.interp1d(t,x,kind='linear',bounds_error=False,fill_value="extrapolate")
    x_=f(t__)
    return x_

PA_=[]
Rg_=[]
Dlt_=[]
S_=[]
for i in range(nsw):
    PA_+=[np.apply_along_axis(itplt_chnn,0,PA[i]).reshape(int(nt/3),ntr[i],3)]
    Rg_+=[np.apply_along_axis(itplt_chnn,0,Rg[i]).reshape(int(nt/3),ntr[i])]
    Dlt_+=[np.apply_along_axis(itplt_chnn,0,Dlt[i]).reshape(int(nt/3),ntr[i])]
    S_+=[np.apply_along_axis(itplt_chnn,0,S[i]).reshape(int(nt/3),ntr[i])]

"""Plot"""
fig,ax=plt.subplots(2,3,figsize=(25,14))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=[[1.0,1.0,0.5],[0.1,0.1,0.02]]
ytick=[[1.0,1.0,0.5],[1.0,1.0,0.5]]
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=mcolors.LogNorm(vmin=10**ts,vmax=10**te)
ann=['(A)','(B)','(C)','(D)','(E)','(F)']
for axs,anns in zip(ax.flat,ann):
    axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

xmid1=[5.0,5.0]
ymid1=[3.25,3.0]
axwth1=[3.0,3.0]
xlim2=[[-0.02,0.3],[-0.02,0.3]]
ylim2=[[5.0,9.4],[5.0,8.6]]
for i in range(nsw):
    scatters0=[]
    for j in range(ntr[i]):
        ax[0,i].scatter(PA_[i][:,j,0],PA_[i][:,j,2],c=t__,cmap=cmap,norm=norm,s=wth/3,alpha=0.1)
    ax[0,i].scatter(itplt(t,PA_mean[i,:,0],t_),itplt(t,PA_mean[i,:,2],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
    scatter=ax[0,i].scatter(itplt(t,PA_mean[i,:,0],t_),itplt(t,PA_mean[i,:,2],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
    for j in range(nsw):
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[0,i].scatter(PA0[j,k,0],PA0[j,k,2],color=(0.5,0.5,0.5),marker='^',s=150*wth)
                scatters0+=[ax[0,i].scatter(PA0[j,k,0],PA0[j,k,2],color=color0[j][k],marker='^',s=100*wth,label=label0[j][k])]
    ax[0,i].minorticks_on()
    ax[0,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[0,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # xmax=max(np.max(PA_[i][:,:,0]),np.max(itplt(t,PA_mean[i,:,0],t_)),np.max(PA0[:,:,0]))
    # xmin=min(np.min(PA_[i][:,:,0]),np.min(itplt(t,PA_mean[i,:,0],t_)),np.min(PA0[:,:,0]))
    # xmid=(xmin+xmax)/2
    # ymax=max(np.max(PA_[i][:,:,2]),np.max(itplt(t,PA_mean[i,:,2],t_)),np.max(PA0[:,:,2]))
    # ymin=min(np.min(PA_[i][:,:,2]),np.min(itplt(t,PA_mean[i,:,2],t_)),np.min(PA0[:,:,2]))
    # ymid=(ymin+ymax)/2
    # axwth=max(xmax-xmin,ymax-ymin)
    # xmax=xmid+0.5*1.2*axwth
    # xmin=xmid-0.5*1.2*axwth
    # ymax=ymid+0.5*1.2*axwth
    # ymin=ymid-0.5*1.2*axwth
    # ax[0,i].set_xlim(xmin,xmax)
    # ax[0,i].set_ylim(ymin,ymax)
    bbox=ax[0,i].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axwth10,axht10=bbox.width,bbox.height
    axht=axwth1[i]*axht10/axwth10
    xmin1=xmid1[i]-0.5*axwth1[i]
    xmax1=xmid1[i]+0.5*axwth1[i]
    ymin1=ymid1[i]-0.5*axht
    ymax1=ymid1[i]+0.5*axht
    ax[0,i].set_xlim(xmin1,xmax1)
    ax[0,i].set_ylim(ymin1,ymax1)
    # ax[0,i].axis('equal')
    ax[0,i].xaxis.set_major_locator(MultipleLocator(xtick[0][i]))
    ax[0,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0,i].yaxis.set_major_locator(MultipleLocator(ytick[0][i]))
    ax[0,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0,i].spines['bottom'].set_linewidth(wth)
    ax[0,i].spines['top'].set_linewidth(wth)
    ax[0,i].spines['left'].set_linewidth(wth)
    ax[0,i].spines['right'].set_linewidth(wth)
    ax[0,i].set_xlabel('$R_{\\mathrm{PA1}}~(\\sigma)$ in %s'%label[i],fontsize=size)
    ax[0,i].set_ylabel('$R_{\\mathrm{PA3}}~(\\sigma)$ in %s'%label[i],fontsize=size)

    scatters0=[]
    for j in range(ntr[i]):
        ax[1,i].scatter(Dlt_[i][:,j],Rg_[i][:,j],c=t__,cmap=cmap,norm=norm,s=wth/3,alpha=0.1)
    ax[1,i].scatter(itplt(t,Dlt_mean[i,:],t_),itplt(t,Rg_mean[i,:],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
    scatter=ax[1,i].scatter(itplt(t,Dlt_mean[i,:],t_),itplt(t,Rg_mean[i,:],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
    for j in range(nsw):
        for k in range(ncl):
            if (j,k)!=(0,1):
                ax[1,i].scatter(Dlt0[j,k],Rg0[j,k],color=(0.5,0.5,0.5),marker='^',s=150*wth)
                scatters0+=[ax[1,i].scatter(Dlt0[j,k],Rg0[j,k],color=color0[j][k],marker='^',s=100*wth,label=label0[j][k])]
    ax[1,i].minorticks_on()
    ax[1,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[1,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[1,i].set_xlim(xlim2[i][0],xlim2[i][1])
    ax[1,i].set_ylim(ylim2[i][0],ylim2[i][1])
    ax[1,i].xaxis.set_major_locator(MultipleLocator(xtick[1][i]))
    ax[1,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1,i].yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
    ax[1,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[1,i].spines['bottom'].set_linewidth(wth)
    ax[1,i].spines['top'].set_linewidth(wth)
    ax[1,i].spines['left'].set_linewidth(wth)
    ax[1,i].spines['right'].set_linewidth(wth)
    ax[1,i].set_xlabel('$\\mathit{\\Delta}$ in %s'%label[i],fontsize=size)
    ax[1,i].set_ylabel('$R_{\\mathrm{g}}~(\\sigma)$ in %s'%label[i],fontsize=size)

scatters0=[]
for i in range(nsw):
    ax[0,2].scatter(itplt(t,PA_mean[i,:,0],t_),itplt(t,PA_mean[i,:,2],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
    scatter=ax[0,2].scatter(itplt(t,PA_mean[i,:,0],t_),itplt(t,PA_mean[i,:,2],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
    for j in range(ncl):
        if (i,j)!=(0,1):
            ax[0,2].scatter(PA0[i,j,0],PA0[i,j,2],color=(0.5,0.5,0.5),marker='^',s=150*wth)
            scatters0+=[ax[0,2].scatter(PA0[i,j,0],PA0[i,j,2],color=color0[i][j],marker='^',s=100*wth,label=label0[i][j])]
ax[0,2].minorticks_on()
ax[0,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[0,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# xmax=max(np.max(itplt(t,PA_mean[:,:,0],t_)),np.max(PA0[:,:,0]))
# xmin=min(np.min(itplt(t,PA_mean[:,:,0],t_)),np.min(PA0[:,:,0]))
# xmid=(xmin+xmax)/2
# ymax=max(np.max(itplt(t,PA_mean[:,:,2],t_)),np.max(PA0[:,:,2]))
# ymin=min(np.min(itplt(t,PA_mean[:,:,2],t_)),np.min(PA0[:,:,2]))
# ymid=(ymin+ymax)/2
# axwth=max(xmax-xmin,ymax-ymin)
# xmax=xmid+0.5*1.2*axwth
# xmin=xmid-0.5*1.2*axwth
# ymax=ymid+0.5*1.2*axwth
# ymin=ymid-0.5*1.2*axwth
# ax[0,2].set_xlim(xmin,xmax)
# ax[0,2].set_ylim(ymin,ymax)
ax[0,2].axis('equal')
ax[0,2].xaxis.set_major_locator(MultipleLocator(xtick[0][2]))
ax[0,2].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[0,2].yaxis.set_major_locator(MultipleLocator(ytick[0][2]))
ax[0,2].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[0,2].spines['bottom'].set_linewidth(wth)
ax[0,2].spines['top'].set_linewidth(wth)
ax[0,2].spines['left'].set_linewidth(wth)
ax[0,2].spines['right'].set_linewidth(wth)
ax[0,2].set_xlabel('$R_{\\mathrm{PA1}}~(\\sigma)$',fontsize=size)
ax[0,2].set_ylabel('$R_{\\mathrm{PA3}}~(\\sigma)$',fontsize=size)

scatters0=[]
for i in range(nsw):
    ax[1,2].scatter(itplt(t,Dlt_mean[i,:],t_),itplt(t,Rg_mean[i,:],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
    scatter=ax[1,2].scatter(itplt(t,Dlt_mean[i,:],t_),itplt(t,Rg_mean[i,:],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
    for j in range(ncl):
        if (i,j)!=(0,1):
            ax[1,2].scatter(Dlt0[i,j],Rg0[i,j],color=(0.5,0.5,0.5),marker='^',s=150*wth)
            scatters0+=[ax[1,2].scatter(Dlt0[i,j],Rg0[i,j],color=color0[i][j],marker='^',s=100*wth,label=label0[i][j])]
ax[1,2].minorticks_on()
ax[1,2].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax[1,2].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.set_xlim(115,225)
# ax.set_ylim(2,32)
ax[1,2].xaxis.set_major_locator(MultipleLocator(xtick[1][2]))
ax[1,2].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[1,2].yaxis.set_major_locator(MultipleLocator(ytick[1][2]))
ax[1,2].yaxis.set_minor_locator(AutoMinorLocator(2))
ax[1,2].spines['bottom'].set_linewidth(wth)
ax[1,2].spines['top'].set_linewidth(wth)
ax[1,2].spines['left'].set_linewidth(wth)
ax[1,2].spines['right'].set_linewidth(wth)
ax[1,2].set_xlabel('$\\mathit{\\Delta}$',fontsize=size)
ax[1,2].set_ylabel('$R_{\\mathrm{g}}~(\\sigma)$',fontsize=size)
empty_handles=[Line2D([],[],color='none') for _ in scatters0]
label0_=[label0[j][i] for i in range(nsw) for j in range(ncl) if (i,j)!=(0,1)]
legend=ax[0,0].legend(empty_handles,label0_,loc='best',handlelength=0,handletextpad=0,fontsize=0.8*size)
for text,line in zip(legend.get_texts(),scatters0):
    text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
    text.set_color(line.get_facecolor()[0])

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