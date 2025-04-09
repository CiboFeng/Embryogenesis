"""Import Modules"""
import numpy as np
from scipy import interpolate
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Arguments"""
file1='../simi_cont_prob/r2_cont_prob.npy'
file2='../dim_rdc_obs_exp/pca_obs_exp.pyw'
file3='../dim_rdc_insul_score/pca_insul_score.pyw'
figname='valid'
cl1=['Z','8C','ESC']
cl2=['ZP_','ZM_','8CP','8CM','ESC_']
sw=['ZP-ESC','ZM-ESC']
nsw=len(sw)
ncl1=len(cl1)
ncl2=len(cl2)
nfr=56
idx_cl=[0,3,4]
ts=-2
te=3
nt=1000
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
t_=np.logspace(ts,te,nt)

"""Read Data"""
r2=np.load(file1)
# r2_cross=np.zeros((nsw,nsw-1,nfr))
# r2_cross[0,0]=np.diagonal(np.load(file2)[0,1])
# r2_cross[1,0]=np.diagonal(np.load(file2)[1,0])

perct1=pyw(file2,'The Percentage of the First Component:')
perct2=pyw(file3,'The Percentage of the First Component:')
pc10=[]
pc20=[]
for i in range(ncl2):
    pc10.append(pyw(file2,f'Component of {cl2[i]}:')[0])
    pc20.append(pyw(file3,f'Component of {cl2[i]}:')[0])
pc1=[]
pc2=[]
for i in range(nsw):
    pc1.append(pyw(file2,f'Component of {sw[i]}:')[0])
    pc2.append(pyw(file3,f'Component of {sw[i]}:')[0])

"""Define a Function for Interpolating"""
def itplt(t,x,t_):
   f=interpolate.interp1d(t,x,kind='cubic')
   x_=f(t_)
   return x_

"""Plot"""
fig=plt.figure(figsize=(25,9))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
subp=[(0,0),(1,0),(0,1)]
ann=['(B)','(C)','(D)']
markerstyle=[]
color0=[[(1.0,0.8,0.6),(0.9,0.9,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,0.8,0.8),(0.8,1.0,0.6)]]
label0=[['ZP','8CP','ESC'],['ZM','8CM','ESC']]
label=['ZP$\\rightarrow$ESC','ZM$\\rightarrow$ESC']
coloredge=['k',(0.8,0.0,0.0),'k']

for i in range(nsw):
    ax1=plt.subplot2grid((2,2),subp[i],colspan=1,rowspan=1)
    ax1.annotate(ann[i],xy=(0.0,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    ax2=ax1.twinx()
    lines=[]
    for j in [0,2]:
        ax1.plot(t[1:-9],r2[i,idx_cl[j],i,1:-9],'-',color=coloredge[j],linewidth=4*wth,linestyle='-')
        ax1.plot(t[1:-9],r2[i,idx_cl[j],i,1:-9],'-',color=color0[i][j],linewidth=2*wth,linestyle='-')
    for j in range(ncl1):
        if j in [0,2]:
            ax2.plot([],[],'-',color=coloredge[j],linewidth=4*wth,linestyle='-')
            lines+=ax2.plot([],[],'-',color=color0[i][j],linewidth=2*wth,linestyle='-',label=label0[i][j])
        if j in [1]:
            ax2.plot(t[1:-9],r2[i,idx_cl[j],i,1:-9],'-',color=coloredge[j],linewidth=4*wth,linestyle='-')
            lines+=ax2.plot(t[1:-9],r2[i,idx_cl[j],i,1:-9],'-',color=color0[i][j],linewidth=2*wth,linestyle='-',label=label0[i][j])
    # plt.plot(t[:-9],r2_cross[i,0,:-9],'-',color='k',linewidth=wth,label=sw[-i-1])
    plt.xscale('log')
    ax1.minorticks_on()
    ax1.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax1.tick_params(axis='x',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax1.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # plt.xlim(115,225)
    # plt.ylim(0.92,1)
    # ax1.xaxis.set_major_locator(MultipleLocator(xtick))
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax1.yaxis.set_major_locator(MultipleLocator(ytick))
    ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
    axis=plt.gca()
    axis.spines['bottom'].set_linewidth(wth)
    axis.spines['top'].set_linewidth(wth)
    axis.spines['left'].set_linewidth(wth)
    axis.spines['right'].set_linewidth(wth)
    axis.spines['right'].set_color(coloredge[1])
    if i==1:
        ax1.set_xlabel(f'$t~(\\tau)$',fontsize=size)
    ax1.set_ylabel('$R^2(P,P^{\\mathrm{\\{%s,ESC\\}}})$\nin %s'%(label0[i][0],label[i]),fontsize=size)
    ax2.minorticks_on()
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size,colors=coloredge[1])
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size,colors=coloredge[1])
    # ax2.xaxis.set_major_locator(MultipleLocator(xtick))
    # ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax2.yaxis.set_major_locator(MultipleLocator(ytick))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.set_ylabel('$R^2(P,P^{\\mathrm{%s}})$\nin %s'%(label0[i][1],label[i]),fontsize=size)
    ax2.yaxis.label.set_color(coloredge[1])
    legend=ax2.legend(bbox_to_anchor=(1.0,0.9),loc='upper right',handlelength=0.0,handletextpad=0.0,columnspacing=0.5,fontsize=0.8*size,ncol=len(lines))
    j=0
    for text,line in zip(legend.get_texts(),lines):
        text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=coloredge[j])])
        text.set_color(line.get_color())
        j+=1

color0=[(1.0,0.8,0.6),(0.8,0.6,1.0),(0.9,0.9,0.6),(0.8,0.8,0.8),(0.8,1.0,0.6),(0.8,1.0,0.6)]
label0=['ZP','ZM','8CP','8CM','ESC']
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=mcolors.LogNorm(vmin=10**ts,vmax=10**te)
ax=plt.subplot2grid((2,2),subp[nsw],colspan=1,rowspan=2)
ax.annotate(ann[nsw],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
scatters0=[]
for i in range(nsw):
    plt.scatter(itplt(t,pc2[i],t_),itplt(t,pc1[i],t_),c='k',s=20*wth,cmap=cmap,norm=norm)
    scatter=plt.scatter(itplt(t,pc2[i],t_),itplt(t,pc1[i],t_),c=t_,s=10*wth,cmap=cmap,norm=norm)
for i in range(ncl2):
    plt.scatter(pc20[i],pc10[i],color=(0.5,0.5,0.5),marker='^',s=150*wth)
    scatters0+=[plt.scatter(pc20[i],pc10[i],color=color0[i],marker='^',s=100*wth,label=label0[i])]
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
plt.xlabel(f'PC1 of $S_{{\\mathrm{{i}}}}$',fontsize=size)
plt.ylabel(f'PC1 of $S_{{\\mathrm{{C}}}}$',fontsize=size)
empty_handles=[Line2D([],[],color='none') for _ in scatters0]
legend=ax.legend(empty_handles,label0,loc='best',handlelength=0,handletextpad=0,fontsize=0.8*size)
for text,line in zip(legend.get_texts(),scatters0):
    text.set_path_effects([PathEffects.withStroke(linewidth=wth,foreground=(0.5,0.5,0.5))])
    text.set_color(line.get_facecolor()[0])

cbar_ax=fig.add_axes([0.9,0.125,0.02,0.774])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.set_label('$t~(\\tau)$',fontsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0,0.9,1])
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()