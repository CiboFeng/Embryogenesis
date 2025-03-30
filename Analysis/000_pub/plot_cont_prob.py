"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from split_data import split_data
from color_bar import color_bar

"""Set Arguments"""
npyfile='../cont_prob/cont_prob_rdc.npy'
reffile=[[],[]]
reffile[0]+=['../cont_prob/ref/zp.dat']
reffile[0]+=['../cont_prob/ref/zm.dat']
reffile[1]+=['../cont_prob/ref/esc.dat']
reffile[1]+=['../cont_prob/ref/esc.dat']
figname='cont_prob'
cl=[['ZP','ZM'],['ESC','ESC']]
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
cl_=['zp','zm','esc']
sw_=['zp_esc','zm_esc']
nfr=56
nat=600
binsize=1000000
t=[0,0.01,0.1,1.0,10.0,100.0,1000.0,10000.0]
t_show=[0.1,1.0,10.0,100.0]

"""Read Data, Make Some Modifications in Order to Make Plot More Convenient"""
P0=np.ones((ncl,nat,nat))
for i in range(nsw):
        if i==0:
                for j in range(ncl):
                        f=open(reffile[j][i],'r')
                        lsf=f.readlines()
                        for k in range(len(lsf)):
                                ls_f=lsf[k].strip('\n').split()
                                l,m=(int(ls_f[0])-1,int(ls_f[1])-1)
                                if l>m:
                                        P0[j,l,m]=float(ls_f[2])
        else:
                for j in range(ncl):
                        f=open(reffile[j][i],'r')
                        lsf=f.readlines()
                        for k in range(len(lsf)):
                                ls_f=lsf[k].strip('\n').split()
                                l,m=(int(ls_f[0])-1,int(ls_f[1])-1)
                                if l<m:
                                        P0[j,l,m]=float(ls_f[2])

p=np.load(npyfile)
P=np.ones((len(t_show),nat,nat))
for i in range(nsw):
        if i==0:
                for j in range(len(t_show)):
                        for k in range(nat):
                                for l in range(k):
                                        P[j,k,l]=p[i,t.index(t_show[j])+1,k,l]
        else:
                for j in range(len(t_show)):
                        for k in range(nat):
                                for l in range(k+1,nat):
                                        P[j,k,l]=p[i,t.index(t_show[j])+1,k,l]

"""Plot"""
# fig,ax=plt.subplots(1,2,figsize=(14,8))
# wth=2
# size=20
# lenmaj=12
# lenmin=8
# lenbar=8
# xtick=100
# ytick=100
# color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
# nodes=[0/2,1/2,2/2]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.LogNorm(vmin=1e-3,vmax=1)
#
# for i in range(ncl):
#         img=ax[i].imshow(P0[i],cmap=cmap,norm=norm)
#         x_positions=[0.5,200.5,400.5,599.5]
#         y_positions=[0.5,200.5,400.5,599.5]
#         x_labels=[0,20,40,60]
#         y_labels=[0,20,40,60]
#         ax[i].xaxis.set_minor_locator(MultipleLocator(xtick))
#         ax[i].yaxis.set_minor_locator(MultipleLocator(ytick))
#         ax[i].set_xticks(x_positions,x_labels)
#         ax[i].set_yticks(y_positions,y_labels)
#         ax[i].set_xlabel('Base Pairs in Mb',fontsize=size)
#         ax[i].set_ylabel('Base Pairs in Mb',fontsize=size)
#         ax[i].invert_yaxis()
#         ax[i].minorticks_on()
#         ax[i].tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
#         ax[i].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
#         for spine in ax[i].spines.values():
#                 spine.set_linewidth(wth)
#         ax[i].figure.canvas.draw()
#
# fig.subplots_adjust(right=0.85)
# cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
# sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
# cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
# cbar.outline.set_linewidth(wth)
#
# plt.tight_layout(rect=[0,0,0.85,1])
# plt.savefig(f'{figname}_ref.png',format='png')
# plt.savefig(f'{figname}_ref.pdf',format='pdf')
#
# fig,ax=plt.subplots(1,4,figsize=(14,4))
# wth=2
# size=20
# lenmaj=12
# lenmin=8
# lenbar=8
# color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
# nodes=[0/2,1/2,2/2]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.LogNorm(vmin=1e-3,vmax=1)
#
# for i in range(len(t_show)):
#         img=ax[i].imshow(P[i],cmap=cmap,norm=norm)
#         ax[i].set_xticks([],[])
#         ax[i].set_yticks([],[])
#         ax[i].invert_yaxis()
#         ax[i].minorticks_on()
#         ax[i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#         ax[i].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
#         for spine in ax[i].spines.values():
#                 spine.set_linewidth(wth)
#
# fig.subplots_adjust(right=0.85)
# cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
# sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
# sm.set_array([])
# cbar=plt.colorbar(sm,cax=cbar_ax)
# # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
# cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
# cbar.outline.set_linewidth(wth)
#
# plt.tight_layout(rect=[0,0,0.85,1])
# plt.savefig(f'{figname}.png',format='png')
# plt.savefig(f'{figname}.pdf',format='pdf')

fig,ax=plt.subplots(1,6,figsize=(25,5))
wth=3
size=30
lenmaj=15
lenmin=8
lenbar=8
xtick=100
ytick=100
# color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
# nodes=[0/2,1/2,2/2]
color=[(0.5,1.0,1.0),(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
nodes=[0,1/6,4/6,1]
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/16,4/16,9/16,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=1e-3,vmax=1)
# ann=['(A)','(B)','(C)','(D)','(E)','(F)']
# for axs,anns in zip(ax.flat,ann):
#     axs.annotate(anns,xy=(0.2,1.1),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
ax[0].annotate('(A)',xy=(0.2,1.1),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

top_label=['Z','ESC']
for i in range(ncl):
        I=[0,-1]
        img=ax[I[i]].imshow(P0[i],cmap=cmap,norm=norm)
        x_positions=[0.5,200.5,400.5,599.5]
        y_positions=[0.5,200.5,400.5,599.5]
        x_labels=[0,20,40,60]
        y_labels=[0,20,40,60]
        ax[I[i]].xaxis.set_minor_locator(MultipleLocator(xtick))
        ax[I[i]].set_xticks(x_positions,x_labels)
        ax[I[i]].set_xlabel('Loci (Mb)',fontsize=size)
        if I[i]==0:
                ax[I[i]].yaxis.set_minor_locator(MultipleLocator(ytick))
                ax[I[i]].set_yticks(y_positions,y_labels)
                ax[I[i]].set_ylabel('Base Pairs in Mb',fontsize=size)
        else:
                # ax[I[i]].yaxis.set_minor_locator(MultipleLocator(ytick))
                ax[I[i]].set_yticks([],[])
                # ax[I[i]].set_ylabel('Base Pairs in Mb',fontsize=size)
        ax[I[i]].invert_yaxis()
        ax[I[i]].minorticks_on()
        ax[I[i]].tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
        ax[I[i]].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
        for spine in ax[I[i]].spines.values():
                spine.set_linewidth(wth)
        ax[I[i]].text(0.05,0.95,'Pat',transform=ax[I[i]].transAxes,verticalalignment='top',horizontalalignment='left',fontsize=size)
        ax[I[i]].text(0.95,0.05,'Mat',transform=ax[I[i]].transAxes,verticalalignment='bottom',horizontalalignment='right',fontsize=size)
        ax[I[i]].text(0.5,1.05,f'{top_label[i]}',transform=ax[I[i]].transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=size)

for i in range(len(t_show)):
        img=ax[i+1].imshow(P[i],cmap=cmap,norm=norm)
        x_positions=[0.5,200.5,400.5,599.5]
        y_positions=[0.5,200.5,400.5,599.5]
        x_labels=[0,20,40,60]
        y_labels=[0,20,40,60]
        ax[i+1].xaxis.set_minor_locator(MultipleLocator(xtick))
        ax[i+1].set_xticks(x_positions,x_labels)
        ax[i+1].set_xlabel('Loci (Mb)',fontsize=size)
        # ax[i+1].yaxis.set_minor_locator(MultipleLocator(ytick))
        ax[i+1].set_yticks([],[])
        # ax[i+1].set_ylabel('Base Pairs in Mb',fontsize=size)
        ax[i+1].invert_yaxis()
        ax[i+1].minorticks_on()
        ax[i+1].tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
        ax[i+1].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
        for spine in ax[i+1].spines.values():
                spine.set_linewidth(wth)
        ax[i+1].text(0.05,0.95,'Pat',transform=ax[i+1].transAxes,verticalalignment='top',horizontalalignment='left',fontsize=size)
        ax[i+1].text(0.95,0.05,'Mat',transform=ax[i+1].transAxes,verticalalignment='bottom',horizontalalignment='right',fontsize=size)
        ax[i+1].text(0.5,1.05,'$10^{%s}\\tau$'%int(np.log10(t_show[i])),transform=ax[i+1].transAxes,verticalalignment='bottom',horizontalalignment='center',fontsize=size)

cbar_ax=fig.add_axes([0.9,0.25,0.02,0.59])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.set_label('$P$',fontsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0.1,0.9,0.9])
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')

