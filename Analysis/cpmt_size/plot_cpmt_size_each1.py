"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from tot_len import tot_len

"""Set Aruments"""
file='cpmt_size_each1.npz'
figname=file[:-4]
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nfr=56
nhist=50
sgm=4.0
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]

"""Read Data"""
data=np.load(file)
Rga0=data['Rga0']/sgm
P_Rga0=data['P_Rga0']*sgm
Rgb0=data['Rgb0']/sgm
P_Rgb0=data['P_Rgb0']*sgm
Rga0_mean=data['Rga0_mean']/sgm
Rga0_std=data['Rga0_std']/sgm
Rgb0_mean=data['Rgb0_mean']/sgm
Rgb0_std=data['Rgb0_std']/sgm
Rga=data['Rga']/sgm
P_Rga=data['P_Rga']*sgm
Rgb=data['Rgb']/sgm
P_Rgb=data['P_Rgb']*sgm
Rga_mean=data['Rga_mean']/sgm
Rga_std=data['Rga_std']/sgm
Rgb_mean=data['Rgb_mean']/sgm
Rgb_std=data['Rgb_std']/sgm

na=np.shape(Rga0_mean)[2]
nb=np.shape(Rgb0_mean)[2]
la=[]
for i in range(na):
    la+=[np.sum(Rga0_mean[0,0,i]!=0)+1]
lb=[]
for i in range(nb):
    lb+=[np.sum(Rgb0_mean[0,0,i]!=0)+1]
la=np.array(la)
lb=np.array(lb)

xia_mean=(Rga_mean-Rga_mean[:,-10:-9])/(Rga_mean[:,1:2]-Rga_mean[:,-10:-9])
xia_std=(Rga_std-Rga_std[:,-10:-9])/(Rga_std[:,1:2]-Rga_std[:,-10:-9])
xib_mean=(Rgb_mean-Rgb_mean[:,-10:-9])/(Rgb_mean[:,1:2]-Rgb_mean[:,-10:-9])
xib_std=(Rgb_std-Rgb_std[:,-10:-9])/(Rgb_std[:,1:2]-Rgb_std[:,-10:-9])
Xia_mean=np.sum(xia_mean[:,:,1:]*la[np.newaxis,np.newaxis,1:],axis=-1)/np.sum(la[1:])
Xia_std=np.sum(xia_std[:,:,1:]*la[np.newaxis,np.newaxis,1:],axis=-1)/np.sum(la[1:])
Xib_mean=np.sum(xib_mean[:,:,1:]*lb[np.newaxis,np.newaxis,1:],axis=-1)/np.sum(lb[1:])
Xib_std=np.sum(xib_std[:,:,1:]*lb[np.newaxis,np.newaxis,1:],axis=-1)/np.sum(lb[1:])

"""Plot"""
fig,ax=plt.subplots(4,max(na,nb),figsize=(25,18))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.LogNorm(vmin=min(t[1:-9]),vmax=max(t[1:-9]))

for i in range(nsw):
    for j in range(na):
        for k in range(1,nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[2*i,j].plot(Rga[i,k,j],P_Rga[i,k,j],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[2*i,j].plot(Rga0[i,k,j],P_Rga0[i,k,j],'-',color='k',linewidth=2*wth)
            ax[2*i,j].plot(Rga0[i,k,j],P_Rga0[i,k,j],'-',color=color0[i][k],linewidth=wth,label=label0[i][k])
        ax[2*i,j].minorticks_on()
        ax[2*i,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[2*i,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax[2*i,j].xaxis.set_major_locator(MultipleLocator(xtick))
        # ax[2*i,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax[2*i,j].yaxis.set_major_locator(MultipleLocator(ytick))
        # ax[2*i,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[2*i,j].spines['bottom'].set_linewidth(wth)
        ax[2*i,j].spines['top'].set_linewidth(wth)
        ax[2*i,j].spines['left'].set_linewidth(wth)
        ax[2*i,j].spines['right'].set_linewidth(wth)
        ax[2*i,j].legend(loc='upper right',fontsize=0.8*size)
        ax[2*i,j].set_xlabel('$R_\\mathrm{g}~(\\sigma)$ of A%s'%(j+1),fontsize=size)
        if j==0:
            ax[2*i,j].set_ylabel('$P~(\\sigma^{-1})$ in %s'%sw[i],fontsize=size)

    for j in range(nb):
        for k in range(1,nfr-9):
            normc=norm(t[k])
            rgb=cmap(normc)
            ax[2*i+1,j].plot(Rgb[i,k,j],P_Rgb[i,k,j],'-',color=rgb,linewidth=wth)
        for k in range(ncl):
            ax[2*i+1,j].plot(Rgb0[i,k,j],P_Rgb0[i,k,j],'-',color='k',linewidth=2*wth)
            ax[2*i+1,j].plot(Rgb0[i,k,j],P_Rgb0[i,k,j],'-',color=color0[i][k],linewidth=wth,label=label0[i][k])
        ax[2*i+1,j].minorticks_on()
        ax[2*i+1,j].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[2*i+1,j].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax[2*i+1,j].xaxis.set_major_locator(MultipleLocator(xtick))
        ax[2*i+1,j].xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax[2*i+1,j].yaxis.set_major_locator(MultipleLocator(ytick))
        ax[2*i+1,j].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[2*i+1,j].spines['bottom'].set_linewidth(wth)
        ax[2*i+1,j].spines['top'].set_linewidth(wth)
        ax[2*i+1,j].spines['left'].set_linewidth(wth)
        ax[2*i+1,j].spines['right'].set_linewidth(wth)
        ax[2*i+1,j].legend(loc='upper right',fontsize=0.8*size)
        ax[2*i+1,j].set_xlabel('$R_\\mathrm{g}~(\\sigma)$ of B%s'%(j+1),fontsize=size)
        if j==0:
            ax[2*i+1,j].set_ylabel('$P~(\\sigma^{-1})$ in %s'%sw[i],fontsize=size)

fig.subplots_adjust(right=0.85)
cbar_ax=fig.add_axes([0.88,0.15,0.03,0.7])
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
cbar=plt.colorbar(sm,cax=cbar_ax)
cbar.set_label('$t~(\\tau)$',fontsize=size)
# cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
# cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
cbar.outline.set_linewidth(wth)

plt.tight_layout(rect=[0,0,0.85,1])
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()

fig,ax=plt.subplots(2,max(na,nb),figsize=(25,10))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,1.0,0.6)],[(0.8,0.6,1.0),(0.8,1.0,0.6)]]
label0=[['ZP','ESC'],['ZM','ESC']]
color=[[(0,0.5,0),(0.5,0,0)],[(0,1,0),(1,0,0)]]

for i in range(na):
    ax2=ax[0,i].twinx()
    for j in range(nsw):
        # ax[0,i].fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
        # for k in range(ncl):
        #     if (i,k)!=(0,1):
        #         ax[0,i].plot([min(t[1:-9]),max[0,i](t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
        ax2.plot(t[1:-9],Rga_mean[j,1:-9,i],'-',color=color[j][0],linewidth=wth,label=sw[j])
    ax[0,i].minorticks_on()
    ax[0,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax[0,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax[0,i].set_xscale('log')
    # ax[0,i].xax[0,i]is.set_major_locator(MultipleLocator(xtick))
    ax[0,i].xaxis.set_minor_locator(AutoMinorLocator(2))
    # ax[0,i].yax[0,i]is.set_major_locator(MultipleLocator(ytick))
    ax[0,i].yaxis.set_minor_locator(AutoMinorLocator(2))
    ax[0,i].spines['bottom'].set_linewidth(wth)
    ax[0,i].spines['top'].set_linewidth(wth)
    ax[0,i].spines['left'].set_linewidth(wth)
    ax[0,i].spines['right'].set_linewidth(wth)
    # ymax[0,i]=max[0,i](np.max[0,i](r_mean[:,1:-9,:]),np.max[0,i](r0_mean))
    # ymin=min(np.min(r_mean[:,1:-9,:]),np.min(r0_mean))
    # ax[0,i].set_ylim(ymin-0.1*(ymax[0,i]-ymin),ymax[0,i]+0.1*(ymax[0,i]-ymin))
    # ax[0,i].legend(loc='lower right',fontsize=0.8*size)
    ax[0,i].set_ylim(ax2.get_ylim())
    ax[0,i].set_xlabel(f'Time $(\\tau)$',fontsize=size)
    ax[0,i].set_ylabel('$\\langle{}R_\\mathrm{g}\\rangle{} ~(\\sigma)$ of A%s'%(i+1),fontsize=size)
    ax2.set_yticks([])
    # ax2.set_ylim(ax[0,i].get_ylim())
    ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
    ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
    ax2.legend(loc='upper right',fontsize=0.6*size)

for i in range(nb):
    ax2=ax[1,i].twinx()
    for j in range(nsw):
        # ax[1,i].fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
        # for k in range(ncl):
        #     if (i,k)!=(0,1):
        #         ax[1,i].plot([min(t[1:-9]),max[1,i](t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
        ax2.plot(t[1:-9],Rgb_mean[j,1:-9,i],'-',color=color[j][1],linewidth=wth,label=sw[j])
        ax[1,i].minorticks_on()
        ax[1,i].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax[1,i].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax[1,i].set_xscale('log')
        # ax[1,i].xax[1,i]is.set_major_locator(MultipleLocator(xtick))
        ax[1,i].xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax[1,i].yax[1,i]is.set_major_locator(MultipleLocator(ytick))
        ax[1,i].yaxis.set_minor_locator(AutoMinorLocator(2))
        ax[1,i].spines['bottom'].set_linewidth(wth)
        ax[1,i].spines['top'].set_linewidth(wth)
        ax[1,i].spines['left'].set_linewidth(wth)
        ax[1,i].spines['right'].set_linewidth(wth)
        # ymax[1,i]=max[1,i](np.max[1,i](r_mean[:,1:-9,:]),np.max[1,i](r0_mean))
        # ymin=min(np.min(r_mean[:,1:-9,:]),np.min(r0_mean))
        # ax[1,i].set_ylim(ymin-0.1*(ymax[1,i]-ymin),ymax[1,i]+0.1*(ymax[1,i]-ymin))
        # ax[1,i].legend(loc='lower right',fontsize=0.8*size)
        ax[1,i].set_ylim(ax2.get_ylim())
        ax[1,i].set_xlabel(f'Time $(\\tau)$',fontsize=size)
        ax[1,i].set_ylabel('$\\langle{}R_\\mathrm{g}\\rangle{} ~(\\sigma)$ of B%s'%(i+1),fontsize=size)
        ax2.set_yticks([])
        # ax2.set_ylim(ax[1,i].get_ylim())
        ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
        ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
        ax2.legend(loc='upper right',fontsize=0.6*size)

plt.tight_layout()
plt.savefig(f'{figname}_mean_each.png', format='png')
plt.savefig(f'{figname}_mean_each.pdf', format='pdf')
# plt.show()

# fig,ax=plt.subplots(1,1,figsize=(9,5))
# wth=3
# size=30
# lenmaj=15
# lenmin=8
# xtick=1
# ytick=0.02
# color0=[[(1.0,0.8,0.6),(0.8,0.6,1.0)],[(0.8,1.0,0.6),(0.8,1.0,0.6)]]
# label0=[['ZP','ZM'],['ESC','ESC']]
# color=[[(0,0.5,0),(0.5,0,0)],[(0,1,0),(1,0,0)]]
#
# ax2=ax.twinx()
# for i in range(nsw):
#         # ax.fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
#         # for k in range(ncl):
#         #     if (i,k)!=(0,1):
#         #         ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
#         ax2.plot(t[1:-9],RGa_mean[i,1:-9],'-',color=color[i][0],linewidth=wth,label='A in '+sw[i])
#         # ax.fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
#         # for k in range(ncl):
#         #     if (i,k)!=(0,1):
#         #         ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
#         ax2.plot(t[1:-9],RGb_mean[i,1:-9],'-',color=color[i][1],linewidth=wth,label='B in '+sw[i])
# ax.minorticks_on()
# ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
# ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
# ax.set_xscale('log')
# # ax.xaxis.set_major_locator(MultipleLocator(xtick))
# ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# # ax.yaxis.set_major_locator(MultipleLocator(ytick))
# ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax.spines['bottom'].set_linewidth(wth)
# ax.spines['top'].set_linewidth(wth)
# ax.spines['left'].set_linewidth(wth)
# ax.spines['right'].set_linewidth(wth)
# # ymax=max(np.max(r_mean[:,1:-9,:]),np.max(r0_mean))
# # ymin=min(np.min(r_mean[:,1:-9,:]),np.min(r0_mean))
# # ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
# # ax.legend(loc='lower right',fontsize=0.8*size)
# ax.set_ylim(ax2.get_ylim())
# ax.set_xlabel(f'Time $(\\tau)$',fontsize=size)
# ax.set_ylabel('$\\langle{}R_\\mathrm{g}\\rangle{} ~(\\sigma)$',fontsize=size)
# ax2.set_yticks([])
# # ax2.set_ylim(ax.get_ylim())
# ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
# ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
# ax2.legend(loc='upper right',fontsize=0.6*size)
#
# plt.tight_layout()
# plt.savefig(f'{figname}_mean.png', format='png')
# plt.savefig(f'{figname}_mean.pdf', format='pdf')
# # plt.show()

fig,ax=plt.subplots(1,1,figsize=(9,5))
wth=3
size=30
lenmaj=15
lenmin=8
xtick=1
ytick=0.02
color0=[[(1.0,0.8,0.6),(0.8,0.6,1.0)],[(0.8,1.0,0.6),(0.8,1.0,0.6)]]
label0=[['ZP','ZM'],['ESC','ESC']]
color=[[(0,0.5,0),(0.5,0,0)],[(0,1,0),(1,0,0)]]

ax2=ax.twinx()
for i in range(nsw):
        # ax.fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
        # for k in range(ncl):
        #     if (i,k)!=(0,1):
        #         ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
        ax2.plot(t[1:-9],Xia_mean[i,1:-9],'-',color=color[i][0],linewidth=wth,label='A in '+sw[i])
        # ax.fill_between(t[1:-9],r_mean[i,1:-9,j]-r_std[i,1:-9,j]/10,r_mean[i,1:-9,j]+r_std[i,1:-9,j]/10,color=color[i][j],alpha=0.1)
        # for k in range(ncl):
        #     if (i,k)!=(0,1):
        #         ax.plot([min(t[1:-9]),max(t[1:-9])],[r0_mean[i,k,j],r0_mean[i,k,j]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
        ax2.plot(t[1:-9],Xib_mean[i,1:-9],'-',color=color[i][1],linewidth=wth,label='B in '+sw[i])
ax.minorticks_on()
ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
ax.set_xscale('log')
# ax.xaxis.set_major_locator(MultipleLocator(xtick))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax.yaxis.set_major_locator(MultipleLocator(ytick))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.spines['bottom'].set_linewidth(wth)
ax.spines['top'].set_linewidth(wth)
ax.spines['left'].set_linewidth(wth)
ax.spines['right'].set_linewidth(wth)
# ymax=max(np.max(r_mean[:,1:-9,:]),np.max(r0_mean))
# ymin=min(np.min(r_mean[:,1:-9,:]),np.min(r0_mean))
# ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
# ax.legend(loc='lower right',fontsize=0.8*size)
ax.set_ylim(ax2.get_ylim())
ax.set_xlabel(f'Time $(\\tau)$',fontsize=size)
ax.set_ylabel('$\\frac{\\langle{}R_\\mathrm{g}\\rangle-\\langle{}R_\\mathrm{g}^\\mathrm{end}\\rangle}'
              '{\\langle{}R_\\mathrm{g}^\\mathrm{start}\\rangle-\\langle{}R_\\mathrm{g}^\\mathrm{end}\\rangle}$',fontsize=size)
ax2.set_yticks([])
# ax2.set_ylim(ax.get_ylim())
ax2.tick_params(axis='both',which='major',direction='in',width=wth,length=0,labelsize=size)
ax2.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
ax2.legend(loc='upper right',fontsize=0.6*size)

plt.tight_layout()
plt.savefig(f'{figname}_mean1.png', format='png')
plt.savefig(f'{figname}_mean1.pdf', format='pdf')
# plt.show()