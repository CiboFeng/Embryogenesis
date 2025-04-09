"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib.collections import LineCollection
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
from pyw import pyw
from dist_prob import dist_prob
from cpmt_sgnl import cpmt_sgnl
from insul_score import insul_score

"""Set Arguments"""
file1='../clst/rmsd.pyw'
file2='../clst/clst1.pyw'
file4='../clst/clst2.pyw'
file5=['../clst/cont_prob_zp.dat','../clst/cont_prob_zm.dat','../clst/cont_prob_esc.dat']
figname='clst'
cls={'ZP':0.0,'ZM':1.0}
nsw=2
sgm=4
nat=600
nhist=20

"""Read Data"""
q=pyw(file1,dir=True)
t=q[:,:,:,0]
d=q[:,:,:,1]/sgm
dm=np.mean(d,axis=1)

q=pyw(file2,dir=True,mat=False)
Z=[]
for i in range(nsw):
    z=np.zeros((int(np.shape(q)[1]/4),4))
    for j in range(int(len(q[i])/4)):
        for k in range(4):
            z[j,k]=q[i][j*4+k][0]
    z[:,2]/=sgm
    Z+=[z]

idx=[]
DD=np.zeros((nsw,nhist))
PP=np.zeros((nsw,nhist))
P_rdc=np.ones((nsw,nat,nat))
DP_rdc=np.zeros((nsw,nat))
CS_rdc=np.zeros((nsw,nat))
P_nor_rdc=np.ones((nsw,nat,nat))
IS_rdc=np.zeros((nsw,nat))
for i in range(nsw):
    idx+=pyw(file4,f'Index of {list(cls.keys())[i]}:',fmt='int')
    DD[i]=np.array(pyw(file4,f'{list(cls.keys())[i]},')[0])/sgm
    PP[i]=np.array(pyw(file4,f'{list(cls.keys())[i]},')[1])
    q=np.array(pyw(file4,f'Conformations of {list(cls.keys())[i]}:')[0])
    for j in range(nat):
        for k in range(j+2,nat):
            P_rdc[i,j,k]=q[j*nat+k-int(j*(j+5)/2+2)]
            P_rdc[i,k,j]=q[j*nat+k-int(j*(j+5)/2+2)]
    DP_rdc[i]=dist_prob(P_rdc[i])
    CS_rdc[i]=cpmt_sgnl(P_rdc[i])[0]
    P_nor0=cpmt_sgnl(P_rdc[i])[1]
    P_nor0[P_nor0==0]=np.min(P_nor0[P_nor0>0])/2
    P_nor0=np.log2(P_nor0)
    mx=max(np.max(P_nor0),-np.min(P_nor0))
    P_nor0[P_nor0==np.max(P_nor0)]=mx
    P_nor0[P_nor0==np.min(P_nor0)]=-mx
    P_nor_rdc[i]=P_nor0
    IS_rdc[i]=insul_score(P_rdc[i],bin_size=10**5)[0]

P_ful=np.ones((nsw,nat,nat))
DP_ful=np.zeros((nsw,nat))
CS_ful=np.zeros((nsw,nat))
P_nor_ful=np.ones((nsw,nat,nat))
IS_ful=np.zeros((nsw,nat))
for i in range(nsw):
    f=open(file5[i],'r')
    lsf=f.readlines()
    f.close()
    for j in range(len(lsf)):
        ls_f=lsf[j].strip('\n').split()
        P_ful[i,int(float(ls_f[0])-1),int(float(ls_f[1])-1)]=float(ls_f[2])
        P_ful[i,int(float(ls_f[1])-1),int(float(ls_f[0])-1)]=float(ls_f[2])
    DP_ful[i]=dist_prob(P_ful[i])
    CS_ful[i]=cpmt_sgnl(P_ful[i])[0]
    P_nor0=cpmt_sgnl(P_ful[i])[1]
    P_nor0[P_nor0==0]=np.min(P_nor0[P_nor0>0])/2
    P_nor0=np.log2(P_nor0)
    mx=max(np.max(P_nor0),-np.min(P_nor0))
    P_nor0[P_nor0==np.max(P_nor0)]=mx
    P_nor0[P_nor0==np.min(P_nor0)]=-mx
    P_nor_ful[i]=P_nor0
    IS_ful[i]=insul_score(P_ful[i],bin_size=10**5)[0]

vec=[]
for i in range(nsw):
    vec.append([])
    vec[-1].append([P_ful[i].reshape(-1),P_rdc[i].reshape(-1)])
    vec[-1].append([DP_ful[i],DP_rdc[i]])
    x=P_nor_ful[i].reshape(-1)
    y=P_nor_rdc[i].reshape(-1)
    mx=np.max(x)
    x[x==mx]=y[x==mx]
    x[x==-mx]=y[x==-mx]
    mx=np.max(y)
    y[y==mx]=x[y==mx]
    y[y==-mx]=x[y==-mx]
    vec[-1].append([x,y])
    vec[-1].append([CS_ful[i],CS_rdc[i]])
    vec[-1].append([IS_ful[i],IS_rdc[i]])

r2=np.zeros((nsw,5))
for i in range(nsw):
    for j in range(5):
        r2[i,j]=r2_score(vec[i][j][1],vec[i][j][0])

"""Plot"""
# fig,ax=plt.subplots(2,2,figsize=(25,10))
# wth=3
# size=30
# lenmaj=15
# lenmin=8
# lenbar=8
# xtick=100
# ytick=100
# ann=['(A)','(C)','(B)','(D)']
# for axs,anns in zip(ax.flat,ann):
#     axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
#
# ws=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
# nw=len(ws)
# color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
# nodes=[0.00,1/4,2/4,3/4,1.00]
# cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
# norm=colors.Normalize(vmin=min(ws),vmax=max(ws))
# for i in range(nsw):
#     for j in range(nw):
#         normc=norm(ws[j])
#         rgb=cmap(normc)
#         ax[i,0].plot(t[i,j],d[i,j],color=rgb,linewidth=wth,alpha=0.5)
#     ax[i,0].plot(t[i,0],dm[i],color='k',linewidth=wth)
#     ax[i,0].autoscale()
#     ax[i,0].minorticks_on()
#     ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax[i,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i,0].spines['bottom'].set_linewidth(wth)
#     ax[i,0].spines['top'].set_linewidth(wth)
#     ax[i,0].spines['left'].set_linewidth(wth)
#     ax[i,0].spines['right'].set_linewidth(wth)
#     ax[i,0].set_ylabel('$d_\\mathrm{rms}^\\mathrm{%s}~(\\sigma)$'%list(cls.keys())[i],fontsize=size)
#     ax[i,0].set_xscale('log')
#     ax[i,0].set_yscale('log')
#     # ax[i,0].yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))
#     ticks=[2,3,4,5,6]
#     ax[i,0].yaxis.set_major_locator(ticker.FixedLocator(ticks))
#     ax[i,0].yaxis.set_major_formatter(ticker.FixedFormatter([str(tick) for tick in ticks]))
#     # ax[i,0].yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))
#
#     dendrogram(Z[i],ax=ax[i,1])
#     ax[i,1].autoscale()
#     # ax[i,1].minorticks_on()
#     ax[i,1].tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax[i,1].tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     # ax[i,1].xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i,1].yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i,1].spines['bottom'].set_linewidth(wth)
#     ax[i,1].spines['top'].set_linewidth(wth)
#     ax[i,1].spines['left'].set_linewidth(wth)
#     ax[i,1].spines['right'].set_linewidth(wth)
#     ax[i,1].set_xticks([])
#     ax[i,1].set_ylabel('$d_\\mathrm{rms}^\\mathrm{%s}~(\\sigma)$'%list(cls.keys())[i],fontsize=size)
# ax[1,0].set_xlabel('$t~(\\tau)$',fontsize=size)
# ax[1,1].set_xlabel('Structures',fontsize=size)
#
# plt.tight_layout()
# plt.savefig(f'{figname1}.png',format='png')
# plt.savefig(f'{figname1}.pdf',format='pdf')
# # plt.show()
#
# fig,ax=plt.subplots(2,6,figsize=(25,8))
# wth=3
# size=30
# lenmaj=15
# lenmin=8
# lenbar=8
# xtick=100
# ytick=100
# ann=['(E)','(G)','(I)','(K)','(M)','(O)','(F)','(H)','(J)','(L)','(N)','(P)']
# for axs,anns in zip(ax.flat,ann):
#     axs.annotate(anns,xy=(0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
#
# for i in range(nsw):
#     ax[i,0].bar(DD[i],PP[i]/10000,color='skyblue',edgecolor='black',width=DD[i,1]-DD[i,0])
#     ax[i,0].autoscale()
#     ax[i,0].minorticks_on()
#     ax[i,0].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#     ax[i,0].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#     ax[i,0].xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i,0].yaxis.set_minor_locator(AutoMinorLocator(2))
#     ax[i,0].spines['bottom'].set_linewidth(wth)
#     ax[i,0].spines['top'].set_linewidth(wth)
#     ax[i,0].spines['left'].set_linewidth(wth)
#     ax[i,0].spines['right'].set_linewidth(wth)
#     ax[i,0].set_ylabel('$N_\\mathrm{samp}\\times10^4$ (%s)'%list(cls.keys())[i],fontsize=size)
#     if i==1:
#         ax[i,0].set_xlabel('$d_\\mathrm{rms}~(\\sigma)$',fontsize=size)
#
#     xlabels=['$P^\\mathrm{ful}$','$P_\\xi^\\mathrm{ful}$','$S_\\mathrm{C}^\\mathrm{ful}$','$S_\\mathrm{c}^\\mathrm{ful}$','$S_\\mathrm{i}^\\mathrm{ful}$']
#     ylabels=['$P^\\mathrm{rdc,%s}$','$P_\\xi^\\mathrm{rdc,%s}$','$S_\\mathrm{C}^\\mathrm{rdc,%s}$','$S_\\mathrm{c}^\\mathrm{rdc,%s}$','$S_\\mathrm{i}^\\mathrm{rdc,%s}$']
#     for j in range(5):
#         ax[i,j+1].plot(vec[i][j][0],vec[i][j][1],'.',color='c',linewidth=wth)
#         ax[i,j+1].plot([np.min(vec[i][j][0]),np.max(vec[i][j][0])],[np.min(vec[i][j][0]),np.max(vec[i][j][0])],color='k',linewidth=wth,label='$R^2=%.2f$'%r2[i,j])
#         ax[i,j+1].autoscale()
#         ax[i,j+1].minorticks_on()
#         ax[i,j+1].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#         ax[i,j+1].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#         ax[i,j+1].xaxis.set_minor_locator(AutoMinorLocator(2))
#         ax[i,j+1].yaxis.set_minor_locator(AutoMinorLocator(2))
#         ax[i,j+1].spines['bottom'].set_linewidth(wth)
#         ax[i,j+1].spines['top'].set_linewidth(wth)
#         ax[i,j+1].spines['left'].set_linewidth(wth)
#         ax[i,j+1].spines['right'].set_linewidth(wth)
#         ax[i,j+1].set_ylabel(ylabels[j]%list(cls.keys())[i],fontsize=size)
#         ax[i,j+1].legend(loc='upper left',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
#     if i==1:
#         for j in range(5):
#             ax[i,j+1].set_xlabel(xlabels[j],fontsize=size)
#
# plt.tight_layout()
# plt.savefig(f'{figname2}.png',format='png')
# plt.savefig(f'{figname2}.pdf',format='pdf',dpi=10)
# print(plt.rcParams['savefig.dpi'])
# plt.show()

fig=plt.figure(figsize=(25,18))
wth=3
size=30
lenmaj=15
lenmin=8
lenbar=8
xtick=100
ytick=100
ann=[['(A)','(A)','(A)','(C)','(C)','(C)'],
     ['(B)','(B)','(B)','(D)','(D)','(D)'],
     ['(E)','(G)','(I)','(K)','(M)','(O)'],
     ['(F)','(H)','(J)','(L)','(N)','(P)']]
# for axs,anns in zip(ax.flat,ann):
#     axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')

ws=[0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
nw=len(ws)
color=[(0,0,1),(0,1,1),(0,1,0),(1,1,0),(1,0,0)]
nodes=[0.00,1/4,2/4,3/4,1.00]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=min(ws),vmax=max(ws))
for i in range(nsw):
    ax=plt.subplot2grid((4,6),(i,0),colspan=3,rowspan=1)
    ax.annotate(ann[i][0],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    for j in range(nw):
        normc=norm(ws[j])
        rgb=cmap(normc)
        ax.plot(t[i,j],d[i,j],color=rgb,linewidth=wth,alpha=0.5)
    ax.plot(t[i,0],dm[i],color='k',linewidth=wth)
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))
    ticks=[2,3,4,5,6]
    ax.yaxis.set_major_locator(ticker.FixedLocator(ticks))
    ax.yaxis.set_major_formatter(ticker.FixedFormatter([str(tick) for tick in ticks]))
    # ax.yaxis.set_major_locator(LogLocator(base=10.0,numticks=10))
    ax.set_xlabel('$t~(\\tau)$',fontsize=size)
    ax.set_ylabel('$d_\\mathrm{rms}^\\mathrm{%s}~(\\sigma)$'%list(cls.keys())[i],fontsize=size)

    ax=plt.subplot2grid((4,6),(i,3),colspan=3,rowspan=1)
    ax.annotate(ann[i][3],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    dendrogram(Z[i],ax=ax)
    ax.autoscale()
    # ax.minorticks_on()
    ax.tick_params(axis='y',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='y',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    # ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xticks([])
    ax.set_xlabel('Structures',fontsize=size)
    ax.set_ylabel('$d_\\mathrm{rms}^\\mathrm{%s}~(\\sigma)$'%list(cls.keys())[i],fontsize=size)

    ax=plt.subplot2grid((4,6),(i+nsw,0),colspan=1,rowspan=1)
    ax.annotate(ann[i+nsw][0],xy=(0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
    ax.bar(DD[i],PP[i]/10000,color='skyblue',edgecolor='black',width=DD[i,1]-DD[i,0])
    ax.autoscale()
    ax.minorticks_on()
    ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
    ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.spines['bottom'].set_linewidth(wth)
    ax.spines['top'].set_linewidth(wth)
    ax.spines['left'].set_linewidth(wth)
    ax.spines['right'].set_linewidth(wth)
    ax.set_xlabel('$d_\\mathrm{rms}~(\\sigma)$',fontsize=size)
    ax.set_ylabel('$N_\\mathrm{samp}\\times10^4$ (%s)'%list(cls.keys())[i],fontsize=size)

    xlabels=['$P^\\mathrm{ful,%s}$','$P_\\xi^\\mathrm{ful,%s}$','$S_\\mathrm{C}^\\mathrm{ful,%s}$','$S_\\mathrm{c}^\\mathrm{ful,%s}$','$S_\\mathrm{i}^\\mathrm{ful,%s}$']
    ylabels=['$P^\\mathrm{rdc,%s}$','$P_\\xi^\\mathrm{rdc,%s}$','$S_\\mathrm{C}^\\mathrm{rdc,%s}$','$S_\\mathrm{c}^\\mathrm{rdc,%s}$','$S_\\mathrm{i}^\\mathrm{rdc,%s}$']
    for j in range(5):
        ax=plt.subplot2grid((4,6),(i+nsw,j+1),colspan=1,rowspan=1)
        ax.annotate(ann[i+nsw][j+1],xy=(0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
        ax.plot(vec[i][j][0],vec[i][j][1],'.',color='c',linewidth=wth)
        ax.plot([np.min(vec[i][j][0]),np.max(vec[i][j][0])],[np.min(vec[i][j][0]),np.max(vec[i][j][0])],color='k',linewidth=wth,label='$R^2=%.2f$'%r2[i,j])
        ax.autoscale()
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel(xlabels[j]%list(cls.keys())[i],fontsize=size)
        ax.set_ylabel(ylabels[j]%list(cls.keys())[i],fontsize=size)
        ax.legend(loc='upper left',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
# plt.show()