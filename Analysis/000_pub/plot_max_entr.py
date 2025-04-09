"""Import Modules"""
import numpy as np
from sklearn.metrics import r2_score
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
from dist_prob import dist_prob
from cpmt_sgnl import cpmt_sgnl
from insul_score import insul_score

"""Set Arguments"""
npyfile='../cont_prob/cont_prob.npy'
reffile=[[],[]]
reffile[0]+=['../cont_prob/ref/zp.dat']
reffile[0]+=['../cont_prob/ref/zm.dat']
reffile[1]+=['../cont_prob/ref/esc.dat']
reffile[1]+=['../cont_prob/ref/esc.dat']
figname='max_entr'
cl=[['ZP','ZM'],['ESC','ESC']]
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
nfr=56
nat=600
binsize=1000000
idx=[0,-1]

"""Read Data, Make Some Modifications in Order to Make Plot More Convenient"""
P_exp=np.zeros((nsw,ncl,nat,nat))
P_nor1_exp=np.zeros((nsw,ncl,nat,nat))
P_nor2_exp=np.zeros((nsw,ncl,nat,nat))
DP_exp=np.zeros((nsw,ncl,nat))
CS_exp=np.zeros((nsw,ncl,nat))
IS_exp=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile[j][i])
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            P_exp[i,j,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
        DP_exp[i,j]=dist_prob(P_exp[i,j])
        CS_exp[i,j]=cpmt_sgnl(P_exp[i,j])[0]
        P_nor10=cpmt_sgnl(P_exp[i,j])[1]
        P_nor10[P_nor10==0]=np.min(P_nor10[P_nor10>0])/2
        P_nor10=np.log2(P_nor10)
        mx=max(np.max(P_nor10),-np.min(P_nor10))
        P_nor10[P_nor10==np.max(P_nor10)]=mx
        P_nor10[P_nor10==np.min(P_nor10)]=-mx
        P_nor1_exp[i,j]=P_nor10
        P_nor2_exp[i,j]=np.ma.corrcoef(P_nor1_exp[i,j])
        IS_exp[i,j]=insul_score(P_exp[i,j],bin_size=10**5)[0]
        
P_all_sim=np.load(npyfile)
P_sim=np.zeros((nsw,ncl,nat,nat))
P_nor1_sim=np.zeros((nsw,ncl,nat,nat))
P_nor2_sim=np.zeros((nsw,ncl,nat,nat))
DP_sim=np.zeros((nsw,ncl,nat))
CS_sim=np.zeros((nsw,ncl,nat))
IS_sim=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    P_temp=P_all_sim[i]
    for j in range(ncl):
        P_sim[i,j]=P_temp[idx[j]]
        DP_sim[i,j]=dist_prob(P_sim[i,j])
        CS_sim[i,j]=cpmt_sgnl(P_sim[i,j])[0]
        if np.corrcoef(CS_sim[i,j],CS_exp[i,j])[0,1]<0:
            CS_sim[i,j]*=-1
        P_nor10=cpmt_sgnl(P_sim[i,j])[1]
        P_nor10[P_nor10==0]=np.min(P_nor10[P_nor10>0])/2
        P_nor10=np.log2(P_nor10)
        mx=max(np.max(P_nor10),-np.min(P_nor10))
        P_nor10[P_nor10==np.max(P_nor10)]=mx
        P_nor10[P_nor10==np.min(P_nor10)]=-mx
        P_nor1_sim[i,j]=P_nor10
        P_nor2_sim[i,j]=np.ma.corrcoef(P_nor1_sim[i,j])
        IS_sim[i,j]=insul_score(P_sim[i,j],bin_size=10**5)[0]

P_exp_=np.zeros(np.shape(P_exp))
P_exp_+=P_exp
P_sim_=np.zeros(np.shape(P_sim))
P_sim_+=P_sim
mat=np.zeros((nsw,ncl,2,nat,nat))
for i in range(nsw):
    for j in range(ncl):
        P_exp_[i,j][P_exp_[i,j]==0]=np.min(P_exp_[i,j][P_exp_[i,j]>0])/2
        P_sim_[i,j][P_sim_[i,j]==0]=np.min(P_sim_[i,j][P_sim_[i,j]>0])/2
        for k in range(nat):
            for l in range(k):
                mat[i,j,0,k,l]=np.log10(P_exp_[i,j,k,l])
                mat[i,j,0,l,k]=np.log10(P_sim_[i,j,l,k])
                mat[i,j,1,k,l]=P_nor1_exp[i,j,k,l]
                mat[i,j,1,l,k]=P_nor1_sim[i,j,l,k]
        # mx=max(np.max(mat[i,j,1]),-np.min(mat[i,j,1]))
        mx=5
        mat[i,j,1][mat[i,j,1]==np.max(mat[i,j,1])]=mx
        mat[i,j,1][mat[i,j,1]==np.min(mat[i,j,1])]=-mx
        mat[i,j,1][mat[i,j,1]>mx]=mx
        mat[i,j,1][mat[i,j,1]<-mx]=-mx

vec=[]
for i in range(nsw):
    vec.append([])
    for j in range(ncl):
        vec[-1].append([])
        vec[-1][-1].append([P_exp[i,j].reshape(-1),P_sim[i,j].reshape(-1)])
        x=P_nor1_exp[i,j].reshape(-1)
        y=P_nor1_sim[i,j].reshape(-1)
        mx=np.max(x)
        x[x==mx]=y[x==mx]
        x[x==-mx]=y[x==-mx]
        mx=np.max(y)
        y[y==mx]=x[y==mx]
        y[y==-mx]=x[y==-mx]
        vec[-1][-1].append([x,y])
        vec[-1][-1].append([DP_exp[i,j],DP_sim[i,j]])
        vec[-1][-1].append([CS_exp[i,j],CS_sim[i,j]])
        # vec[-1][-1].append([IS_exp[i,j,:int(nat/3)],IS_sim[i,j,:int(nat/3)]])
        # vec[-1][-1].append([IS_exp[i,j,int(nat/3):int(2*nat/3)],IS_sim[i,j,int(nat/3):int(2*nat/3)]])
        # vec[-1][-1].append([IS_exp[i,j,int(2*nat/3):],IS_sim[i,j,int(2*nat/3):]])
        vec[-1][-1].append([IS_exp[i,j],IS_sim[i,j]])

r2=np.zeros((nsw,ncl,5))
for i in range(nsw):
    for j in range(ncl):
        for k in range(4):
            r2[i,j,k]=r2_score(vec[i][j][k][1],vec[i][j][k][0])
        r2[i,j,4]=r2_score(IS_sim[i,j],IS_exp[i,j])

"""Plot"""
# for i in range(nsw):
#     for j in range(ncl):
#         fig,ax=plt.subplots(3,4,figsize=(25,15))
#         wth=3
#         size=30
#         lenmaj=15
#         lenmin=8
#         lenbar=8
#         xtick=100
#         ytick=100
#         ann=['(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)','(K)','(L)']
#         for axs,anns in zip(ax.flat,ann):
#             axs.annotate(anns,xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
# 
#         idx1=[(0,0),(1,0)]
#         cmap=[]
#         norm=[]
#         color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
#         nodes=[0/2,1/2,2/2]
#         cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
#         norm+=[colors.LogNorm(vmin=1e-3,vmax=1)]
#         color=[(0.0,0.0,1.0),(1.0,1.0,1.0),(1.0,0.0,0.0)]
#         nodes=[0/2,1/2,2/2]
#         cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
#         norm+=[colors.Normalize()]
#         clabels=['$P$','$S_\\mathrm{C}$']
#         for k in range(len(idx1)):
#             I=[0,-1]
#             img=ax[idx1[k]].imshow(mat[i,j,k],cmap=cmap[k],norm=norm[k])
#             x_positions=[0.5,200.5,400.5,599.5]
#             y_positions=[0.5,200.5,400.5,599.5]
#             x_labels=[0,20,40,60]
#             y_labels=[0,20,40,60]
#             ax[idx1[k]].xaxis.set_minor_locator(MultipleLocator(xtick))
#             ax[idx1[k]].set_xticks(x_positions,x_labels)
#             ax[idx1[k]].set_xlabel('Loci (Mb)',fontsize=size)
#             ax[idx1[k]].yaxis.set_minor_locator(MultipleLocator(ytick))
#             ax[idx1[k]].set_yticks(y_positions,y_labels)
#             ax[idx1[k]].set_ylabel('Loci (Mb)',fontsize=size)
#             ax[idx1[k]].invert_yaxis()
#             ax[idx1[k]].minorticks_on()
#             ax[idx1[k]].tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
#             ax[idx1[k]].tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
#             for spine in ax[idx1[k]].spines.values():
#                 spine.set_linewidth(wth)
#             ax[idx1[k]].text(0.05,0.95,'Exp',transform=ax[idx1[k]].transAxes,verticalalignment='top',horizontalalignment='left',fontsize=size)
#             ax[idx1[k]].text(0.95,0.05,'Sim',transform=ax[idx1[k]].transAxes,verticalalignment='bottom',horizontalalignment='right',fontsize=size)
# 
#             divider=make_axes_locatable(ax[idx1[k]])
#             cax=divider.append_axes("right",size="5%",pad=0.1)
#             cbar=plt.colorbar(img,cax=cax)
#             # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
#             # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
#             cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
#             cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
#             cbar.set_label(clabels[k],fontsize=size)
#             cbar.outline.set_linewidth(wth)
#         
#         idx2=[(0,2),(1,2),(2,0),(2,1),(2,2)]
#         xlabels=['$\\xi$ (Mb)','Loci (Mb)','Loci (Mb)','Loci (Mb)','Loci (Mb)']
#         ylabels=['$P_\\xi$','$S_\\mathrm{c}$','$S_\\mathrm{i}$','$S_\\mathrm{i}$','$S_\\mathrm{i}$']
#         x=[np.arange(0,60,0.1),np.arange(0,60,0.1),np.arange(0,20,0.1),np.arange(20,40,0.1),np.arange(40,60,0.1)]
#         for k in range(len(idx2)):
#             lines=[]
#             lines+=ax[idx2[k]].plot(x[k],vec[i][j][k+2][0],linewidth=wth,label='Exp')
#             lines+=ax[idx2[k]].plot(x[k],vec[i][j][k+2][1],linewidth=wth,label='Sim')
#             ax[idx2[k]].minorticks_on()
#             ax[idx2[k]].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#             ax[idx2[k]].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#             if k==0:
#                 ax[idx2[k]].set_xscale('log')
#                 ax[idx2[k]].set_yscale('log')
#             # ax[idx2[k]].xaxis.set_major_locator(MultipleLocator(xtick))
#             ax[idx2[k]].xaxis.set_minor_locator(AutoMinorLocator(2))
#             # ax[idx2[k]].yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
#             ax[idx2[k]].yaxis.set_minor_locator(AutoMinorLocator(2))
#             ax[idx2[k]].spines['bottom'].set_linewidth(wth)
#             ax[idx2[k]].spines['top'].set_linewidth(wth)
#             ax[idx2[k]].spines['left'].set_linewidth(wth)
#             ax[idx2[k]].spines['right'].set_linewidth(wth)
#             ax[idx2[k]].set_xlabel(xlabels[k],fontsize=size)
#             ax[idx2[k]].set_ylabel(ylabels[k],fontsize=size)
#             legend=ax[idx2[k]].legend(loc='best',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
#             for text,line in zip(legend.get_texts(),lines):
#                 text.set_color(line.get_color())
#             
#         idx3=[(0,1),(1,1),(0,3),(1,3),(2,3)]
#         xlabels=['$P^\\mathrm{exp}$','$S_\\mathrm{C}^\\mathrm{exp}$','$P_\\xi^\\mathrm{exp}$','$S_\\mathrm{c}^\\mathrm{exp}$','$S_\\mathrm{i}^\\mathrm{exp}$']
#         ylabels=['$P^\\mathrm{sim}$','$S_\\mathrm{C}^\\mathrm{sim}$','$P_\\xi^\\mathrm{sim}$','$S_\\mathrm{c}^\\mathrm{sim}$','$S_\\mathrm{i}^\\mathrm{sim}$']
#         for k in range(len(idx3)):
#             if k==4:
#                 ax[idx3[k]].plot(IS_exp[i,j],IS_sim[i,j],'.',color='c',linewidth=wth)
#                 ax[idx3[k]].plot([np.min(IS_exp[i,j]),np.max(IS_exp[i,j])],
#                                  [np.min(IS_exp[i,j]),np.max(IS_exp[i,j])],
#                                  color='k',linewidth=wth,label='$R^2=%.2f$'%r2[i,j,k])
#             else:
#                 ax[idx3[k]].plot(vec[i][j][k][0],vec[i][j][k][1],'.',color='c',linewidth=wth)
#                 ax[idx3[k]].plot([np.min(vec[i][j][k][0]),np.max(vec[i][j][k][0])],
#                                  [np.min(vec[i][j][k][0]),np.max(vec[i][j][k][0])],
#                                  color='k',linewidth=wth,label='$R^2=%.2f$'%r2[i,j,k])
#             ax[idx3[k]].minorticks_on()
#             ax[idx3[k]].tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
#             ax[idx3[k]].tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
#             if k in [0,2]:
#                 ax[idx3[k]].set_xscale('log')
#                 ax[idx3[k]].set_yscale('log')
#             # ax[idx3[k]].xaxis.set_major_locator(MultipleLocator(xtick))
#             ax[idx3[k]].xaxis.set_minor_locator(AutoMinorLocator(2))
#             # ax[idx3[k]].yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
#             ax[idx3[k]].yaxis.set_minor_locator(AutoMinorLocator(2))
#             ax[idx3[k]].spines['bottom'].set_linewidth(wth)
#             ax[idx3[k]].spines['top'].set_linewidth(wth)
#             ax[idx3[k]].spines['left'].set_linewidth(wth)
#             ax[idx3[k]].spines['right'].set_linewidth(wth)
#             ax[idx3[k]].set_xlabel(xlabels[k],fontsize=size)
#             ax[idx3[k]].set_ylabel(ylabels[k],fontsize=size)
#             ax[idx3[k]].legend(loc='upper left',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
#         
#         plt.tight_layout()
#         plt.savefig(f'{figname}_{cl_[j][i]}.png',format='png')
#         plt.savefig(f'{figname}_{cl_[j][i]}.pdf',format='pdf')

for i in range(nsw):
    for j in range(ncl):
        fig=plt.figure(figsize=(25,17))
        wth=3
        size=30
        lenmaj=15
        lenmin=8
        lenbar=8
        xtick=100
        ytick=100
        ctick=5
        ann=[['(A)','(B)','(C)','(D)'],['(E)','(F)','(G)','(H)'],['(I)','(I)','(I)','(J)']]

        subp1=[(0,0),(1,0)]
        cmap=[]
        norm=[]
        color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
        nodes=[0/2,1/2,2/2]
        cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
        norm+=[colors.Normalize(vmin=-3,vmax=0)]
        color=[(0.0,0.0,1.0),(1.0,1.0,1.0),(1.0,0.0,0.0)]
        nodes=[0/2,1/2,2/2]
        cmap+=[LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))]
        norm+=[colors.Normalize()]
        ctitle=['lg$P$','$S_\\mathrm{C}$']
        for k in range(len(subp1)):
            ax=plt.subplot2grid((3,4),subp1[k],colspan=1,rowspan=1)
            ax.annotate(ann[subp1[k][0]][subp1[k][1]],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
            img=ax.imshow(mat[i,j,k],cmap=cmap[k],norm=norm[k])
            x_positions=[0.5,200.5,400.5,599.5]
            y_positions=[0.5,200.5,400.5,599.5]
            x_labels=[0,20,40,60]
            y_labels=[0,20,40,60]
            ax.xaxis.set_minor_locator(MultipleLocator(xtick))
            ax.set_xticks(x_positions,x_labels)
            ax.set_xlabel('Loci (Mb)',fontsize=size)
            ax.yaxis.set_minor_locator(MultipleLocator(ytick))
            ax.set_yticks(y_positions,y_labels)
            ax.set_ylabel('Loci (Mb)',fontsize=size)
            ax.invert_yaxis()
            ax.minorticks_on()
            ax.tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
            ax.tick_params(axis='both',which='minor',direction='out',width=wth,length=lenmin,labelsize=size)
            for spine in ax.spines.values():
                spine.set_linewidth(wth)
            ax.text(0.05,0.95,'Exp',transform=ax.transAxes,verticalalignment='top',horizontalalignment='left',fontsize=size)
            ax.text(0.95,0.05,'Sim',transform=ax.transAxes,verticalalignment='bottom',horizontalalignment='right',fontsize=size)

            divider=make_axes_locatable(ax)
            cax=divider.append_axes("right",size="5%",pad=0.1)
            cbar=plt.colorbar(img,cax=cax)
            # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
            # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
            cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
            cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
            cbar.ax.set_title(ctitle[k],fontsize=size,pad=20)
            if k==1:
                cbar.ax.yaxis.set_major_locator(MultipleLocator(ctick))
            cbar.outline.set_linewidth(wth)

        subp2=[(0,2),(1,2),(2,0)]
        colspan=[1,1,3]
        xy=[(-0.1,1.05),(-0.1,1.05),(-0.02,1.05)]
        xlabels=['$\\xi$ (Mb)','Loci (Mb)','Loci (Mb)']
        ylabels=['$P_\\xi$','$S_\\mathrm{c}$','$S_\\mathrm{i}$']
        x=np.arange(0,60,0.1)
        for k in range(len(subp2)):
            ax=plt.subplot2grid((3,4),subp2[k],colspan=colspan[k],rowspan=1)
            ax.annotate(ann[subp2[k][0]][subp2[k][1]],xy=xy[k],xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
            lines=[]
            lines+=ax.plot(x,vec[i][j][k+2][0],linewidth=wth,label='Exp')
            lines+=ax.plot(x,vec[i][j][k+2][1],linewidth=wth,label='Sim')
            ax.minorticks_on()
            ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
            if k==0:
                ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
                ax.set_xscale('log')
                ax.set_yscale('log')
            else:
                ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
                # ax.xaxis.set_major_locator(MultipleLocator(xtick))
                ax.xaxis.set_minor_locator(AutoMinorLocator(2))
                # ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
                ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.spines['bottom'].set_linewidth(wth)
            ax.spines['top'].set_linewidth(wth)
            ax.spines['left'].set_linewidth(wth)
            ax.spines['right'].set_linewidth(wth)
            ax.set_xlabel(xlabels[k],fontsize=size)
            ax.set_ylabel(ylabels[k],fontsize=size)
            legend=ax.legend(loc='best',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
            for text,line in zip(legend.get_texts(),lines):
                text.set_color(line.get_color())

        subp3=[(0,1),(1,1),(0,3),(1,3),(2,3)]
        xlabels=['$P^\\mathrm{exp}$','$S_\\mathrm{C}^\\mathrm{exp}$','$P_\\xi^\\mathrm{exp}$',
                 '$S_\\mathrm{c}^\\mathrm{exp}$','$S_\\mathrm{i}^\\mathrm{exp}$']
        ylabels=['$P^\\mathrm{sim}$','$S_\\mathrm{C}^\\mathrm{sim}$','$P_\\xi^\\mathrm{sim}$',
                 '$S_\\mathrm{c}^\\mathrm{sim}$','$S_\\mathrm{i}^\\mathrm{sim}$']
        for k in range(len(subp3)):
            ax=plt.subplot2grid((3,4),subp3[k],colspan=1,rowspan=1)
            ax.annotate(ann[subp3[k][0]][subp3[k][1]],xy=(-0.1,1.05),xycoords='axes fraction',fontsize=size,fontweight='bold',va='bottom',ha='right')
            ax.plot(vec[i][j][k][0],vec[i][j][k][1],'.',color='c',linewidth=wth)
            ax.plot([np.min(vec[i][j][k][0]),np.max(vec[i][j][k][0])],[np.min(vec[i][j][k][0]),np.max(vec[i][j][k][0])],color='k',linewidth=wth,label='$R^2=%.2f$'%r2[i,j,k])
            ax.minorticks_on()
            ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
            if k in [0,2]:
                ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=0,labelsize=size)
                ax.set_xscale('log')
                ax.set_yscale('log')
            else:
                ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
                # ax.xaxis.set_major_locator(MultipleLocator(xtick))
                ax.xaxis.set_minor_locator(AutoMinorLocator(2))
                # ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
                ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.spines['bottom'].set_linewidth(wth)
            ax.spines['top'].set_linewidth(wth)
            ax.spines['left'].set_linewidth(wth)
            ax.spines['right'].set_linewidth(wth)
            ax.set_xlabel(xlabels[k],fontsize=size)
            ax.set_ylabel(ylabels[k],fontsize=size)
            ax.legend(loc='upper left',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)

        plt.tight_layout()
        plt.savefig(f'{figname}_{cl_[j][i]}.png',format='png')
        plt.savefig(f'{figname}_{cl_[j][i]}.pdf',format='pdf')

