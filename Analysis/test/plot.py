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
reffile1=[[],[]]
reffile1[0]+=[['../max_entr/zp/dist_prob_exp.dat','../max_entr/zp/insul_score_exp.dat','../max_entr/zp/cpmt_sgnl_exp.dat','../max_entr/zp/Po_e.dat']]
reffile1[0]+=[['../max_entr/zm/dist_prob_exp.dat','../max_entr/zm/insul_score_exp.dat','../max_entr/zm/cpmt_sgnl_exp.dat','../max_entr/zm/Po_e.dat']]
reffile1[1]+=[['../max_entr/esc/dist_prob_exp.dat','../max_entr/esc/insul_score_exp.dat','../max_entr/esc/cpmt_sgnl_exp.dat','../max_entr/esc/Po_e.dat']]
reffile1[1]+=[['../max_entr/esc/dist_prob_exp.dat','../max_entr/esc/insul_score_exp.dat','../max_entr/esc/cpmt_sgnl_exp.dat','../max_entr/esc/Po_e.dat']]
reffile2=[[],[]]
reffile2[0]+=['../cont_prob/ref/zp.dat']
reffile2[0]+=['../cont_prob/ref/zm.dat']
reffile2[1]+=['../cont_prob/ref/esc.dat']
reffile2[1]+=['../cont_prob/ref/esc.dat']
figname='test'
cl=[['ZP','ZM'],['ESC','ESC']]
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
nat=600
binsize=1000000

"""Read Data, Make Some Modifications in Order to Make Plot More Convenient"""
P_nor1_he=np.zeros((nsw,ncl,nat,nat))
# P_nor2_he=np.zeros((nsw,ncl,nat,nat))
DP_he=np.zeros((nsw,ncl,nat))
CS_he=np.zeros((nsw,ncl,int(nat/10)))
IS_he=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile1[j][i][0])
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            DP_he[i,j,k]=float(ls_f[1])
        f=open(reffile1[j][i][1])
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            IS_he[i,j,k+5]=float(ls_f[1])
        f=open(reffile1[j][i][2])
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            CS_he[i,j,k]=float(ls_f[0])
        mat=np.loadtxt(reffile1[j][i][3])
        for k in range(int(nat/10)):
            for l in range(int(nat/10)):
                P_nor1_he[i,j,k*10+5,l*10+5]=mat[k,l]
        # P_nor2_he[i,j]=np.ma.corrcoef(P_nor1_he[i,j])

P_me=np.zeros((nsw,ncl,nat,nat))
P_nor1_me=np.zeros((nsw,ncl,nat,nat))
# P_nor2_me=np.zeros((nsw,ncl,nat,nat))
DP_me=np.zeros((nsw,ncl,nat))
CS_me=np.zeros((nsw,ncl,nat))
IS_me=np.zeros((nsw,ncl,nat))
for i in range(nsw):
    for j in range(ncl):
        f=open(reffile2[j][i])
        lsf=f.readlines()
        for k in range(len(lsf)):
            ls_f=lsf[k].strip('\n').split()
            P_me[i,j,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
        DP_me[i,j]=dist_prob(P_me[i,j])
        IS_me[i,j],pos=insul_score(P_me[i,j],bin_size=10**5)
        CS_me[i,j],P_nor10=cpmt_sgnl(P_me[i,j])
        if np.corrcoef(CS_me[i,j][::10],CS_he[i,j])[1,0]<0:
            CS_me[i,j]*=-1
        P_nor10[P_nor10==0]=np.min(P_nor10[P_nor10>0])/2
        P_nor10=np.log2(P_nor10)
        mx=max(np.max(P_nor10),-np.min(P_nor10))
        P_nor10[P_nor10==np.max(P_nor10)]=mx
        P_nor10[P_nor10==np.min(P_nor10)]=-mx
        P_nor1_me[i,j]=P_nor10
        # P_nor2_me[i,j]=np.ma.corrcoef(P_nor1_me[i,j])

"""Plot"""
for i in range(nsw):
    for j in range(ncl):
        fig,ax=plt.subplots(3,4,figsize=(25,15))
        wth=3
        size=30
        lenmaj=15
        lenmin=8
        lenbar=8
        xtick=100
        ytick=100

        color=[(1.0,1.0,1.0),(1.0,0.0,0.0),(0.0,0.0,0.0)]
        nodes=[0/2,1/2,2/2]
        cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
        norm=colors.LogNorm(vmin=1e-3,vmax=1)
        ax=plt.subplot2grid((3,2),(0,0),colspan=1,rowspan=1)
        img=ax.imshow(P_me[i,j],cmap=cmap,norm=norm)
        for k in range(len(pos)):
            ax.plot([pos[k],pos[k]],[0,nat],color='c')
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

        divider=make_axes_locatable(ax)
        cax=divider.append_axes("right",size="5%",pad=0.1)
        cbar=plt.colorbar(img,cax=cax)
        # cbar.ax.yaxis.set_major_locator(LogLocator(subs='all'))  ### For color parameters with small range, without spanning multiple orders of magnitude.
        # cbar.ax.yaxis.set_major_formatter(LogFormatter(minor_thresholds=(2,1)))  ### For color parameters with small range, without spanning multiple orders of magnitude.
        cbar.ax.tick_params(which='major',direction='in',width=wth,length=lenmin,labelsize=size)
        cbar.ax.tick_params(which='minor',direction='in',width=wth,length=0,labelsize=size)
        cbar.set_label('$P$',fontsize=size)
        cbar.outline.set_linewidth(wth)

        ax=plt.subplot2grid((3,2),(0,1),colspan=1,rowspan=1)
        ax.plot(P_nor1_he[i][j].reshape(-1),P_nor1_me[i][j].reshape(-1),'.',color='c',linewidth=wth)
        ax.plot([np.min(P_nor1_he[i][j]),np.max(P_nor1_he[i][j])],[np.min(P_nor1_he[i][j]),np.max(P_nor1_he[i][j])],color='k',linewidth=wth)
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax.xaxis.set_major_locator(MultipleLocator(xtick))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel('He',fontsize=size)
        ax.set_ylabel('Me',fontsize=size)

        x=np.arange(0,60,0.1)
        ax=plt.subplot2grid((3,2),(1,0),colspan=1,rowspan=1)
        lines=[]
        lines+=ax.plot(x,DP_he[i][j],linewidth=wth,label='He')
        lines+=ax.plot(x,DP_me[i][j],linewidth=wth,label='Me')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.xaxis.set_major_locator(MultipleLocator(xtick))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel('$\\xi$',fontsize=size)
        ax.set_ylabel('$P_\\xi$',fontsize=size)
        legend=ax.legend(loc='best',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
        for text,line in zip(legend.get_texts(),lines):
            text.set_color(line.get_color())

        x=np.arange(0,60,0.1)
        ax=plt.subplot2grid((3,2),(1,1),colspan=1,rowspan=1)
        lines=[]
        lines+=ax.plot(x[5:][::10],CS_he[i][j]*60,linewidth=wth,label='He')
        lines+=ax.plot(x,CS_me[i][j],linewidth=wth,label='Me')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax.xaxis.set_major_locator(MultipleLocator(xtick))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel('$\\xi$',fontsize=size)
        ax.set_ylabel('$P_\\xi$',fontsize=size)
        legend=ax.legend(loc='best',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
        for text,line in zip(legend.get_texts(),lines):
            text.set_color(line.get_color())

        x=np.arange(0,60,0.1)
        ax=plt.subplot2grid((3,2),(2,0),colspan=2,rowspan=1)
        lines=[]
        lines+=ax.plot(x,IS_he[i][j],linewidth=wth,label='He')
        lines+=ax.plot(x,IS_me[i][j],linewidth=wth,label='Me')
        ax.minorticks_on()
        ax.tick_params(axis='both',which='major',direction='in',width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis='both',which='minor',direction='in',width=wth,length=lenmin,labelsize=size)
        # ax.xaxis.set_major_locator(MultipleLocator(xtick))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        # ax.yaxis.set_major_locator(MultipleLocator(ytick[1][i]))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.spines['bottom'].set_linewidth(wth)
        ax.spines['top'].set_linewidth(wth)
        ax.spines['left'].set_linewidth(wth)
        ax.spines['right'].set_linewidth(wth)
        ax.set_xlabel('$\\xi$',fontsize=size)
        ax.set_ylabel('$P_\\xi$',fontsize=size)
        legend=ax.legend(loc='best',fontsize=0.8*size,handlelength=0.0,handletextpad=0.0)
        for text,line in zip(legend.get_texts(),lines):
            text.set_color(line.get_color())

        plt.tight_layout()
        plt.savefig(f'{figname}_{cl_[j][i]}_noice.png',format='png')
        plt.savefig(f'{figname}_{cl_[j][i]}_noice.pdf',format='pdf')

