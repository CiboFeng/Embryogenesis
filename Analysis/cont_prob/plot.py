"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
from split_data import split_data
from color_bar import color_bar

"""Set Arguments"""
file='seq.pyw'
cl=['ZP','ZM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
cl_=['zp','zm','esc']
sw_=['zp_esc','zm_esc']
nfr=56
nat=600
binsize=1000000
t=[0,0.01,0.1,1.0,10.0,100.0,1000.0,10000.0]
t=[0,0.01]

"""Define Funtions for Plot"""
def plot2(P,P_nor1,P_nor2,figname):
        wth=2
        size=20
        lenmaj=12
        lenmin=8
        lenbar=8
        tick=200
        ### Diaplay labels every #tick matrix indexes.
        unit='mb'
        if unit=='kb':
            scale=1000
            ### Unit conversion scale.
        if unit=='mb':
            scale=1000000
        x_positions=np.arange(0.5,nat+0.5+1,tick)
        y_positions=np.arange(0.5,nat+0.5+1,tick)
        x_labels=np.arange(0,int(binsize*nat/scale)+1,int(binsize*tick/scale))
        y_labels=np.arange(0,int(binsize*nat/scale)+1,int(binsize*tick/scale))

        plt.figure(figsize=(14,4))

        ax=plt.subplot(131)
        color=[(1,1,1),(1,0,0)]
        nodes=[0.0,1.0]
        cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
        norm=colors.LogNorm()
        img=plt.imshow(P,cmap=cmap,norm=norm)
        plt.xticks(x_positions,x_labels)
        plt.yticks(y_positions,y_labels)
        plt.gca().invert_yaxis()
        # img.tick_params(axis='x',labeltop=True,labelbottom=False)
        plt.xlabel('Base pairs in '+unit,fontsize=size)
        plt.ylabel('Base pairs in '+unit,fontsize=size)
        plt.title('Contact Probability',fontsize=size)
        plt.minorticks_on()
        plt.tick_params(axis="both",which="major",direction="out",width=wth,length=lenmaj,labelsize=size)
        plt.tick_params(axis="both",which="minor",direction="out",width=wth,length=lenmin,labelsize=size)
        for spine in ax.spines.values():
                spine.set_linewidth(wth)
        divider=make_axes_locatable(ax)
        cax=divider.append_axes("right",size="5%",pad=0.05)
        cbar=plt.colorbar(img,cax=cax)
        cbar.ax1.tick_params(direction="out",width=wth,length=lenbar,labelsize=size,labelcolor='black')
        cbar.outline.set_linewidth(wth)

        ax=plt.subplot(132)
        color=[(0,0,1),(1,1,1),(1,0,0)]
        nodes=[0.0,0.5,1.0]
        cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
        norm=colors.Normalize()
        img=plt.imshow(P_nor1,cmap=cmap,norm=norm)
        plt.xticks(x_positions,x_labels)
        plt.yticks(y_positions,[])
        plt.gca().invert_yaxis()
        plt.xlabel('Base pairs in '+unit,fontsize=size)
        plt.title('$\mathrm{log}_2(obs/exp)$',fontsize=size)
        plt.minorticks_on()
        plt.tick_params(axis="both",which="major",direction="out",width=wth,length=lenmaj,labelsize=size)
        plt.tick_params(axis="both",which="minor",direction="out",width=wth,length=lenmin,labelsize=size)
        for spine in ax.spines.values():
                spine.set_linewidth(wth)
        divider=make_axes_locatable(ax)
        cax=divider.append_axes("right",size="5%",pad=0.05)
        cbar=plt.colorbar(img,cax=cax)
        cbar.ax1.tick_params(direction="out",width=wth,length=lenbar,labelsize=size,labelcolor='black')
        cbar.outline.set_linewidth(wth)

        ax=plt.subplot(133)
        color,nodes=color_bar(np.min(P_nor2),np.max(P_nor2))
        cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
        norm=colors.Normalize()
        img=plt.imshow(P_nor2,cmap=cmap,norm=norm)
        plt.xticks(x_positions,x_labels)
        plt.yticks(y_positions,[])
        plt.gca().invert_yaxis()
        plt.xlabel('Base pairs in '+unit,fontsize=size)
        # plt.title('Compared with average',fontsize=size)
        plt.title('Correlation coefficient',fontsize=size)
        plt.minorticks_on()
        plt.tick_params(axis="both",which="major",direction="out",width=wth,length=lenmaj,labelsize=size)
        plt.tick_params(axis="both",which="minor",direction="out",width=wth,length=lenmin,labelsize=size)
        for spine in ax.spines.values():
                spine.set_linewidth(wth)
        divider=make_axes_locatable(ax)
        cax=divider.append_axes("right",size="5%",pad=0.05)
        cbar=plt.colorbar(img,cax=cax)
        cbar.ax1.tick_params(direction="out",width=wth,length=lenbar,labelsize=size,labelcolor='black')
        cbar.outline.set_linewidth(wth)

        plt.tight_layout()
        plt.savefig(f'{figname}.png',format='png')
        plt.savefig(f'{figname}.pdf',format='pdf')
        # plt.show()

def plot1(x,CD,x_,CP_,CP,IS,figname):
        wth=3
        size=30
        lenmaj=15
        lenmin=8
        xtick=200

        plt.figure(figsize=(10,15))

        ax=plt.subplot(311)
        plt.plot(x,CD,'-',color='c',linewidth=wth)
        plt.xscale('log')
        plt.yscale('log')
        ax.minorticks_on()
        ax.tick_params(axis="both",which="major",direction="in",width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis="both",which="minor",direction="in",width=wth,length=lenmin,labelsize=size)
        # plt.xlim(115,225)
        # plt.ylim(2,32)
        axis=plt.gca()
        # axis.xaxis.set_major_locator(MultipleLocator(xtick))
        # axis.xaxis.set_minor_locator(MultipleLocator(xtick/2))
        # axis.yaxis.set_major_locator(MultipleLocator(0.5))
        # axis.yaxis.set_minor_locator(MultipleLocator(0.25))
        axis.spines['bottom'].set_linewidth(wth)
        axis.spines['top'].set_linewidth(wth)
        axis.spines['left'].set_linewidth(wth)
        axis.spines['right'].set_linewidth(wth)
        plt.xlabel('Sequence distance',fontsize=size)
        plt.ylabel('Contact probability',fontsize=size)

        ax=plt.subplot(312)
        for i in range(len(CP_)):
            if np.mean(CP_[i])>0:
                plt.plot(np.array(x_[i])+1,CP_[i],'-',color='r',linewidth=wth)
            else:
                plt.plot(np.array(x_[i])+1,CP_[i],'-',color='b',linewidth=wth)
        ax.minorticks_on()
        ax.tick_params(axis="both",which="major",direction="in",width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis="both",which="minor",direction="in",width=wth,length=lenmin,labelsize=size)
        # plt.xlim(115,225)
        plt.ylim(np.min(CP)-((np.max(CP)-np.min(CP))/10),np.max(CP)+((np.max(CP)-np.min(CP))/10))
        axis=plt.gca()
        axis.xaxis.set_major_locator(MultipleLocator(xtick))
        axis.xaxis.set_minor_locator(MultipleLocator(xtick/5))
        # axis.yaxis.set_major_locator(MultipleLocator(10))
        # axis.yaxis.set_minor_locator(MultipleLocator(5))
        axis.spines['bottom'].set_linewidth(wth)
        axis.spines['top'].set_linewidth(wth)
        axis.spines['left'].set_linewidth(wth)
        axis.spines['right'].set_linewidth(wth)
        plt.xlabel('Index of particles',fontsize=size)
        plt.ylabel('Compartment profile',fontsize=size)

        ax=plt.subplot(313)
        plt.plot(x+1,IS,'-',color='c',linewidth=wth)
        ax.minorticks_on()
        ax.tick_params(axis="both",which="major",direction="in",width=wth,length=lenmaj,labelsize=size)
        ax.tick_params(axis="both",which="minor",direction="in",width=wth,length=lenmin,labelsize=size)
        # plt.xlim(115,225)
        # plt.ylim(2,32)
        axis=plt.gca()
        axis.xaxis.set_major_locator(MultipleLocator(xtick))
        axis.xaxis.set_minor_locator(MultipleLocator(xtick/5))
        # axis.yaxis.set_major_locator(MultipleLocator(10))
        # axis.yaxis.set_minor_locator(MultipleLocator(5))
        axis.spines['bottom'].set_linewidth(wth)
        axis.spines['top'].set_linewidth(wth)
        axis.spines['left'].set_linewidth(wth)
        axis.spines['right'].set_linewidth(wth)
        plt.xlabel('Index of particles',fontsize=size)
        plt.ylabel('Insulation score',fontsize=size)

        plt.tight_layout()
        plt.savefig(f'{figname}.png', format='png')
        plt.savefig(f'{figname}.pdf', format='pdf')
        # plt.show()

"""Read Data, Make Some Modifications in Order to Make Plot More Convenient, and Plot"""
for i in range(ncl):
        P0=pyw(file,f'The Contact Probability of {cl[i]}:')[0]
        P_nor10=pyw(file,f'The Obs/Exp Matrix of {cl[i]}:')[0]
        P_nor20=pyw(file,f'The Correlation Coefficient Matrix of {cl[i]}:')[0]
        DP=pyw(file,f'The Contact Probability v.s. Sequence Distance of {cl[i]}:')[0]
        CS=pyw(file,f'The Compartment Signal of {cl[i]}:')[0]
        IS=pyw(file,f'The Insulation Score of {cl[i]}:')[0]

        P=np.zeros((nat,nat))
        for j in range(nat):
                for k in range(nat):
                        P[j,k]=P0[j*nat+k]
        P=np.where(P==np.max(P),1,P)
        P_nor1=np.zeros((nat,nat))
        for j in range(nat):
                for k in range(nat):
                        P_nor1[j,k]=P_nor10[j*nat+k]
        mm=max([np.abs(np.max(P_nor1)),np.abs(np.min(P_nor1))])
        P_nor1=np.where(P_nor1==np.max(P_nor1),mm,P_nor1)
        P_nor1=np.where(P_nor1==np.min(P_nor1),-mm,P_nor1)
        P_nor2=np.zeros((nat,nat))
        for j in range(nat):
                for k in range(nat):
                        P_nor2[j,k]=P_nor20[j*nat+k]
        mm=max([np.abs(np.max(P_nor2)),np.abs(np.min(P_nor2))])
        P_nor2=np.where(P_nor2==np.max(P_nor2),mm,P_nor2)
        P_nor2=np.where(P_nor2==np.min(P_nor2),-mm,P_nor2)
        x=np.linspace(0,nat-1,nat)
        x_,CP_=split_data(x,CS,0)

        plot2(P,P_nor1,P_nor2,f'mat_ref_{cl_[i]}')
        plot1(x,DP,x_,CP_,CS,IS,f'seq_ref_{cl_[i]}')

for i in range(nsw):
        for j in range(len(t)):
                P0=pyw(file,f'The Contact Probability of {sw[i]} at t={t[j]}:')[0]
                P_nor10=pyw(file,f'The Obs/Exp Matrix of {sw[i]} at t={t[j]}:')[0]
                P_nor20=pyw(file,f'The Correlation Coefficient Matrix of {sw[i]} at t={t[j]}:')[0]
                DP=pyw(file,f'The Contact Probability v.s. Sequence Distance of {sw[i]} at t={t[j]}:')[0]
                CS=pyw(file,f'The Compartment Signal of {sw[i]} at t={t[j]}:')[0]
                IS=pyw(file,f'The Insulation Score of {sw[i]} at t={t[j]}:')[0]

                P=np.zeros((nat,nat))
                for k in range(nat):
                        for l in range(nat):
                                P[k,l]=P0[k*nat+l]
                P=np.where(P==np.max(P),1,P)
                P_nor1=np.zeros((nat,nat))
                for k in range(nat):
                        for l in range(nat):
                                P_nor1[k,l]=P_nor10[k*nat+l]
                mm=max([np.abs(np.max(P_nor1)),np.abs(np.min(P_nor1))])
                P_nor1=np.where(P_nor1==np.max(P_nor1),mm,P_nor1)
                P_nor1=np.where(P_nor1==np.min(P_nor1),-mm,P_nor1)
                P_nor2=np.zeros((nat,nat))
                for k in range(nat):
                        for l in range(nat):
                                P_nor2[k,l]=P_nor20[k*nat+l]
                mm=max([np.abs(np.max(P_nor2)),np.abs(np.min(P_nor2))])
                P_nor2=np.where(P_nor2==np.max(P_nor2),mm,P_nor2)
                P_nor2=np.where(P_nor2==np.min(P_nor2),-mm,P_nor2)
                x=np.linspace(0,nat-1,nat)
                x_,CP_=split_data(x,CS,0)

                plot2(P,P_nor1,P_nor2,f'mat_{sw_[i]}_{t[j]}')
                plot1(x,DP,x_,CP_,CS,IS,f'seq_{sw_[i]}_{t[j]}')



