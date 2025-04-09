"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter
from matplotlib.pyplot import MultipleLocator
import matplotlib.ticker as ticker
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw

"""Set Arguments"""
file='r2_cont_prob_cross.npy'
figname=file[:-4]
cl=['ZP','ZM','ESC']
sw=['ZP-ESC','ZM-ESC']
ncl=len(cl)
nsw=len(sw)
nfr=56

"""Read Data"""
r2=np.load(file)

"""Plot"""
plt.figure(figsize=(6,5))
wth=2
size=20
lenmaj=12
lenmin=8
lenbar=8
tick=9
subp=[111]
x_positions=np.arange(1.5,nfr-9+0.5+1,tick)
y_positions=np.arange(1.5,nfr-9+0.5+1,tick)
# x_labels=np.array([-2,-1,0,1,2,3,4])
x_labels=np.array(['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])
# y_labels=np.array([-2,-1,0,1,2,3,4])
y_labels=np.array(['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])
color=[(1,1,1),(1,1,0),(1,0,0),(0,0,0)]
nodes=[0/3,1/3,2/3,3/3]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize()

for i in [0]:
        for j in [1]:
                ax=plt.subplot(subp[i])
                img=plt.imshow(r2[i,j,:-9,:-9],cmap=cmap,norm=norm)
                plt.xticks(x_positions,x_labels)
                plt.yticks(y_positions,y_labels)
                plt.gca().invert_yaxis()
                # img.tick_params(axis='x',labeltop=True,labelbottom=False)
                plt.xlabel(f'Time in {sw[i]}',fontsize=size)
                plt.ylabel(f'Time in {sw[j]}',fontsize=size)
                plt.minorticks_on()
                plt.tick_params(axis='both',which='major',direction='out',width=wth,length=lenmaj,labelsize=size)
                plt.tick_params(axis='both',which='minor',direction='out',width=wth,length=0,labelsize=size)
                for spine in ax.spines.values():
                        spine.set_linewidth(wth)
                divider=make_axes_locatable(ax)
                cax=divider.append_axes('right',size='5%',pad=0.05)
                cbar=plt.colorbar(img,cax=cax)
                cbar.set_label('$R^2$',fontsize=size)
                cbar.ax.tick_params(direction='in',width=wth,length=lenbar,labelsize=size,labelcolor='black')
                cbar.outline.set_linewidth(wth)
                plt.tight_layout()
plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()



