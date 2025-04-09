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
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
from pyw import pyw

"""Set Arguments"""
file='shp_seg.npy'
figname=file[:-4]
nat=600

"""Read Data"""
q=np.load(file)
for i in range(nat):
    for j in range(i):
        q[:,:,j,i,:,:]=q[:,:,i,j,:,:]

"""Plot"""
wth=2
size=20
lenmaj=12
lenmin=8
lenbar=8
tick=200
### Diaplay labels every #tick matrix indexes.
unit='Mb'
if unit=='kb':
    scale=1000
    ### Unit conversion scale.
if unit=='Mb':
    scale=1000000
binsize=100000
x_positions=np.arange(0.5,nat+0.5+1,tick)
y_positions=np.arange(0.5,nat+0.5+1,tick)
x_labels=np.arange(0,int(binsize*nat/scale)+1,int(binsize*tick/scale))
y_labels=np.arange(0,int(binsize*nat/scale)+1,int(binsize*tick/scale))
subp=[231,232,233,234,235,236]
title=['R_{\\mathrm{PA1}}~(\\sigma)',
       'R_{\\mathrm{PA2}}~(\\sigma)',
       'R_{\\mathrm{PA3}}~(\\sigma)',
       'R_{\\mathrm{g}}~(\\sigma)',
       '\\Delta',
       'S']

plt.figure(figsize=(14,8))
for i in range(np.shape(q)[-1]):
    ax=plt.subplot(subp[i])
    color=[(1,1,1),(1,0,0)]
    nodes=[0.0,1.0]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=colors.Normalize()
    img=plt.imshow(q[0,0,:,:,0,i],cmap=cmap,norm=norm)
    plt.gca().invert_yaxis()
    # img.tick_params(axis='x',labeltop=True,labelbottom=False)
    if subp[i] in [234,235,236]:
        plt.xticks(x_positions,x_labels)
        plt.xlabel('Base pairs in '+unit,fontsize=size)
    else:
        plt.xticks([],[])
    if subp[i] in [231,234]:
        plt.ylabel('Base pairs in '+unit,fontsize=size)
        plt.yticks(y_positions,y_labels)
    else:
        plt.yticks([],[])
    plt.minorticks_on()
    plt.tick_params(axis='both',which='major',direction="out",width=wth,length=lenmaj,labelsize=size)
    plt.tick_params(axis='both',which='minor',direction="out",width=wth,length=lenmin,labelsize=size)
    for spine in ax.spines.values():
            spine.set_linewidth(wth)
    divider=make_axes_locatable(ax)
    cax=divider.append_axes("right",size="5%",pad=0.05)
    cbar=plt.colorbar(img,cax=cax)
    cbar.set_label(f'${title[i]}$',fontsize=size)
    cbar.ax1.tick_params(direction='in',width=wth,length=lenbar,labelsize=size,labelcolor='black')
    cbar.outline.set_linewidth(wth)

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()