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
from dist_prob import dist_prob

reffile=[]
reffile+=['e2cp.dat']
reffile+=['e2cm.dat']
reffile+=['l2cp.dat']
reffile+=['l2cm.dat']
reffile+=['8cp.dat']
reffile+=['8cm.dat']
figname='plot'
nat=600

P_ref=np.zeros((len(reffile),nat,nat))
for i in range(len(reffile)):
    f=open(reffile[i],'r')
    lsf=f.readlines()
    for j in range(len(lsf)):
        ls_f=lsf[j].strip('\n').split()
        P_ref[i,int(ls_f[0])-1,int(ls_f[1])-1]=float(ls_f[2])
# for i in range(len(reffile)):
#     dp=dist_prob(P_ref[i])
#     norm=np.zeros((nat,nat))
#     for j in range(nat):
#         for k in range(nat):
#             norm[j,k]=dp[abs(j-k)]
#     P_ref[i]-=norm

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
binsize=100000
x_positions=np.arange(0.5,nat+0.5+1,tick)
y_positions=np.arange(0.5,nat+0.5+1,tick)
x_labels=np.arange(0,int(binsize*nat/scale)+1,int(binsize*tick/scale))
y_labels=np.arange(0,int(binsize*nat/scale)+1,int(binsize*tick/scale))
subp=[231,232,233,234,235,236]

plt.figure(figsize=(14,10))
for i in range(len(reffile)):
    ax=plt.subplot(subp[i])
    color=[(1,1,1),(1,0,0)]
    nodes=[0.0,1.0]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=colors.LogNorm()
    img=plt.imshow(P_ref[i],cmap=cmap,norm=norm)
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

plt.tight_layout()
plt.savefig(f'{figname}.png',format='png')
plt.savefig(f'{figname}.pdf',format='pdf')
plt.show()