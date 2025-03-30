"""Import Modules"""
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.pyplot import MultipleLocator
from matplotlib import colors
from matplotlib.ticker import LogLocator,LogFormatter,AutoMinorLocator

file='test_zp_cor.npz'
figname=file[:-4]
data=np.load(file)
CS=data['CS']
P=data['P0_nor']
m=max(np.max(P),-np.min(P))
P[P==np.max(P)]=m
P[P==np.min(P)]=-m

CSs=np.zeros(3)
CSs[0]=np.mean(P[:10,:10])
CSs[1]=np.mean(P[-10:,-10:])
CSs[2]=np.mean(P[-10:,:10])
print(CSs)

# CS=np.array([1,0,0,1,0,1,1,1,0,1])
# P=np.array([[1,0,0,1,0,1,1,1,0,1],
#             [0,1,1,0,1,0,0,0,1,0],
#             [0,1,1,0,1,0,0,0,1,0],
#             [1,0,0,1,0,1,1,1,0,1],
#             [0,1,1,0,1,0,0,0,1,0],
#             [1,0,0,1,0,1,1,1,0,1],
#             [1,0,0,1,0,1,1,1,0,1],
#             [1,0,0,1,0,1,1,1,0,1],
#             [0,1,1,0,1,0,0,0,1,0],
#             [1,0,0,1,0,1,1,1,0,1]])
# CS[CS==0]=-1
# P[P==0]=-1
#
# sort_idx=np.argsort(CS)
# P=P[sort_idx][:,sort_idx]

plt.figure(figsize=(6,6))
wth=2
size=20
lenmaj=12
lenmin=8
lenbar=8
tick=200

ax=plt.subplot(111)
color=[(0,0,1),(1,1,1),(1,0,0)]
nodes=[0.0,0.5,1.0]
cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
norm=colors.Normalize(vmin=-1,vmax=1)
ax.plot(CS-np.min(CS),color=(0,1,1),linewidth=wth)
img=plt.imshow(P,cmap=cmap,norm=norm)
plt.gca().invert_yaxis()
plt.minorticks_on()
plt.tick_params(axis="both",which="major",direction="out",width=wth,length=lenmaj,labelsize=size)
plt.tick_params(axis="both",which="minor",direction="out",width=wth,length=lenmin,labelsize=size)
for spine in ax.spines.values():
        spine.set_linewidth(wth)
divider=make_axes_locatable(ax)
cax=divider.append_axes("right",size="5%",pad=0.05)
cbar=plt.colorbar(img,cax=cax)
cbar.ax.tick_params(direction="out",width=wth,length=lenbar,labelsize=size,labelcolor='black')
cbar.outline.set_linewidth(wth)

plt.savefig(f'{figname}.png')
plt.show()
