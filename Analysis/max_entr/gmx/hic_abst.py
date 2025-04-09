# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib import colors
# from matplotlib.colors import LinearSegmentedColormap
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# import sys
# sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
# from freq_pbty import freq_pbty
#
# ifile='/hpc2hdd/home/chu-amat/cbfengphy/hic_fert/sperm_pwk_rep123_100000_iced.matrix'
# sysname=ifile.split('/')[len(ifile.split('/'))-1].split('_')[0]
# def f(idx):
#     chr_idx=idx
#     bin_size=100000
#     ofile=f'/hpc2hdd/home/chu-amat/cbfengphy/hic_fert/sperm_fig/{chr_idx}th_sperm_pwk_rep123_100000_iced.matrix'
#     lfile=ofile.split('/')
#     dir=''
#     for i in range(len(lfile)-1):
#         dir+='/'+lfile[i]
#     chr_size=[197195432,
#               181748087,
#               159599783,
#               155630120,
#               152537259,
#               149517037,
#               152524553,
#               131738871,
#               124076172,
#               129993255,
#               121843856,
#               121257530,
#               120284312,
#               125194864,
#               103494974,
#               98319150,
#               95272651,
#               90772031,
#               61342430]
#
#     f=open(ifile)
#     lsf=f.readlines()
#     f.close()
#     ls_f=[]
#     for lf in lsf:
#         l_f=lf.strip('\n').split()
#         ls_f.append(l_f)
#     ls_f=[list(map(float,sublist)) for sublist in ls_f]
#     N=int(ls_f[len(ls_f)-1][0])
#     P=np.zeros((N,N),dtype=float)
#     for i in range(len(ls_f)):
#         #if len(ls_f[i])>0: #In order to make sure the list index in next line does not exceed the range.
#         P[int(ls_f[i][0])-1,int(ls_f[i][1])-1]=ls_f[i][2]
#         P[int(ls_f[i][1])-1,int(ls_f[i][0])-1]=ls_f[i][2]
#
#     start=0
#     for i in range(chr_idx-1):
#         start+=chr_size[i]
#     end=start+chr_size[chr_idx-1]
#     start=int(start/bin_size)+6
#     end=int(end/bin_size)+20
#     P=P[start:end,start:end]
#     P=freq_pbty(P)
#
#     f=open(ofile,'w')
#     N=len(P)
#     for i in range(N):
#         for j in range(N):
#             f.write(f'{i+1} {j+1} {P[i,j]}\n')
#     f.close()
#
#     P[P<0.0001]=0.0001
#     ax=plt.subplot(111)
#     color=[(1,1,1),(1,0,0)]
#     nodes=[0.0,1.0]
#     cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
#     norm=colors.LogNorm(vmin=P.min(),vmax=P.max())
#     img=plt.imshow(P,cmap=cmap,norm=norm)
#     plt.title(f'{sysname} {chr_idx}-th chromosome')
#     plt.gca().invert_yaxis()
#     divider=make_axes_locatable(ax)
#     cax=divider.append_axes("right",size="5%",pad=0.05)
#     cbar=plt.colorbar(img,cax=cax)
#     plt.tight_layout()
#     plt.savefig(f'{dir}/{sysname}_{chr_idx}.png',format='png')
#     plt.savefig(f'{dir}/{sysname}_{chr_idx}.pdf',format='pdf')
#
# for i in range(19):
#     f(i+1)


import matplotlib.pyplot as plt
import numpy as np

# 创建一个随机矩阵
matrix = np.random.rand(500, 500)

# 显示矩阵
plt.imshow(matrix, cmap='viridis')

# 标记特定的矩阵元素
row = 20
col = 30
plt.annotate('.', xy=(col, row), color='red', fontsize=30, ha='center', va='center')

plt.show()