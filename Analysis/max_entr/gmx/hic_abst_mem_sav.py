import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from freq_pbty import freq_pbty

ifile='/hpc2hdd/home/cfeng593/proj/hicpro/rbl/output/hic_results/matrix/SRR22859374/iced/100000/SRR22859374_100000_iced.matrix'
sysname=ifile.split('/')[len(ifile.split('/'))-1].split('_')[0]
def f(idx):
    chr_idx=idx
    bin_size=100000
    ofile=f'/hpc2hdd/home/cfeng593/proj/hicpro/rbl/output/hic_results/matrix/fig/rbl_100000_iced.matrix'
    lfile=ofile.split('/')
    dir=''
    for i in range(len(lfile)-1):
        dir+='/'+lfile[i]
    chr_size=[249250621,
            243199373,
            198022430,
            191154276,
            180915260,
            171115067,
            159138663,
            146364022,
            141213431,
            135534747,
            135006516,
            133851895,
            115169878,
            107349540,
            102531392,
            90354753,
            81195210,
            78077248,
            59128983,
            63025520,
            48129895,
            51304566,
            155270560,
            59373566,
            16571]

    start=0
    for i in range(chr_idx-1):
        start+=chr_size[i]
    start+=116000000
    end=start+2000000
    start=int(start/bin_size)+6-6
    end=int(end/bin_size)+20-20
    N=end-start
    P=np.zeros((N,N),dtype=float)

    f=open(ifile)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for lf in lsf:
        l_f=lf.strip('\n').split()
        ls_f.append(l_f)
    ls_f=[list(map(float,sublist)) for sublist in ls_f]
    for i in range(len(ls_f)):
        if int(ls_f[i][0])-1>=start and int(ls_f[i][0])-1<end and int(ls_f[i][1])-1>=start and int(ls_f[i][1])-1<end:
            P[int(ls_f[i][0])-1-start,int(ls_f[i][1])-1-start]=ls_f[i][2]
            P[int(ls_f[i][1])-1-start,int(ls_f[i][0])-1-start]=ls_f[i][2]

    P=freq_pbty(P)

    f=open(ofile,'w')
    N=len(P)
    for i in range(N):
        for j in range(N):
            f.write(f'{i+1} {j+1} {P[i,j]}\n')
    f.close()

    ax=plt.subplot(111)
    color=[(1,1,1),(1,0,0)]
    nodes=[0.0,1.0]
    cmap=LinearSegmentedColormap.from_list('custom_cmap',list(zip(nodes,color)))
    norm=colors.LogNorm()
    img=plt.imshow(P,cmap=cmap,norm=norm)
    plt.title(f'{sysname} {chr_idx}-th chromosome')
    plt.gca().invert_yaxis()
    divider=make_axes_locatable(ax)
    cax=divider.append_axes("right",size="5%",pad=0.05)
    cbar=plt.colorbar(img,cax=cax)
    plt.tight_layout()
    plt.savefig(f'{dir}/{sysname}.png',format='png')
    plt.savefig(f'{dir}/{sysname}.pdf',format='pdf')

f(1)
