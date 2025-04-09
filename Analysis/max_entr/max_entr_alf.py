"""Import Modules"""
import argparse
import os
import numpy as np
import datetime
import scipy.io
import subprocess
import h5py
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from matplotlib import colors
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw
from dist_prob import dist_prob
from cpmt_sgnl import cpmt_sgnl
from insul_score import insul_score
from tot_len import tot_len

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-s",type=str,help="The sysname name")
parser.add_argument("-w",type=str,help="The work directory")
parser.add_argument("-d",type=str,help="The directory of .npz files")
parser.add_argument("-i",type=str,help="The inputting alpha file")
parser.add_argument("-o",type=str,help="The outputting alpha file")
args=parser.parse_args()
sysname=args.s
wkdir=args.w
dir=args.d
alfifile=args.i
alfofile=args.o
mpirfile=wkdir+f'/mtn_pair_{sysname}.pyw'
tpirfile=wkdir+f'/trn_pair_{sysname}.pyw'
pexpfile=wkdir+f'/exp_cont_prob_{sysname}.pyw'
chkfile=wkdir+f'/chk_{sysname}.pyw'
ipathls=alfifile.split('/')[:len(alfifile.split('/'))-1]
alyspath=''
for i in range(len(ipathls)):
    alyspath+=ipathls[i]+'/'
alysfile=alyspath+f'alys_{sysname}.pyw'
kB=8.314
T=1/kB
ord=int(dir.split('_')[len(dir.split('_'))-1])
ord_1=20
scl_1=0.01

"""Read Maintained Pairs, Trained Pairs, Experimental Contact Probabilities, and the Alfs of the Last Step, from the File"""
if os.path.exists(mpirfile):
    f=open(mpirfile)
    lsf=f.readlines()
    f.close()
    ls_f=[]
    for i in range(1,len(lsf)):
        l_f=lsf[i].strip('\n').split()
        ls_f.append(l_f)
    mpir=[list(map(int,subls)) for subls in ls_f]
else:
    mpir=[]
nmpir=len(mpir)

f=open(tpirfile)
lsf=f.readlines()
f.close()
ls_f=[]
for i in range(1,len(lsf)):
    l_f=lsf[i].strip('\n').split()
    ls_f.append(l_f)
tpir=[list(map(int,subls)) for subls in ls_f]
ntpir=len(tpir)

p_exp_all=pyw(pexpfile,'Probability:')[0]
alf_i_all=pyw(alfifile,'Alpha')[0]
N=int((len(p_exp_all))**0.5)
P_exp=np.zeros((N,N))
Alf_i=np.zeros((N,N))
for i in range(N):
    for j in range(N):
        P_exp[i,j]=p_exp_all[i*N+j]
        Alf_i[i,j]=alf_i_all[i*N+j]
for i in range(N-1):
    P_exp[i,i+1]=1
    P_exp[i+1,i]=1

p_exp=[]
alf_i=[]
alf_0=[]
for i in range(ntpir):
    p_exp+=[P_exp[tpir[i][0],tpir[i][1]]]
    alf_i+=[Alf_i[tpir[i][0],tpir[i][1]]]
for i in range(nmpir):
    alf_0+=[Alf_i[mpir[i][0],mpir[i][1]]]
p_exp=np.array(p_exp)
alf_i=np.array(alf_i)
alf_0=np.array(alf_0)

"""Read f_sim_ave, C, and P_sim from .npz Files"""
files=[f for f in os.listdir(dir) if f.endswith('.npz')]
print('The Number of Trajectories:',len(files))
p_sim=np.zeros((1,ntpir))
C=np.zeros((ntpir,ntpir))
P_sim=np.zeros((N,N))
for i in range(len(files)):
    data=np.load(dir+f'/{files[i]}')
    p_sim+=data['p']
    C+=data['C']
    P_sim+=data['P']
    os.remove(dir+f'/{files[i]}')
p_sim/=len(files)
C/=len(files)
P_sim/=len(files)
for i in range(N-1):
    P_sim[i,i+1]=1
    P_sim[i+1,i]=1

"""Update Alphas"""
D=np.matmul(p_sim.T,p_sim)
B=C-D
now=datetime.datetime.now()
print(f'Complete B',now)
tf1=wkdir+'/temp_file1.mat'
tf2=wkdir+'/temp_file2.mat'
# scipy.io.savemat(tf1,{'B':B})
with h5py.File(tf1,'w') as f:
    f.create_dataset('B',data=B)
#matlab_command=f"/hpc2ssd/softwares/matlab/bin/matlab -r \"data=load(\'{tf1}\');B=data.B;B_=pinv(B);save(\'{tf2}\',\'B_\',\'-v7.3\');exit;\""
matlab_command=f"/hpc2ssd/softwares/matlab/bin/matlab -r \"B=h5read(\'{tf1}\',\'/B\');B_=pinv(B);save(\'{tf2}\',\'B_\',\'-v7.3\');exit;\""
subprocess.check_output(matlab_command,shell=True,text=True)
# mat=scipy.io.loadmat(tf2)
# B_=mat['B_']
with h5py.File(tf2,'r') as f:
    B_=f['B_'][:]
os.remove(tf1)
os.remove(tf2)
now=datetime.datetime.now()
print(f'Complete B^T',now)
p_dlt=(p_sim-p_exp).reshape(-1,1)
### For one dimentional array, neither "#.T" nor "np.transpose(#)" can turn it into a column vector. For
### this case, they all can do.
alf_dlt=kB*T*np.matmul(B_,p_dlt)
if ord<ord_1:
    scl=scl_1**((ord-ord_1)/(1-ord_1))
else:
    scl=1
print('Scaling Factor:',scl)
alf_dlt*=scl
alf_dlt=alf_dlt.reshape(-1)
alf_o=alf_i+alf_dlt

Alf_o=np.zeros((N,N))
for i in range(nmpir):
    Alf_o[mpir[i][0],mpir[i][1]]=alf_0[i]
    Alf_o[mpir[i][1],mpir[i][0]]=alf_0[i]
for i in range(ntpir):
    Alf_o[tpir[i][0],tpir[i][1]]=alf_o[i]
    Alf_o[tpir[i][1],tpir[i][0]]=alf_o[i]

"""Calculate Three Analysis Quantities: Sequence-distance-dependent Contact Probability, Compartment Signal, and Insulation Score"""
d=np.linspace(0,N-1,N)
dp_exp=dist_prob(P_exp)
dp_sim=dist_prob(P_sim)

cs_exp=cpmt_sgnl(P_exp)[0]
cs_sim=cpmt_sgnl(P_sim)[0]
cs_cor=np.corrcoef(cs_sim,cs_exp)[0,1]
if cs_cor<0:
    cs_sim*=-1

is_exp=insul_score(P_exp,bin_size=10**5)[0]
is_sim=insul_score(P_sim,bin_size=10**5)[0]

"""Compare Simulation Results with Experiment and Record"""
p_dlt_sum=0
p_exp_sum=0
for i in range(ntpir):
    p_dlt_sum+=abs(p_dlt[i,0])
    p_exp_sum+=p_exp[i]
rel_err1=p_dlt_sum/p_exp_sum

P_dlt=P_sim-P_exp
rel_err2=[]
for i in range(N):
    for j in range(N):
        if P_exp[i,j]!=0:
            rel_err2+=[abs(P_dlt[i,j])/P_exp[i,j]]
rel_err2=np.mean(rel_err2)

p_exp_all=P_exp.reshape(-1)
p_sim_all=P_sim.reshape(-1)
r2_p=r2_score(p_exp_all,p_sim_all)

r2_dp=r2_score(dp_exp,dp_sim)
r2_cs=r2_score(cs_exp,cs_sim)
r2_is=r2_score(is_exp,is_sim)

"""Output Data"""
if not os.path.exists(chkfile):
    chk=open(chkfile,'w')
    chk.write(f'The Number of Iteration Steps, the Relative Error of Contact Probabilities among the Trained Pairs, and '
              f'among the Whole Matrix, the Coefficient of Determination of Contact Probabilities, Dequence-distance-'
              f'dependent Contact Probability, Compartment Signal, and Insulation Score:\n')
    chk.close()
chk=open(chkfile,'a')
chk.write(f'{ord}\t{tot_len(rel_err1,8)}\t{tot_len(rel_err2,8)}\t{tot_len(r2_p,8)}\t{tot_len(r2_dp,8)}\t{tot_len(r2_cs,8)}\t{tot_len(r2_is,8)}\n')
chk.close()

pyw=open(alfofile,'w')
pyw.write('Alpha of Data-driven Term:\n')
for i in range(N):
    for j in range(N):
        pyw.write(f'{Alf_o[i,j]}\n')

pyw=open(alysfile,'w')
pyw.write('Contact Probability:\n')
for i in range(N):
    for j in range(N):
        pyw.write(f'{P_sim[i,j]}\n')
pyw.write('\n')
pyw.write('Sequence-distance-dependent Contact Probability, Compartment Signal, Insulation Score:\n')
for i in range(N):
    pyw.write(f'{dp_sim[i]}\t{cs_sim[i]}\t{is_sim[i]}\n')

"""Plot"""
Compare=np.zeros((N,N))
for i in range(N):
    for j in range(N):
        if i>j:
            Compare[i,j]=P_exp[i,j]
        else:
            Compare[i,j]=P_sim[i,j]
Compare[Compare==0]=np.min(Compare[Compare>0])/2

plt.figure()
norm=colors.SymLogNorm(0.001)
plt.imshow(Alf_i,norm=norm)
plt.gca().invert_yaxis()
plt.xlabel('Index of Particles')
plt.ylabel('Index of Particles')
plt.colorbar()
plt.savefig(alyspath+f'alf_{sysname}.png',format='png')
plt.savefig(alyspath+f'alf_{sysname}.pdf',format='pdf')

plt.figure()
norm=colors.LogNorm()
plt.imshow(Compare,norm=norm)
plt.gca().invert_yaxis()
plt.title('Upper Left: Exp, Lower Right: Sim')
plt.xlabel('Index of Particles')
plt.ylabel('Index of Particles')
plt.colorbar()
plt.savefig(alyspath+f'cont_prob_{sysname}.png',format='png')
plt.savefig(alyspath+f'cont_prob_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(p_exp_all,p_sim_all,'.')
plt.plot([min(p_exp_all),max(p_exp_all)],[min(p_exp_all),max(p_exp_all)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Contact Probabilities of Experiment')
plt.ylabel('Contact Probabilities of Simulation')
plt.savefig(alyspath+f'cont_prob_cor_{sysname}.png',format='png')
plt.savefig(alyspath+f'cont_prob_cor_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(d,dp_exp,label='Exp')
plt.plot(d,dp_sim,label='Sim')
plt.legend(loc='best')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Sequence Distance')
plt.ylabel('Sequence-distance-dependent Contact Probability')
plt.savefig(alyspath+f'dist_prob_{sysname}.png',format='png')
plt.savefig(alyspath+f'dist_prob_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(dp_exp,dp_sim,'.')
plt.plot([min(dp_exp),max(dp_exp)],[min(dp_exp),max(dp_exp)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Sequence-distance-dependent\nContact Probability of Experiment')
plt.ylabel('Sequence-distance-dependent\nContact Probability of Simulation')
plt.savefig(alyspath+f'dist_prob_cor_{sysname}.png',format='png')
plt.savefig(alyspath+f'dist_prob_cor_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(d+1,cs_exp,label='Exp')
plt.plot(d+1,cs_sim,label='Sim')
plt.legend(loc='best')
plt.xlabel('Index of Particles')
plt.ylabel('Compartment Signal')
plt.savefig(alyspath+f'cpmt_sgnl_{sysname}.png',format='png')
plt.savefig(alyspath+f'cpmt_sgnl_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(cs_exp,cs_sim,'.')
plt.plot([min(cs_exp),max(cs_exp)],[min(cs_exp),max(cs_exp)])
plt.xlabel('Compartment Signal of Experiment')
plt.ylabel('Compartment Signal of Simulation')
plt.savefig(alyspath+f'cpmt_sgnl_cor_{sysname}.png',format='png')
plt.savefig(alyspath+f'cpmt_sgnl_cor_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(d+1,is_exp,label='Exp')
plt.plot(d+1,is_sim,label='Sim')
plt.legend(loc='best')
plt.xlabel('Index of Particles')
plt.ylabel('Insulation Score')
plt.savefig(alyspath+f'insul_score_{sysname}.png',format='png')
plt.savefig(alyspath+f'insul_score_{sysname}.pdf',format='pdf')

plt.figure()
plt.plot(is_exp,is_sim,'.')
plt.plot([min(is_exp),max(is_exp)],[min(is_exp),max(is_exp)])
plt.xlabel('Insulation Score of Experiment')
plt.ylabel('Insulation Score of Simulation')
plt.savefig(alyspath+f'insul_score_cor_{sysname}.png',format='png')
plt.savefig(alyspath+f'insul_score_cor_{sysname}.pdf',format='pdf')
