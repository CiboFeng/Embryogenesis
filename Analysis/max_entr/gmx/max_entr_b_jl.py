# """Import Modules"""
# import datetime
# import numpy as np
# import multiprocessing
# import subprocess
# import time
# import tempfile
# import os
# import sys
# sys.path.append('D:\Biophysics\MyPython\\functions')
# sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
# from file_par import file_par
# from pdbs import pdbs
# from cont import cont
# from pyw import pyw
# from cont_dist import cont_dist
# from comprtmt import comprtmt
# from insul_score import insul_score
# import scipy.io
# import matplotlib
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib import colors
# from matplotlib.colors import LinearSegmentedColormap
# from mpl_toolkits.axes_grid1 import make_axes_locatable
#
# def PC(i,q1,q2,dir,file,dt,r0,hipair,wkdir):
#     """Calculate the contact probability matrix as a flat vector and the cross-averaged matrix."""
#
#     tf1=wkdir+f'temp_file_{i}.txt'
#     tf2=wkdir+f'temp_file_{i}.pyw'
#
#     with open(tf1,'w') as f:
#         for j in range(len(hipair)):
#             f.write(f'{hipair[j][0]} {hipair[j][1]}\n')
#
#     # result=subprocess.run(['julia',wkdir+'max_entr_cont.jl',dir,files[i],str(dt),str(r0),'cos/sqr',str(3),tf1,tf2],capture_output=True,text=True)
#     # print("STDOUT:", result.stdout)
#     # print("STDERR:", result.stderr)
#
#     # result=subprocess.run(['matlab','-nodisplay','-nosplash','-r',f"addpath('{wkdir}'); max_entr_cont('{dir}','{files[i]}','{str(dt)}','{str(r0)}','cos/sqr','{str(3)}','{tf1}','{tf2}'); exit;"],capture_output=True,text=True)
#     # print("STDOUT:", result.stdout)
#     # print("STDERR:", result.stderr)
#     # matlab_command = f"/hpc2ssd/softwares/matlab/bin/matlab -r 'addpath(\"{wkdir}\"); max_entr_cont(\"{dir}\",\"{file}\",\"{dt}\",\"{r0}\",\"cos/sqr\",\"3\",\"{tf1}\",\"{tf2}\"); exit;'"
#     matlab_command = f"/hpc2ssd/softwares/matlab/bin/matlab -r \"addpath(\'{wkdir}\'); max_entr_cont(\'{dir}\',\'{file}\',\'{dt}\',\'{r0}\',\'cos/sqr\',\'3\',\'{tf1}\',\'{tf2}\'); exit;\""
#     subprocess.check_output(matlab_command, shell=True, text=True)
#
#     while not os.path.exists(tf2):
#         time.sleep(1)
#
#     while True:
#         try:
#             with open(tf2, 'r') as f:
#                 print('OK')
#             break
#         except IOError:
#             time.sleep(1)
#
#     P=pyw(tf2,'P:')
#     c=pyw(tf2,'C:')
#     nhipair=len(hipair)
#     C=np.zeros((nhipair,nhipair))
#     for i in range(nhipair):
#         for j in range(nhipair):
#             C[i,j]=c[i*nhipair+j]
#     os.remove(tf1)
#     os.remove(tf2)
#     now = datetime.datetime.now()
#     print(f'Complete the calculation of P and C ({i}):', now)
#
#     q1.put(P)
#     q2.put(C)
#
# def max_entr_b(sysname,imapfile,dir,checkfile,matfile,omapfile,pmin,dt,r0,tol):
#     """Read .trj files and calculate B matrix according to maximum entropy principle."""
#     """
#         Inputs:
#         imapfile is the inputting experimental Hi-C map file.
#         path is the inputting path to .trj files.
#         checkfile is the outputting check.pyw file.
#         matfile is the outputting B file.
#         omapfile is the outputting simulated Hi-C map file.
#         pmin is the the minimum contact probability of pairs as restraint of maximun entropy iteration.
#         dt is the simulation step length in .cfg file.
#         rc is the critical distance of the tanh data-driven term.
#         tol is the judgement of whether the simulated contact map is closed enough to the experimental
#         Hi-C map.
#     """
#
#     from pyw import pyw
#
#     now=datetime.datetime.now()
#     print(f'Begin to updata alphas:',now)
#
#     """Find the High-frequency Contacted Pairs"""
#     dirls=imapfile.split('/')
#     wkdir=''
#     for i in range(len(dirls)-1):
#         wkdir+=dirls[i]
#         wkdir+='/'
#     Pa=pyw(imapfile,'Map')[0]
#     ### Even though there is only one column of data, Pa is still a hierarchical list. So "[0]" is
#     ### essential.
#     N=int((len(Pa))**0.5)
#     P=np.zeros((N,N),dtype=float)
#     for i in range(N):
#         for j in range(N):
#             P[i,j]=Pa[i*N+j]
#     hipair=[]
#     for i in range(N):
#         for j in range(i):
#             if P[i,j]>=pmin and abs(i-j)>2:
#                 hipair.append([i,j])
#     nhipair=len(hipair)
#
#     """Enter Experimental Hi-C Data of High-frequency Contacted Pairs into List"""
#     # f_exp_ave=[]
#     # for i in range(len(hipair)):
#     #     f_exp_ave+=[P_exp[hipair[i][0],hipair[i][1]]]
#     f_exp_ave=[]
#     for i in range(N):
#         for j in range(N):
#             if [i,j] in hipair:
#                 f_exp_ave+=[P[i,j]]
#     ### These two strategies are equivalent, but the latter is faster.
#
#     files=file_par(dir)[0]
#     now = datetime.datetime.now()
#     print(f'Begin to open multiprocess:', now)
#     with multiprocessing.Manager() as manager:
#         q1=manager.Queue()
#         q2=manager.Queue()
#         procs=[]
#         for i in range(len(files)):
#             p=multiprocessing.Process(target=PC,args=(i,q1,q2,dir,files[i],dt,r0,hipair,wkdir))
#             p.start()
#             procs.append(p)
#         for p in procs:
#             p.join()
#         f_sim_ave=np.zeros((1,nhipair),dtype=float)
#         C=np.zeros((nhipair,nhipair),dtype=float)
#         while not q1.empty() and not q2.empty():
#             f_sim_ave+=q1.get()
#             C+=q2.get()
#         f_sim_ave/=len(files)
#         C/=len(files)
#
#     f_exp_ave=np.array(f_exp_ave)
#     delta_f_ave=f_sim_ave-f_exp_ave
#     delta_f_ave_sum=0
#     f_exp_ave_sum=0
#     for i in range(nhipair):
#         delta_f_ave_sum+=abs(delta_f_ave[0,i])
#         f_exp_ave_sum+=f_exp_ave[i]
#     rela_err=delta_f_ave_sum/f_exp_ave_sum
#     if rela_err<tol:
#         check=open(checkfile,'a')
#         ### If 'pyw' is defined as a variable here, the following pyw function will encounter error.
#         check.write('End\n')
#         check.close()
#     else:
#         check=open(checkfile,'a')
#         check.write(f'{rela_err}\n')
#         check.close()
#
#     """Update Alphas"""
#     D=np.zeros((nhipair,nhipair),dtype=float)
#     for i in range(nhipair):
#         for j in range(nhipair):
#             D[i,j]=f_sim_ave[0,i]*f_sim_ave[0,j]
#     # B=np.cov(f_sim,rowvar=False)
#     B=C-D
#     now=datetime.datetime.now()
#     print(f'Complete B',now)
#     # temp_filename1=tempfile.mktemp(suffix=".txt")
#     # np.savetxt(temp_filename1,B)
#     # temp_filename2=tempfile.mktemp(suffix=".txt")
#     # matlab_command=f"matlab -r 'data = load(\"{temp_filename1}\"); inverse_a = pinv(data); save(\"{temp_filename2}\", \"inverse_a\", \"-ascii\"); exit;'"
#     # subprocess.check_output(matlab_command, shell=True, text=True)
#     # B_=np.loadtxt(temp_filename2)
#     # os.remove(temp_filename1)
#     # os.remove(temp_filename2)
#     # # B_=np.linalg.pinv(B)
#     # now=datetime.datetime.now()
#     # print(f'Complete B^T',now)
#     # delta_f_ave=delta_f_ave.reshape(-1,1)
#     # ### For array format, neither "#.T" nor "np.transpose(#)" can turn it into a column vector. For
#     # ### this format, they all can do.
#     # delta_alf=kB*T*np.matmul(B_,delta_f_ave)
#     # delta_alf=delta_alf.reshape(1,-1)
#     # alf0_all=pyw(ialffile,'alpha')
#     # alf0=[]
#     # for i in range(len(alf0_all[0])):
#     #     j,k=divmod(i,N)
#     #     if [j,k] in hipair:
#     #         alf0+=[alf0_all[0][i]]
#     # alf0=np.array(alf0)
#     # alf=alf0+delta_alf
#     # alf=alf.tolist()
#     # Alf=np.zeros((N,N),dtype=float)
#     # Cont=np.zeros((N,N),dtype=float)
#     # # Alf+=2
#     # for i in range(N):
#     #     for j in range(N):
#     #         if [i,j] in hipair:
#     #             k=hipair.index([i,j])
#     #             # Alf[i,j]=alf[0][k]
#     #             # Alf[j,i]=alf[0][k]
#     #             if alf[0][k]<=10 and alf[0][k]>=-10:
#     #                 Alf[i,j]=alf[0][k]
#     #                 Alf[j,i]=alf[0][k]
#     #             if alf[0][k]<-10:
#     #                 Alf[i,j]=-10
#     #                 Alf[j,i]=-10
#     #             if alf[0][k]>10:
#     #                 Alf[i,j]=10
#     #                 Alf[j,i]=10
#     #             Cont[i,j]=f_sim_ave[0][k]
#     #             Cont[j,i]=f_sim_ave[0][k]
#     #
#     #
#     # now=datetime.datetime.now()
#     # print(f'Alphas have been updated:',now)
#
#     Cont=np.zeros((N,N),dtype=float)
#     # Dlt_cont=np.zeros((N,N),dtype=float)
#     for i in range(N):
#         for j in range(N):
#             if [i,j] in hipair:
#                 k=hipair.index([i,j])
#                 Cont[i,j]=f_sim_ave[0][k]
#                 Cont[j,i]=f_sim_ave[0][k]
#                 # Dlt_cont[i,j]=delta_f_ave[0][k]
#                 # Dlt_cont[j,i]=delta_f_ave[0][k]
#
#     """Calculate Contact Probability as a Function of Sequence Distance"""
#     d=np.linspace(0,N-1,N)
#     cd_exp=cont_dist(P)
#     cd_sim=cont_dist(Cont)
#
#     """Calculate the Compartment Profile"""
#     c_exp=comprtmt(P)[0]
#     c_sim=comprtmt(Cont)[0]
#
#     """Calculate the Insulation Score"""
#     is_exp=insul_score(P,bin_size=10**5)[0]
#     is_sim=insul_score(Cont,bin_size=10**5)[0]
#
#     """Output Data"""
#     scipy.io.savemat(f'{matfile}',{'B':B})
#     pyw=open(omapfile,'w')
#     pyw.write('Contact probability:\n')
#     for i in range(N):
#         for j in range(N):
#             pyw.write(f'{Cont[i,j]}\n')
#     pyw.write('\n')
#     pyw.write('Simulated contact probability minus experimental one:\n')
#     for i in range(nhipair):
#         pyw.write(f'{delta_f_ave[0][i]}\n')
#     pyw.write('\n')
#     pyw.write('Contact probability in each sequence distance, Compartment signal, Insulation score:\n')
#     for i in range(N):
#         pyw.write(f'{cd_sim[i]}\t{c_sim[i]}\t{is_sim[i]}\n')
#
#     """Plot"""
#     opathlist=omapfile.split('/')[:len(omapfile.split('/'))-1]
#     opath=''
#     for i in range(len(opathlist)):
#         opath+=opathlist[i]+'/'
#     # ax=sns.heatmap(Cont,vmin=0,vmax=1,cmap='seismic',center=0)
#     # ax.invert_yaxis()
#     # plt.xlabel('Index of particles')
#     # plt.ylabel('Index of particles')
#     # ax.minorticks_on()
#     # ax.tick_params(axis="both", which="major", direction="out")
#     # ax.tick_params(axis="both", which="minor", direction="out")
#     # cbar = ax.collections[0].colorbar
#     # cbar.ax.tick_params(direction="out",labelcolor='black')
#     # plt.tight_layout()
#     # plt.savefig(opath+f'cont_map_{sysname}.png',format='png')
#     # plt.savefig(opath+f'cont_map_{sysname}.pdf',format='pdf')
#
#     Compare=np.zeros((N,N))
#     for i in range(N):
#         for j in range(N):
#             if i>j:
#                 Compare[i,j]=P[i,j]
#             else:
#                 Compare[i,j]=Cont[i,j]
#
#     plt.figure()
#     plt.imshow(Compare)
#     plt.gca().invert_yaxis()
#     plt.title('Upper right: Exp, Lower left: Sim')
#     plt.xlabel('Index of particles')
#     plt.ylabel('Index of particles')
#     plt.colorbar()
#     plt.savefig(opath+f'cont_map_{sysname}.png',format='png')
#     plt.savefig(opath+f'cont_map_{sysname}.pdf',format='pdf')
#
#     plt.figure()
#     plt.plot(d,cd_exp,label='Exp')
#     plt.plot(d,cd_sim,label='Sim')
#     plt.legend(loc='best')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.xlabel('Sequence distance')
#     plt.ylabel('Contact probability')
#     plt.savefig(opath+f'cont_dist_{sysname}.png',format='png')
#     plt.savefig(opath+f'cont_dist_{sysname}.pdf',format='pdf')
#
#     plt.figure()
#     plt.plot(d+1,c_exp,label='Exp')
#     plt.plot(d+1,c_sim,label='Sim')
#     plt.legend(loc='best')
#     plt.xlabel('Index of particles')
#     plt.ylabel('Compartment signal')
#     plt.savefig(opath+f'comprtmt_{sysname}.png',format='png')
#     plt.savefig(opath+f'comprtmt_{sysname}.pdf',format='pdf')
#
#     plt.figure()
#     plt.plot(d+1,is_exp,label='Exp')
#     plt.plot(d+1,is_sim,label='Sim')
#     plt.legend(loc='best')
#     plt.xlabel('Index of particles')
#     plt.ylabel('Insulation score')
#     plt.savefig(opath+f'insul_score_{sysname}.png',format='png')
#     plt.savefig(opath+f'insul_score_{sysname}.pdf',format='pdf')


# import datetime
# n=10000
# t1=[]
# for m in range(10):
#     start=datetime.datetime.now()
#     for i in range(n):
#         for j in range(n):
#             x=i+j
#             # y=i*j
#             # z=i-j
#     end=datetime.datetime.now()
#     t1+=[(end-start).total_seconds()]
#
# t2=[]
# for m in range(10):
#     start=datetime.datetime.now()
#     for i in range(n):
#         for j in range(n):
#             x=i+j
#             y=i*j
#             # z=i-j
#     end=datetime.datetime.now()
#     t2+=[(end-start).total_seconds()]
#
# t3=[]
# for m in range(10):
#     start=datetime.datetime.now()
#     for i in range(n):
#         for j in range(n):
#             x=i+j
#             y=i*j
#             z=i-j
#     end=datetime.datetime.now()
#     t3+=[(end-start).total_seconds()]
# time=[1,2,3,4,5,6,7,8,9,10]
# import matplotlib.pyplot as plt
# plt.plot(time,t1,label='t1')
# plt.plot(time,t2,label='t2')
# plt.plot(time,t3,label='t3')
# plt.legend(loc='best')
# plt.show()

# import numpy as np
# import datetime
#
# def f(d,d0=0.6,beta=3):
#     if d < d0 - 1 / beta:
#         f = 1
#     if d >= d0 - 1 / beta and d < d0 + 1 / beta:
#         f = (np.cos(np.pi * beta * (d - d0 + 1 / beta) / 4)) ** 2
#     if d >= d0 + 1 / beta:
#         f = 0
#     return f
#
# F = np.vectorize(f)
#
# time=[]
# t1=[]
# t2=[]
# t3=[]
#
# for n in range(1000,5000,1000):
#     print(n)
#     time+=[n]
#     M = np.random.random((n, n))
#     L=np.random.random((n, n))
#
#     start=datetime.datetime.now()
#     N=F(M)
#     end=datetime.datetime.now()
#     t1+=[(end-start).total_seconds()]
#
#     start=datetime.datetime.now()
#     for i in range(np.shape(M)[0]):
#         for j in range(np.shape(M)[1]):
#             L[i,j]=f(M[i,j])
#     end=datetime.datetime.now()
#     t2+=[(end-start).total_seconds()]
#
#     start=datetime.datetime.now()
#     N=np.tanh(M)
#     end=datetime.datetime.now()
#     t3+=[(end-start).total_seconds()]
#
# import matplotlib.pyplot as plt
# plt.plot(time,t1,label='t1')
# plt.plot(time,t2,label='t2')
# plt.plot(time,t3,label='t3')
# plt.legend(loc='best')
# plt.show()

import numpy as np
def f(x,y,i):
    u=x+y*i
    v=x-y*i
    return u,v

x=np.random.random((3,3))
y=np.random.random((3,3))

a=np.zeros((3,3))
b=np.zeros((3,3))
for i in range(10):
    a,b+=f(x,y,i)