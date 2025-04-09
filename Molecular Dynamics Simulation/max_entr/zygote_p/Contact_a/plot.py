import argparse
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
parser=argparse.ArgumentParser()
parser.add_argument("-f1",type=str,help="The first alpha file")
parser.add_argument("-f2",type=str,help="The second alpha file")
args=parser.parse_args()
file1=args.f1
file2=args.f2
f=open(file1)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
ls_f=[list(map(float,sublist)) for sublist in ls_f]
alf1=[]
for i in range(len(ls_f)):
    alf1+=[ls_f[i][2]]
f=open(file2)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
ls_f=[list(map(float,sublist)) for sublist in ls_f]
alf2=[]
for i in range(len(ls_f)):
    alf2+=[ls_f[i][2]]
plt.figure()
plt.plot(alf1,alf2,'.')
plt.xlabel(file1)
plt.ylabel(file2)
matplotlib.use('QtAgg')
plt.show()
