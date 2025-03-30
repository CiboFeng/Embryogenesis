import argparse
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib import colors
parser=argparse.ArgumentParser()
parser.add_argument("-f",type=str,help="The matrix file")
args=parser.parse_args()
file=args.f
f=open(file)
lsf=f.readlines()
f.close()
ls_f=[]
for lf in lsf:
    l_f=lf.strip('\n').split()
    ls_f.append(l_f)
ls_f=[list(map(float,sublist)) for sublist in ls_f]
N=int(ls_f[len(ls_f)-1][1])
P=np.zeros((N,N),dtype=float)
for i in range(len(ls_f)):
    #if len(ls_f[i])>0: #In order to make sure the list index in next line does not exceed the range.
    P[int(ls_f[i][0])-1,int(ls_f[i][1])-1]=ls_f[i][2]
    P[int(ls_f[i][1])-1,int(ls_f[i][0])-1]=ls_f[i][2]

