import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw
import numpy as np

pywfile="D:\\Biophysics\\MyPython\\inputs\\tse_cont_map_1ubq_143.pyw"
matrixfile='PC_100000_iced_chr14_dense'

p=pyw(pywfile,'Map')[0]
nat=int(len(p)**0.5)
P=np.zeros((nat,nat))
for i in range(nat):
    for j in range(nat):
        P[i,j]=p[i*nat+j]

mat1=open(matrixfile+'.matrix_P','w')
for i in range(nat):
    for j in range(nat):
        mat1.write(f'{P[i,j]} ')
    mat1.write('\n')
mat1.close()

mat2=open(matrixfile+'.matrix_P.dat','w')
for i in range(nat):
    for j in range(nat):
        mat2.write(f'{i+1} {j+1} {P[i,j]}\n')
mat2.close()
