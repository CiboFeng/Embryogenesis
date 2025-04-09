import numpy as np
import sys
sys.path.append('D:\Biophysics\MyPython\\functions')
sys.path.append('/hpc2hdd/home/chu-amat/cbfengphy/functions')
from pyw import pyw

file='../clst2.pyw'
idx=pyw(file,f'Index of ZP:',fmt='int')[0]
idx=sorted(idx)
ndx=open('frame_clst.ndx','w')
ndx.write('[ frames ]\n')
for i in range(len(idx)):
    ndx.write(f'{idx[i]+1}\n')
ndx.close()

idx=list(range(102,10002))
ndx=open('frame_cut.ndx','w')
ndx.write('[ frames ]\n')
for i in range(len(idx)):
    ndx.write(f'{idx[i]}\n')
ndx.close()

t=np.round(np.arange(0.0,100.01,0.01),3).tolist()+np.round(np.arange(101.0,10001.0,1.0),3).tolist()
t_rdc=[0,0.01,
       0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
       0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
       2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
       20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
       200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
       2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
idx=[]
for i in range(len(t_rdc)):
    idx+=[t.index(t_rdc[i])]
ndx=open('frame_rdc.ndx','w')
ndx.write('[ frames ]\n')
for i in range(len(idx)):
    ndx.write(f'{idx[i]+1}\n')
ndx.close()
