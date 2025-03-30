import numpy as np

npyfile='cont_prob.npy'
outname='cont_prob_rdc'
cl=['Z','ESC']
sw=['ZP-ESC','ZM-ESC']
cl_=[['zp','zm'],['esc','esc']]
sw_=['zp_esc','zm_esc']
ncl=len(cl)
nsw=len(sw)
nat=600
nfr=56
t=[0,0.01,
   0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,
   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
   2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,
   20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,
   200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,
   2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0,9000.0,10000.0]
t_rdc=[0,0.01,0.1,1.0,10.0,100.0,1000.0,10000.0]

P=np.load(npyfile)
P_rdc=np.zeros((nsw,len(t_rdc),nat,nat))
for i in range(len(t_rdc)):
    P_rdc[:,i,:,:]=P[:,t.index(t_rdc[i])+1,:,:]

np.save(f'{outname}.npy',P_rdc)