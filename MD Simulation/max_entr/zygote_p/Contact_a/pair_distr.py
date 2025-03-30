import numpy as np
import matplotlib.pyplot as plt
file='Contact_a_0.dat'
f=open(file)
lsf=f.readlines()
f.close()
I=[]
J=[]
for i in range(len(lsf)):
    ls_f=lsf[i].strip('\n').split()
    I+=[int(ls_f[0])]
    J+=[int(ls_f[1])]
I=np.array(I)
J=np.array(J)
d=np.abs(I-J)
hist,edge=np.histogram(d,bins=np.arange(np.min(d),np.max(d)+1,5))
d=(edge[1:]+edge[:-1])/2
plt.plot(d,hist)
plt.show()
