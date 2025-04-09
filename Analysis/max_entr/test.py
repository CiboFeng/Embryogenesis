# import numpy as np
# import datetime
# x=200
# y=200
# z=2000
# f=np.random.rand(x,y,z)
# hipair=[]
# for i in range(int(x*y/10)):
#     hipair.append([])
#     hipair[i]+=[divmod(i,x)[0],divmod(i,x)[1]]
# hipairm=np.array(hipair)
#
# now1=datetime.datetime.now()
# f_sim=[]
# #l=-1
# for i in range(x):
#     for j in range(y):
#         if [i,j] in hipair:
#             #l+=1
#             l=hipair.index([i,j])
#             f_sim.append([])
#             for k in range(z):
#                 f_sim[l]+=[f[i,j,k]]
# f_sim=np.array(f_sim).T #If the "if" is only executed in two layer of loops instead of three layer, the f_sim will not be the normal format. So it should be transposed.
# now2=datetime.datetime.now()
#
# f_sim=[]
# for i in range(len(hipair)):
#     f_sim.append([])
#     for j in range(z):
#         f_sim[i]+=[f[hipair[i][0],hipair[i][1],j]]
# f_sim=np.array(f_sim).T
# now3=datetime.datetime.now()
# print((now3-now2)/(now2-now1))

import numpy as np
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

x=np.array([1,2,3,4,5,6,7,8,9,10])
y=0.5*x+0+np.random.normal(0,0.2,x.shape)

r=np.corrcoef(x,y)[0,1]
r2=r2_score(y,x)
print(r)
print(r2)

x=10**np.linspace(-4,0,200)
y=x+np.random.normal(0,0.1,x.shape)
plt.subplot(121)
plt.plot(x,y,'.')
plt.subplot(122)
plt.plot(x,y,'.')
plt.xscale('log')
plt.yscale('log')
plt.show()