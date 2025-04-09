import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
n=10
a=5
b=1 #小
x=np.linspace(1,n,n+n-1)
y=np.zeros(np.shape(x))
for i in range(0,len(x),2):
    y[i]=-x[i]*a
for i in range(1,len(x),2):
    y[i]=-2*x[i]*b
plt.plot(x,y)
plt.plot(x,y,'.')
ax=plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(1))
plt.grid(True)
plt.show()