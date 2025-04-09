import numpy as np
data = data = np.loadtxt('Po_e.dat')
from iced import normalization
from iced import filter
#counts = filter.filter_low_counts(data, percentage=0.04)
normed = normalization.ICE_normalization(data)
np.savetxt('normal.dat',normed)
