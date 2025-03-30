"""Import Modules"""
import numpy as np
import sys
sys.path.append('D:\\Biophysics\\MyPython\\functions')
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from pyw import pyw

"""Set Arguments"""
file='shp_seg.pyw'
outname=file[:-4]

"""Read Data"""
q=pyw(file,dir=True)

"""Output Data"""
np.save(f'{outname}.npy',q)