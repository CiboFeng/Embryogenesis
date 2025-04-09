"""Import Modules"""
import argparse
import datetime
import numpy as np
import scipy.io as sio
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from max_entr_alf import max_entr_alf

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-w",type=str,help="The work directory")
parser.add_argument("-m",type=str,help="The inputting Binv file")
parser.add_argument("-d",type=str,help="The map of simulated contact probability minus experimental one")
parser.add_argument("-i",type=str,help="The inputting alpha.pyw file to be updated")
parser.add_argument("-o",type=str,help="The outputting updated alpha.pyw file")
parser.add_argument("-T",type=str,help="The simulation temperature")
args=parser.parse_args()
wkdir=args.w
matfile=args.m
dmapfile=args.d
ialffile=args.i
oalffile=args.o
T=float(args.T)
emapfile=wkdir+'/tse_cont_map_1ubq_143.pyw'

sysname='1ubq'
# mapfile='D:\Biophysics\MyPython\inputs\\tse_cont_map_1ubq_143.pyw'
# path='D:\Biophysics\MyPython\inputs\\trj'
# # matfile='cont_samp.mat'
# checkfile='check.pyw'
# ialffile='D:\Biophysics\MyPython\\applications\exp_max_entr_file\\alpha_1ubq.pyw'
# oalffile='alpha.pyw'
# T=120.2717

if __name__=='__main__':
    max_entr_alf(sysname=sysname,
                 emapfile=emapfile,
                 matfile=matfile,
                 dmapfile=dmapfile,
                 ialffile=ialffile,
                 oalffile=oalffile,
                 pmin=0,
                 kB=1,
                 T=T)