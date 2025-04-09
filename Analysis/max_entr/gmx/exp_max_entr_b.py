"""Import Modules"""
import argparse
import sys
sys.path.append('/hpc2hdd/home/cfeng593/proj/mypython/functions')
from max_entr_b import max_entr_b

"""Set Arguments"""
parser=argparse.ArgumentParser()
parser.add_argument("-d",type=str,help="The directory of .trj files")
parser.add_argument("-w",type=str,help="The work directory")
parser.add_argument("-m",type=str,help="The .mat file")
parser.add_argument("-o",type=str,help="The outputting contact map file")
args=parser.parse_args()
dir=args.d
wkdir=args.w
matfile=args.m
omapfile=args.o
imapfile=wkdir+'/tse_cont_map_1ubq_143.pyw'
checkfile=wkdir+'/check_1ubq.pyw'

sysname='1ubq'
# imapfile='D:\Biophysics\MyPython\inputs\\tse_cont_map_1ubq_143.pyw'
# path='D:\Biophysics\MyPython\inputs\\trj'
# # matfile='cont_samp.mat'
# checkfile='check.pyw'
# ialffile='D:\Biophysics\MyPython\\applications\exp_max_entr_file\\alpha_1ubq.pyw'
# oalffile='alpha.pyw'
# T=120.2717

if __name__=='__main__':
    max_entr_b(sysname=sysname,
               imapfile=imapfile,
               dir=dir,
               checkfile=checkfile,
               matfile=matfile,
               omapfile=omapfile,
               pmin=0,
               dt=0.01,
               r0=2.5,
               tol=0.01)