#!/bin/bash
export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:$PATH
export PATH=/hpc2hdd/home/chu-amat/miniconda3/bin:$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:$MODULEPATH
module load pl25mod
module load matlab/2023b
source /hpc2hdd/home/chu-amat/miniconda3/etc/profile.d/conda.sh
conda activate ice

./Obs_Exp_P.x ../PC_100000_iced_chr14_dense.matrix_P.dat 600
python ice_do.py
#octave Po_e.m
/hpc2ssd/softwares/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Po_e.m');exit;"
