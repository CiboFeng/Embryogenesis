#!/bin/bash
export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:$PATH
export PATH=/hpc2hdd/home/chu-amat/miniconda3/bin:$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:$MODULEPATH
module load pl25mod
module load matlab/2023b
#module load gnuplot-5.0.5
source /hpc2hdd/home/chu-amat/miniconda3/etc/profile.d/conda.sh
conda activate ice

RUNROOT=`pwd`
cd ../
WORKDIR=`pwd`
cd $RUNROOT

$RUNROOT/CP_GD.x $WORKDIR/Probability.dat 600
mv CP_GD.dat CP_GD_Sim.dat

#octave Pcompare.m
/hpc2ssd/softwares/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Pcompare.m');exit;"

./Obs_Exp_P.x $WORKDIR/Probability.dat 600
python ice_do.py
#octave Po_e.m
/hpc2ssd/softwares/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Po_e.m');exit;"
./InsulationScore.x $WORKDIR/Probability.dat 600

mv vector.dat CM_Sim.dat
mv InsulationScore.dat InsulationScore_Sim.dat

#gnuplot PLOT.bat 
python plot.py
