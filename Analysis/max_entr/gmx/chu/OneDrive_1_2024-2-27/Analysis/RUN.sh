#!/bin/bash
RUNROOT=`pwd`
cd ../
WORKDIR=`pwd`
cd $RUNROOT

$RUNROOT/CP_GD.x $WORKDIR/Probability.dat 577
mv CP_GD.dat CP_GD_Sim.dat

octave Pcompare.m

./Obs_Exp_P.x $WORKDIR/Probability.dat 577
python ice_do.py
octave Po_e.m
./InsulationScore.x $WORKDIR/Probability.dat 577

mv vector.dat CM_Sim.dat
mv InsulationScore.dat InsulationScore_Sim.dat

gnuplot PLOT.bat 
