#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N b_31_0.10_120.28_9_M_Smooth
#$ -j y
#$ -o RUN.out
#$ -q cpu_short
#$ -P otherprj

echo "0" | trjconv_4.5.7 -f traj.xtc -s md.tpr -o tmp.pdb -b 1000
/fcal.x tmp.pdb /Contact_a/Probability.dat 812
rm tmp.pdb
rm \#*
