#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output ./job_out/%A_%a.out
#SBATCH --error ./job_out/%A_%a.err
#SBATCH --array 0-55

i=$SLURM_ARRAY_TASK_ID
module load mpi/openmpi-4.1.5
export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:$MODULEPATH
module load pl25mod
module load matlab/2023b

mkdir -p num/reshp_$i
touch num/reshp_$i/frame.ndx
echo -e "[ frames ]
$(($i+1))" > num/reshp_$i/frame.ndx

for j in `seq 0 1 1213`
do
	echo "0" | trjconv_mpi -f num/merg_rdc_$j/trj.xtc -s ../diff/zp/Polymer_600_0.4.pdb -o num/reshp_$i/trj.xtc -fr num/reshp_$i/frame.ndx -app
done
echo "0" | trjconv_mpi -f num/reshp_$i/trj.xtc -s ../diff/zp/Polymer_600_0.4.pdb -o num/reshp_$i/trj.pdb
