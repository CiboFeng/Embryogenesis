#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition a128m512u
#SBATCH --output ./job_out/%A_%a.out
#SBATCH --error ./job_out/%A_%a.err
#SBATCH --array 0-9

i=$SLURM_ARRAY_TASK_ID
module load mpi/openmpi-4.1.5
export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:$MODULEPATH
module load pl25mod
module load matlab/2023b

mkdir num
echo "0" | trjconv_mpi -f ../../max_entr/zygote_m/$i/trajout.xtc -s Polymer_600_0.4.pdb -o num/last_$i.gro -fr frame.ndx
grompp_mpi -f 120.28.mdp -c num/last_$i.gro -p chr14_a_11.top -o num/md_$i.tpr -po num/mdout_$i.mdp
mpirun -np 4 mdrun_mpi -s num/md_$i.tpr -x num/trj_$i.xtc -e num/ener_$i.edr -c num/confout_$i.gro -g num/md_$i.log -cpo num/state_$i.cpt -tableb Tables/table.xvg -table Tables/table_PM.xvg -tablep Tables/table_ME.xvg -plumed plumed_SR.dat -noddcheck -pd
echo "0" | trjconv_mpi -f num/trj_$i.xtc -s num/md_$i.tpr -o num/trj_$i.pdb
dl *colvar
