#!/bin/bash
#SBATCH --job-name cfeng593
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output ./job_out/%A_%a.out
#SBATCH --error ./job_out/%A_%a.err
#SBATCH --array 0-1213

i=$SLURM_ARRAY_TASK_ID
module load mpi/openmpi-4.1.5
export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:$MODULEPATH
module load pl25mod
module load matlab/2023b

mkdir -p num/equ_$i
mkdir -p num/swch_$i
mkdir -p num/swch_detl_$i
mkdir -p num/merg_$i
mkdir -p num/merg_rdc_$i
frame=$(awk 'NR=='$(($i+2))'' frame_clst.ndx)
touch num/equ_$i/frame.ndx
echo -e "[ frames ]
$frame" > num/equ_$i/frame.ndx

echo "0" | trjconv_mpi -f ../max_entr/zygote_p/trj_merg.pdb -s ../diff/zp/Polymer_600_0.4.pdb -o num/equ_$i/init.gro -fr num/equ_$i/frame.ndx
grompp_mpi -f equ.mdp -c num/equ_$i/init.gro -p ../diff/zp/chr14_a_13.top -o num/equ_$i/md.tpr -po num/equ_$i/mdout.mdp
mpirun -np 4 mdrun_mpi -s num/equ_$i/md.tpr -x num/equ_$i/trj.xtc -e num/equ_$i/ener.edr -c num/equ_$i/confout.gro -g num/equ_$i/md.log -cpo num/equ_$i/state.cpt -tableb ../diff/zp/Tables/table.xvg -table ../diff/zp/Tables/table_PM.xvg -tablep ../diff/zp/Tables/table_ME.xvg -plumed pld_zp.dat -noddcheck -pd

grompp_mpi -f swch.mdp -c num/equ_$i/confout.gro -t num/equ_$i/state.cpt -p ../clst/esc/chr14_a_9.top -o num/swch_$i/md.tpr -po num/swch_$i/mdout.mdp
mpirun -np 4 mdrun_mpi -s num/swch_$i/md.tpr -x num/swch_$i/trj.xtc -e num/swch_$i/ener.edr -c num/swch_$i/confout.gro -g num/swch_$i/md.log -cpo num/swch_$i/state.cpt -tableb ../clst/esc/Tables/table.xvg -table ../clst/esc/Tables/table_PM.xvg -tablep ../clst/esc/Tables/table_ME.xvg -plumed pld_esc.dat -noddcheck -pd
echo "0" | trjconv_mpi -f num/swch_$i/trj.xtc -s num/swch_$i/md.tpr -o num/swch_$i/trj.pdb -fr frame_cut.ndx

grompp_mpi -f swch_detl.mdp -c num/equ_$i/confout.gro -t num/equ_$i/state.cpt -p ../clst/esc/chr14_a_9.top -o num/swch_detl_$i/md.tpr -po num/swch_detl_$i/mdout.mdp
mpirun -np 4 mdrun_mpi -s num/swch_detl_$i/md.tpr -x num/swch_detl_$i/trj.xtc -e num/swch_detl_$i/ener.edr -c num/swch_detl_$i/confout.gro -g num/swch_detl_$i/md.log -cpo num/swch_detl_$i/state.cpt -tableb ../clst/esc/Tables/table.xvg -table ../clst/esc/Tables/table_PM.xvg -tablep ../clst/esc/Tables/table_ME.xvg -plumed pld_esc.dat -noddcheck -pd
echo "0" | trjconv_mpi -f num/swch_detl_$i/trj.xtc -s num/swch_detl_$i/md.tpr -o num/swch_detl_$i/trj.pdb

trjcat_mpi -f num/swch_detl_$i/trj.pdb num/swch_$i/trj.pdb -cat -o num/merg_$i/trj.xtc
echo "0" | trjconv_mpi -f num/merg_$i/trj.xtc -s ../diff/zp/Polymer_600_0.4.pdb -o num/merg_rdc_$i/trj.xtc -fr frame_rdc.ndx
