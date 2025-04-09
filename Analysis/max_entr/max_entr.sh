#!/bin/bash

sysname=$(echo $(find . -type f -name "exp_cont_prob*") | cut -c 17- | cut -d '.' -f 1)
itrs=1
itre=5
njob=100
nthermo=0
Thi=3
Tlo=1
nstep_annl=500000
nstep_relx=500000
ndump=2000
nstep_equl=2000000
#nthermo=10
#nstep_annl=2000
#nstep_relx=1000
#ndump=100
#nstep_equl=10000
pwd=$(pwd)

touch $pwd/init.sh
echo -e "#!/bin/bash
#SBATCH --job-name init
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition i64m512u
#SBATCH --output $pwd/job_out/%j_init.out
#SBATCH --error $pwd/job_out/%j_init.err

echo Begin" > $pwd/init.sh
jobid=$(sbatch --parsable $pwd/init.sh)
echo init.sh: $jobid

for i in `seq $itrs 1 $itre`
do

	if [ $i == 1 ] 
	then
		mkdir lmp_1
		mkdir alf_1
		touch $pwd/alf_1/max_entr_init.py.sh
		echo -e "#!/bin/bash
#SBATCH --job-name init
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition i64m512u
#SBATCH --output $pwd/job_out/%j_init.out
#SBATCH --error $pwd/job_out/%j_init.err

python $pwd/max_entr_init.py -s $sysname -w $pwd -a $pwd/alf_1/alf_$sysname.pyw -d $pwd/lmp_1 -n $njob" > $pwd/alf_1/max_entr_init.py.sh
    	jobid=$(sbatch --dependency=afterok:$jobid --parsable $pwd/alf_1/max_entr_init.py.sh)
		echo max_entr_init.py.sh: $jobid
	fi

	mkdir lmp_$(($i+1))
	mkdir pc_$i
	mkdir alf_$(($i+1))

	touch $pwd/alf_$i/max_entr_file.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name file
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --partition i64m512u
#SBATCH --output $pwd/job_out/%j_file.out
#SBATCH --error $pwd/job_out/%j_file.err

python $pwd/max_entr_file.py -s $sysname -w $pwd -a $pwd/alf_$i/alf_$sysname.pyw -c $pwd/$sysname.cfg -p $pwd/lmp_$i/${sysname}.par -thermo $nthermo -Th $Thi -Tl $Tlo -na $nstep_annl -nr $nstep_relx -dump $ndump -ne $nstep_equl" > $pwd/alf_$i/max_entr_file.py.sh
    jobid=$(sbatch --dependency=afterok:$jobid --parsable $pwd/alf_$i/max_entr_file.py.sh)
	echo max_entr_file.py.sh: $jobid

	touch $pwd/lmp_$i/$sysname.sh
	echo -e "#!/bin/bash
#SBATCH --job-name lmp
#SBATCH --ntasks 4
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $pwd/job_out/%A_%a_lmp.out
#SBATCH --error $pwd/job_out/%A_%a_lmp.err
#SBATCH --array 1-$njob

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/lammps-9Oct20/build:\$PATH
export PATH=/hpc2hdd/home/chu-amat/cbfengphy/opt/lammps/lammps-9Oct20/build:\$PATH

j=\$SLURM_ARRAY_TASK_ID
cp $sysname.cfg $pwd/lmp_$i/${sysname}_\$j.cfg
sed -i '/log/c\\\log\t\t\t$pwd/lmp_$i/${sysname}_'\$j'.log' $pwd/lmp_$i/${sysname}_\$j.cfg
sed -i '/read_data/c\\\read_data\t\t$pwd/lmp_$i/${sysname}_'\$j'.dat' $pwd/lmp_$i/${sysname}_\$j.cfg
sed -i '/trj/c\\\dump\t\t\tdump all custom $ndump $pwd/lmp_$i/${sysname}_'\$j'.trj id type x y z' $pwd/lmp_$i/${sysname}_\$j.cfg
sed -i '/write_data/c\\\write_data\t\t$pwd/lmp_$(($i+1))/${sysname}_'\$j'.dat' $pwd/lmp_$i/${sysname}_\$j.cfg

mpirun -n 4 lmp -in $pwd/lmp_$i/${sysname}_\$j.cfg" > $pwd/lmp_$i/$sysname.sh
	jobid=$(sbatch --dependency=afterok:$jobid --parsable $pwd/lmp_$i/$sysname.sh)
	echo $sysname.sh: $jobid
		
	touch $pwd/pc_$i/max_entr_pc.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name pc
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $pwd/job_out/%A_%a_pc.out
#SBATCH --error $pwd/job_out/%A_%a_pc.err
#SBATCH --array 1-$njob

j=\$SLURM_ARRAY_TASK_ID
python $pwd/max_entr_pc.py -s $sysname -w $pwd -t $pwd/lmp_$i/${sysname}_\$j.trj -pc $pwd/pc_$i/pc_${sysname}_\$j.npz" > $pwd/pc_$i/max_entr_pc.py.sh
	jobid=$(sbatch --dependency=afterany:$jobid --parsable $pwd/pc_$i/max_entr_pc.py.sh)
	echo max_entr_pc.py.sh: $jobid

	touch $pwd/alf_$i/max_entr_alf.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name alf
#SBATCH --ntasks 1
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $pwd/job_out/%j_alf.out
#SBATCH --error $pwd/job_out/%j_alf.err

python $pwd/max_entr_alf.py -s $sysname -w $pwd -d $pwd/pc_$i -i $pwd/alf_$i/alf_$sysname.pyw -o $pwd/alf_$(($i+1))/alf_$sysname.pyw" > $pwd/alf_$i/max_entr_alf.py.sh
	jobid=$(sbatch --dependency=afterany:$jobid --parsable $pwd/alf_$i/max_entr_alf.py.sh)
    echo max_entr_alf.py.sh: $jobid

done
