#!/bin/bash

sysname=400
itrs=1
itre=2
njob=200
nthermo=0
Thi=3
Tlo=1
nstep_annl=500000
nstep_relx=500000
ndump=2000
nstep_equl=1000000
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
#SBATCH --output $pwd/jobout/%j_init.out
#SBATCH --error $pwd/jobout/%j_init.err

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
#SBATCH --output $pwd/jobout/%j_init.out
#SBATCH --error $pwd/jobout/%j_init.err

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
#SBATCH --output $pwd/jobout/%j_file.out
#SBATCH --error $pwd/jobout/%j_file.err

python $pwd/max_entr_file.py -s $sysname -w $pwd -a $pwd/alf_$i/alf_$sysname.pyw -c $pwd/$sysname.cfg -p $pwd/lmp_$i/${sysname}_i.par -thermo $nthermo -Th $Thi -Tl $Tlo -na $nstep_annl -nr $nstep_relx -dump $ndump -ne $nstep_equl" > $pwd/alf_$i/max_entr_file.py.sh
    jobid=$(sbatch --dependency=afterok:$jobid --parsable $pwd/alf_$i/max_entr_file.py.sh)
	echo max_entr_file.py.sh: $jobid

	jobids=""
	j=1
	while  [ $j -le $njob ]
	do

		touch $pwd/lmp_$i/${sysname}_$j.sh
		echo -e "#!/bin/bash
#SBATCH --job-name lmp
#SBATCH --ntasks 5
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $pwd/jobout/%j_lmp.out
#SBATCH --error $pwd/jobout/%j_lmp.err

module load mpi/mpich-4.1.2
export PATH=/hpc2hdd/home/cfeng593/opt/lammps/lammps-9Oct20/build:\$PATH

cp $sysname.cfg $pwd/lmp_$i/${sysname}_$j.cfg
sed -i '/log/c\\\log\t\t\t$pwd/lmp_$i/${sysname}_$j.log' $pwd/lmp_$i/${sysname}_$j.cfg
sed -i '/read_data/c\\\read_data\t\t$pwd/lmp_$i/${sysname}_$j.dat' $pwd/lmp_$i/${sysname}_$j.cfg
sed -i '/trj/c\\\dump\t\t\tdump all custom $ndump $pwd/lmp_$i/${sysname}_$j.trj id type x y z' $pwd/lmp_$i/${sysname}_$j.cfg
sed -i '/write_data/c\\\write_data\t\t$pwd/lmp_$(($i+1))/${sysname}_$j.dat' $pwd/lmp_$i/${sysname}_$j.cfg

mpirun -n 4 lmp -in $pwd/lmp_$i/${sysname}_$j.cfg" > $pwd/lmp_$i/${sysname}_$j.sh
		jobid_temp=$(sbatch --dependency=afterok:$jobid --parsable $pwd/lmp_$i/${sysname}_$j.sh)
		echo ${sysname}_$j.sh: $jobid_temp
		
		touch $pwd/pc_$i/max_entr_pc_$j.py.sh
		echo -e "#!/bin/bash
#SBATCH --job-name pc
#SBATCH --ntasks 2
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $pwd/jobout/%j_pc.out
#SBATCH --error $pwd/jobout/%j_pc.err

python $pwd/max_entr_pc.py -s $sysname -t $pwd/lmp_$i/${sysname}_$j.trj -w $pwd -pc $pwd/pc_$i/pc_${sysname}_$j.npz" > $pwd/pc_$i/max_entr_pc_$j.py.sh
		jobid_temp=$(sbatch --dependency=afterany:$jobid_temp --parsable $pwd/pc_$i/max_entr_pc_$j.py.sh)
		echo max_entr_pc_$j.py.sh: $jobid_temp
		jobids="$jobids,$jobid_temp"

		((j++))
	done
	jobid=${jobids#*,}

	touch $pwd/alf_$i/max_entr_alf.py.sh
	echo -e "#!/bin/bash
#SBATCH --job-name alf
#SBATCH --ntasks 2
#SBATCH --exclude cpu1-1
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $pwd/jobout/%j_alf.out
#SBATCH --error $pwd/jobout/%j_alf.err

python $pwd/max_entr_alf.py -s $sysname -d $pwd/pc_$i -w $pwd -i $pwd/alf_$i/alf_$sysname.pyw -o $pwd/alf_$(($i+1))/alf_$sysname.pyw" > $pwd/alf_$i/max_entr_alf.py.sh
	jobid=$(sbatch --dependency=afterok:$jobid --parsable $pwd/alf_$i/max_entr_alf.py.sh)
    echo max_entr_alf.py.sh: $jobid

done
