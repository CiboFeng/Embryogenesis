#!/bin/bash
RUNROOT=`pwd`
protein=chr14
MDPfile="$RUNROOT/Polymer.mdp"
GROMACSDIR="/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin"
PDBDIR="$RUNROOT/PDB"
TOPDIR="$RUNROOT/TOP_ME"
ContactDIR="$RUNROOT/Contact_a"
Tablefile="$RUNROOT/Tables/table_PM.xvg"
Tablepfile="$RUNROOT/Tables/table_ME.xvg"

#Variables
#Lambda=0.01
#Lambda=0.05
#Lambda=0.20
#Lambda=0.50
#Lambda=1.00
Lambdas=(0.01 0.01 0.05 0.05 0.20 0.20 0.50 0.50 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00)
MAXNUM=9
BTIME=1000

mkdir -p $RUNROOT/ConstantT
cd $RUNROOT/ConstantT

Nres=600
sigma=0.4
NCPU=4
#Radius=`echo "$Nres $sigma" | awk '{print $2/2*(10*$1)**(1.0/3.0)-$2/2}'`
Radius=`echo "$Nres $sigma" | awk '{print $2/2*(10*$1)**(1.0/3.0)+4*$2}'`
sed "s/NRES/$Nres/g" $RUNROOT/PLUMED/plumed_SR.dat | sed "s/Radius/$Radius/g" > $RUNROOT/plumed_SR.dat
box=`echo "$Radius" | awk '{print $1*2}'`

cat > init.sh <<EOF
#!/bin/bash
#SBATCH --job-name init
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $RUNROOT/jobout/%j.out
#SBATCH --error $RUNROOT/jobout/%j.err
echo OK	
EOF
jobidf=$(sbatch --parsable init.sh)
#jobidf=4280378
echo init: $jobidf

#for alpha in `seq 0 1 1`
#for alpha in `seq 2 1 3`
#for alpha in `seq 4 1 5`
#for alpha in `seq 6 1 7`
#for alpha in `seq 8 1 10`
#for alpha in `seq 11 1 12`
#for alpha in `seq 13 1 15`
for alpha in `seq 0 1 20`
do
	Lambda=${Lambdas[$alpha]}
	alpha_old=$(($alpha-1))
	
	mkdir -p $RUNROOT/ConstantT/a_$alpha
	cd $RUNROOT/ConstantT/a_$alpha
	if [ "$alpha" == "0" ]; then
		cat $ContactDIR/Contact_a_${alpha}.dat $ContactDIR/Contact_a_19_old.dat > $ContactDIR/Contact_a_${alpha}.top.dat
		$TOPDIR/TOP.x $RUNROOT/PDBGeneration/Polymer_${Nres}_$sigma.pdb $ContactDIR/Contact_a_${alpha}.top.dat ${protein}_a_${alpha}
	else
                paste $RUNROOT/ConstantT/a_${alpha_old}/120.28/Probability/Contact_a_abs.dat $RUNROOT/Contact_a/Contact_a_${alpha_old}.dat | awk '{print $1,$2,$3*'$Lambda'+$6}' > $RUNROOT/Contact_a/Contact_a_${alpha}.dat
		cat $ContactDIR/Contact_a_${alpha}.dat $ContactDIR/Contact_a_19_old.dat > $ContactDIR/Contact_a_${alpha}.top.dat
                $TOPDIR/TOP.x $RUNROOT/PDBGeneration/Polymer_${Nres}_$sigma.pdb $ContactDIR/Contact_a_${alpha}.top.dat ${protein}_a_${alpha}
	fi
	TOPfile=$RUNROOT/ConstantT/a_$alpha/${protein}_a_${alpha}.top
	for Temp in 120.28
	do
		mkdir -p $RUNROOT/ConstantT/a_$alpha/$Temp
		cd $RUNROOT/ConstantT/a_$alpha/$Temp
		sed "s/TXTemp/$Temp/g" $MDPfile > $Temp.mdp
		rm list.dat
		for i in `seq 0 1 $MAXNUM`
		do
			mkdir $RUNROOT/ConstantT/a_${alpha}/$Temp/$i
			cd $RUNROOT/ConstantT/a_${alpha}/$Temp/$i
			j=`echo "$i $MAXNUM" | awk '{print $1*($2+1)}'`
			rm -r file 2>/dev/null
			jobids=""
			for k in `seq 0 1 9`
			do
cat > ConstantT.sh <<EOF
#!/bin/bash
#SBATCH --job-name a_${alpha}_${Temp}_${i}_PC_BigMDSteps_1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node $NCPU
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $RUNROOT/jobout/%j.out
#SBATCH --error $RUNROOT/jobout/%j.err

module load mpi/openmpi-4.1.5
export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:\$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:\$MODULEPATH
module load pl25mod
module load matlab/2023b

#for k in \`seq 0 1 4\`
#do
mkdir -p $RUNROOT/ConstantT/a_${alpha}/$Temp/$i/$k
cd $RUNROOT/ConstantT/a_${alpha}/$Temp/$i/$k
l=\$(($j+$k))
PDBfile=$PDBDIR/N\$l.gro	
#EM
$GROMACSDIR/editconf_mpi -f \$PDBfile -o N\${l}_center.gro -box $box $box $box -bt cubic -c
$GROMACSDIR/grompp_mpi -f $RUNROOT/Polymer_EM.mdp -c N\${l}_center.gro -p $TOPfile -o em.tpr > output 2>&1 
mpirun -np $NCPU $GROMACSDIR/mdrun_mpi -noddcheck -deffnm em -tableb $RUNROOT/Tables/table.xvg -table $Tablefile -tablep $Tablepfile -plumed $RUNROOT/plumed_SR.dat -pd
#MD
$GROMACSDIR/grompp_mpi -f ../../$Temp.mdp -c em.gro -p $TOPfile -o md.tpr > output 2>&1 
mpirun -np $NCPU $GROMACSDIR/mdrun_mpi -noddcheck -s md.tpr -tableb $RUNROOT/Tables/table.xvg -table $Tablefile -tablep $Tablepfile -plumed $RUNROOT/plumed_SR.dat -pd

echo "0" | $GROMACSDIR/trjconv_mpi -f traj.xtc -s md.tpr -b $BTIME -o ../tmp$k.xtc
cd $RUNROOT/ConstantT/a_${alpha}/$Temp/$i
if [ -e "tmp$k.xtc" ]; then
	echo "tmp$k.xtc" >> file
fi
#done

#file=\`awk '{printf "%s ", \$1} END {printf "\n"}' file\`
#trjcat_4.5.7 -f \$file -cat -o trajout.xtc
#rm tmp?.xtc
EOF
				# qsub -hold_jid a_${alpha_old}_f_PC_BigMDSteps ConstantT.sh
				jobid1=$(sbatch --dependency=afterok:$jobidf --parsable ConstantT.sh)
				echo a_${alpha}_${Temp}_${i}_PC_BigMDSteps_1: $jobid1
				jobids="$jobids,$jobid1"
			done
			jobid=${jobids#*,}
#cat > ConstantT.sh <<EOF
#!/bin/bash
#SBATCH --job-name a_${alpha}_${Temp}_${i}_PC_BigMDSteps_2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node $NCPU
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $RUNROOT/jobout/%j.out
#SBATCH --error $RUNROOT/jobout/%j.err

#module load mpi/openmpi-4.1.5
#export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:\$PATH
#export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:\$MODULEPATH
#module load pl25mod
#module load matlab/2023b

#for k in \`seq 5 1 9\`
#do
#	mkdir -p $RUNROOT/ConstantT/a_${alpha}/$Temp/$i/\$k
#	cd $RUNROOT/ConstantT/a_${alpha}/$Temp/$i/\$k
#	l=\$(($j+\$k))
#	PDBfile=$PDBDIR/N\$l.gro	
##EM
#	$GROMACSDIR/editconf_mpi -f \$PDBfile -o N\${l}_center.gro -box $box $box $box -bt cubic -c
#	$GROMACSDIR/grompp_mpi -f $RUNROOT/Polymer_EM.mdp -c N\${l}_center.gro -p $TOPfile -o em.tpr > output 2>&1 
#	mpirun -np $NCPU $GROMACSDIR/mdrun_mpi -noddcheck -deffnm em -tableb $RUNROOT/Tables/table.xvg -table $Tablefile -tablep $Tablepfile -plumed $RUNROOT/plumed_SR.dat -pd
##MD
#	$GROMACSDIR/grompp_mpi -f ../../$Temp.mdp -c em.gro -p $TOPfile -o md.tpr > output 2>&1 
#	mpirun -np $NCPU $GROMACSDIR/mdrun_mpi -noddcheck -s md.tpr -tableb $RUNROOT/Tables/table.xvg -table $Tablefile -tablep $Tablepfile -plumed $RUNROOT/plumed_SR.dat -pd
#
#	echo "0" | $GROMACSDIR/trjconv_mpi -f traj.xtc -s md.tpr -b $BTIME -o ../tmp\$k.xtc
#	cd $RUNROOT/ConstantT/a_${alpha}/$Temp/$i
#	if [ -e "tmp\$k.xtc" ]; then
#		echo "tmp\$k.xtc" >> file
#	fi
#done
#
#file=\`awk '{printf "%s ", \$1} END {printf "\n"}' file\`
#trjcat_4.5.7 -f \$file -cat -o trajout.xtc
#rm tmp?.xtc
#EOF
#			# qsub -hold_jid a_${alpha_old}_f_PC_BigMDSteps ConstantT.sh
#			jobid2=$(sbatch --dependency=afterok:$jobidf --parsable ConstantT.sh)
#			echo a_${alpha}_${Temp}_${i}_PC_BigMDSteps_2: $jobid2
cat > Analysis.sh << EOF
#!/bin/bash
#SBATCH --job-name b_${alpha}_${Temp}_${i}_PC_BigMDSteps
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --exclude cpu1-[1,75,76]
#SBATCH --mem 8000
#SBATCH --partition i64m512u
#SBATCH --output $RUNROOT/jobout/%j.out
#SBATCH --error $RUNROOT/jobout/%j.err

export PATH=/hpc2hdd/home/chu-amat/opt/gromacs-4.5.7/bin:\$PATH
export MODULEPATH=/hpc2hdd/home/chu-amat/opt/plumed-2.5.0/lib/plumed/:\$MODULEPATH
module load pl25mod
module load matlab/2023b

file=\`awk '{printf "%s ", \$1} END {printf "\n"}' file\`
trjcat_mpi -f \$file -cat -o trajout.xtc
rm tmp?.xtc
echo "0" | trjconv_mpi -f trajout.xtc -s $RUNROOT/TOP_ME/Polymer_${Nres}_0.4.pdb -o tmp.pdb
rm -r ?
$RUNROOT/fcal.x tmp.pdb $RUNROOT/Contact_a/Probability.dat $Nres
rm tmp.pdb
#rm \#*
EOF
			# qsub -hold_jid a_${alpha}_${Temp}_${i}_PC_BigMDSteps_1,a_${alpha}_${Temp}_${i}_PC_BigMDSteps_2 Analysis.sh
			jobid=$(sbatch --dependency=afterok:$jobid --parsable Analysis.sh)
			echo b_${alpha}_${Temp}_${i}_PC_BigMDSteps: $jobid
			echo $jobid >> $RUNROOT/ConstantT/a_$alpha/$Temp/list.dat
		done
		cd $RUNROOT
		sed "s/ALPHA/$alpha/g" f.sh | sed "s/Lambda/$Lambda/g" | sed "s/MAXNUM/$MAXNUM/g" > f_$alpha.sh
		# hold_jid=`awk '{printf "%s,", $1} END {printf "%s ",$1}' $RUNROOT/ConstantT/a_$alpha/$Temp/list.dat`
		hold_jid=`awk '{if(NR>1) printf ","; printf "%s", $1}' $RUNROOT/ConstantT/a_$alpha/$Temp/list.dat`
		echo $hold_jid
		echo "sbatch --dependency=afterok:$hold_jid f_$alpha.sh" > f_$alpha.submit
		# echo "qsub -hold_jid $hold_jid f_$alpha.sh" > f_$alpha.submit
		# qsub -hold_jid $hold_jid f_$alpha.sh
		jobidf=$(sbatch --dependency=afterok:$hold_jid --parsable f_$alpha.sh)
		echo a_${alpha_old}_f_PC_BigMDSteps: $jobidf
	done
done
