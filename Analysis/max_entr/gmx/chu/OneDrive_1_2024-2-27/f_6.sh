#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N a_6_f_PC_BigMDSteps
#$ -j y
#$ -o f_6.out
#$ -q cpu_short
#$ -P otherprj

RUNROOT=`pwd`
Temp=120.28
Nres=577
sigma=0.4
alpha=6
protein=chr14
	
ContactPrefile="$RUNROOT/Contact_a/Contact_a_${alpha}.dat"
cd $RUNROOT/ConstantT/a_$alpha/$Temp
WORKDIR="$RUNROOT/ConstantT/a_$alpha/$Temp/Probability"
mkdir -p $WORKDIR
rm $WORKDIR/Filelist*

for i in `seq 0 1 9`
do
	DATADIR="$RUNROOT/ConstantT/a_$alpha/$Temp/$i"
	cd $DATADIR
	if [ -e "$DATADIR/f.dat" ];then
	echo "$DATADIR/f.dat" >> $WORKDIR/Filelist_f
	fi
	if [ -e "$DATADIR/B.bin" ];then
	echo "$DATADIR/B.bin" >> $WORKDIR/Filelist_B
	fi
	if [ -e "$DATADIR/Probability.dat" ];then
	echo "$DATADIR/Probability.dat" >> $WORKDIR/Filelist_P
	fi
done

cd $WORKDIR
cp $RUNROOT/Contact_a/Probability.dat Probability_ND.dat
cp $RUNROOT/fsum1.m .
cp $RUNROOT/fsum2_1.m .
cp $RUNROOT/fsum2_2.m .
cp $RUNROOT/fave.x .
octave fsum1.m
./fave.x f.dat
sync
octave-cli fsum2_1.m
sync
octave-cli fsum2_2.m

cp $RUNROOT/Contact_a/Contact_a_${alpha}.dat .
sed "s/TEMP/$Temp/g" $RUNROOT/fsum3.m | sed "s/ContactalphaFILE/Contact_a_${alpha}.dat/g" | sed "s/LAMBDA/0.20/g"> fsum3.m
sync
octave-cli fsum3.m

for i in `seq 0 1 9`
do
	DATADIR="$RUNROOT/ConstantT/a_$alpha/$Temp/$i"
	cd $DATADIR
	rm B.bin f.dat Probability.dat tmp.pdb \#*
done

cd $WORKDIR
alphaNext=$(($alpha+1))
cp Contact_a.dat $RUNROOT/Contact_a/Contact_a_${alphaNext}.dat
rm Contact_a_${alpha}.dat B*bin f.dat fdata.bin

mkdir -p $RUNROOT/ConstantT/a_$alphaNext
cd $RUNROOT/ConstantT/a_$alphaNext
$RUNROOT/TOP_ME/TOP.x $RUNROOT/TOP_ME/Polymer_${Nres}_${sigma}.pdb $RUNROOT/Contact_a/Contact_a_${alphaNext}.dat ${protein}_a_${alphaNext}

cd $WORKDIR	
awk '$2>=$1+2 {print $2,$1,$3}' $RUNROOT/uij/PC_100000_iced_chr14_dense.matrix_P.dat > PU.tmp.dat
awk '$2>=$1+2 {print $1,$2,$3}' Probability.dat    > PD.tmp.dat
seq 1 1 $(($Nres-1)) | awk '{print $1,$1,"\n",$1,$1+1,"\n",$1+1,$1}' > tmp.dat
echo "$Nres $Nres" >> tmp.dat
cat tmp.dat | awk '{print $1,$2,1}' | cat - PU.tmp.dat PD.tmp.dat | sort -k1n -k2n > Probability.Mixed.dat
rm P?.tmp.dat tmp.dat
	
ratio=`cat tol.dat | awk '{printf "%.2f\n",$1*100}'`
sed "s/Alpha/$alpha/g" $RUNROOT/ContactPlot.m | sed "s/RATIO/$ratio/g" > ContactPlot.m
#matlab < ContactPlot.m
#cp Contact_a_${alpha}.pdf $RUNROOT/Contact_a

$RUNROOT/uij/CP_GD.x Probability.dat $Nres
$RUNROOT/uij/TAD.x Probability.dat $Nres

cp -r $RUNROOT/Analysis .
cd Analysis
sh RUN.sh
	
echo "#alpha tolerance" > $RUNROOT/Contact_a/tol.dat
for a in `seq 0 1 $alpha`
do
	tolfile="$RUNROOT/ConstantT/a_$a/$Temp/Probability/tol.dat"
	tol=`awk '{print $1}' $tolfile`
	echo "$a $tol" >> $RUNROOT/Contact_a/tol.dat
done
