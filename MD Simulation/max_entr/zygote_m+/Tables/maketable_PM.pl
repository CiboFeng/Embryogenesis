#!/usr/bin/perl

###############################################################################
# maketable4.pl will make a table of 1/R^10, 1/R^12 and the appropriate       #                           
# negative first derivatives for use with the CA model in gromacs: VERSION 4  #
# Values not provided for distances under 0.1 Angstroms                       #
# WARNING: Always look at the table to ensure no precision issues are present #
# Written by Paul Whitford, 12/6/08                                           #
###############################################################################

  	use Math::Trig;

#what is the length of the table? (in nm. 100 in this example)
        $Rtable=500;
	$sigma=0.4;
#	$sigma_2=$sigma*0.1*(2.0);
	$sigma_r=0.4454493590701696;
#	$sigma_r=0.5345392308842036;
	$sigma_2=0.2*$sigma_r*2**(1.0/6.0);
	$cutoff1=$sigma*1.0/((1+sqrt(2.0))/2)**(1.0/6.0);
	$cutoff2=$sigma*2**(1.0/6.0);
# what is the spacing of the table? (0.005 nm, here)
# NOTE: Gromacs internally uses a spacing of 0.002, so you may want to use 0.002 for consistency.
	$DR=0.002;
        $Ntable=int($Rtable/$DR);
	$R=0.002;
#	$R3 = 4*(($sigma*0.2/$R)**12-($sigma*0.2/$R)**6)+1+4;
#	$R4 = 4*(12*($sigma*0.2)**12/$R**13-6*($sigma*0.2)**6/$R**7);
	$R3 = 4.0;
	$R4 = 0.0;
	print "0.0 0.0 0.0 0.0 0.0 $R3 $R4\n";

for($i=1; $i<$Ntable; $i++){
        $R=$i*$DR;
	
	if($R >= 0.001 && $R < $sigma_2){
#		$R3 = 4*(($sigma*0.1/$R)**2-($sigma*0.1/$R)**1)+1+4;
#		$R4 = 4*(2*($sigma*0.1)**2/$R**3-1*($sigma*0.1)**1/$R**2);
		$R3 = 4*(($sigma_r*0.2/$R)**12-($sigma_r*0.2/$R)**6)+1+4;
		$R4 = 4*(12*($sigma_r*0.2)**12/$R**13-6*($sigma_r*0.2)**6/$R**7);
		print "$R 0.0 0.0 0.0 0.0 $R3 $R4\n";
#		print "$R 0.0 0.0 0.0 0.0 4.0 0.0\n";
	}elsif($R >= $sigma_2 && $R < $cutoff1){
		$R3 = 2*(1+tanh(0.5*(4*(($sigma/$R)**12-($sigma/$R)**6)+1)-1));
		$R4 = 2*(1-tanh(0.5*(4*(($sigma/$R)**12-($sigma/$R)**6)+1)-1)**2)*0.5*(4*(12*$sigma**12/$R**13-6*$sigma**6/$R**7));
		print "$R 0.0 0.0 0.0 0.0 $R3 $R4\n";
	}elsif($R >= $cutoff1 && $R < $cutoff2){
		$R3 = 4*(($sigma/$R)**12-($sigma/$R)**6)+1;
		$R4 = 4*(12*$sigma**12/$R**13-6*$sigma**6/$R**7);
		print "$R 0.0 0.0 0.0 0.0 $R3 $R4\n";
	}elsif($R >= $cutoff2){
		print "$R 0.0 0.0 0.0 0.0 0.0 0.0\n";
		# notes: we put these in so that the spacing is equal throughout the table (gromacs returns an error otherwise).  But, atoms should not get this close.
	}
}
