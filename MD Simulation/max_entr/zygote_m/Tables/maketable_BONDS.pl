#!/usr/bin/perl

###############################################################################
# maketable4.pl will make a table of 1/R^10, 1/R^12 and the appropriate       #                           
# negative first derivatives for use with the CA model in gromacs: VERSION 4  #
# Values not provided for distances under 0.1 Angstroms                       #
# WARNING: Always look at the table to ensure no precision issues are present #
# Written by Paul Whitford, 12/6/08                                           #
###############################################################################


#what is the length of the table? (in nm. 100 in this example)
        $Rtable=10;
	$sigma=0.4;
#	$cutoff1=$sigma*2**(1.0/6.0);
	$cutoff1=$sigma;
	$cutoff2=$sigma*1.5;
# what is the spacing of the table? (0.005 nm, here)
# NOTE: Gromacs internally uses a spacing of 0.002, so you may want to use 0.002 for consistency.
	$DR=0.002;
        $Ntable=int($Rtable/$DR);
#	print "0.0 0.0 0.0\n";
	$R=0.002;
#	$R1 = 4*((0.5*$sigma/$R)**2-(0.5*$sigma/$R)**1)+1-0.5*30*1.5**2*log(1-($R/$cutoff2)**2);
#	$R2 = 4*(2*(0.5*$sigma)**2/$R**3-1*(0.5*$sigma)**1/$R**2)-30/($sigma)**2/(1-($R/$cutoff2)**2)*$R;
	$R1 = 4*(($sigma/$R)**12-($sigma/$R)**6)+1-0.5*30*1.5**2*log(1-($R/$cutoff2)**2);
	$R2 = 4*(12*$sigma**12/$R**13-6*$sigma**6/$R**7)-30/($sigma)**2/(1-($R/$cutoff2)**2)*$R;
	print "0.0 $R1 $R2\n";
for($i=1; $i<$Ntable; $i++){
       $R=$i*$DR;

	if($R >= 0.002 && $R < $cutoff1){
		$R1 = 4*(($sigma/$R)**12-($sigma/$R)**6)+1-0.5*30*1.5**2*log(1-($R/$cutoff2)**2);
		$R2 = 4*(12*$sigma**12/$R**13-6*$sigma**6/$R**7)-30/($sigma)**2/(1-($R/$cutoff2)**2)*$R;
#		$R1 = 4*((0.5*$sigma/$R)**2-(0.5*$sigma/$R)**1)+1-0.5*30*1.5**2*log(1-($R/$cutoff2)**2);
#		$R2 = 4*(2*(0.5*$sigma)**2/$R**3-1*(0.5*$sigma)**1/$R**2)-30/($sigma)**2/(1-($R/$cutoff2)**2)*$R;
		print "$R $R1 $R2\n";
	}elsif($R >= $cutoff1 && $R < $cutoff2-0.001-0.05){
		$R1 = -0.5*30*1.5**2*log(1-($R/$cutoff2)**2);
		$R2 = -30/($sigma)**2/(1-($R/$cutoff2)**2)*$R; # Unit
		print "$R $R1 $R2\n";
	}elsif($R >= $cutoff2-0.001-0.05){
		$R1 = 100.0*($R/$cutoff2)**5;
		$R2 = -100.0*(1.0/$cutoff2)**5*5*($R)**4;
		print "$R $R1 $R2\n";
	}else{
#		print "$R 0.0 0.0\n";
	}
}
