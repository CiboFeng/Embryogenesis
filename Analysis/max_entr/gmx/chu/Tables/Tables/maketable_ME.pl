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
        $Rtable=10;
# what is the spacing of the table? (0.005 nm, here)
# NOTE: Gromacs internally uses a spacing of 0.002, so you may want to use 0.002 for consistency.
	$DR=0.002;
#	$rc=1.76;
#	$miu=3.72*1.0/0.4;
#	$rc=1.50;
#	$miu=2.67*1.0/0.4;
#	$rc=1.76;
#	$miu=3.72*1.0/0.4;
	$rc=1.50;
	$miu=3.00;
        $Ntable=int($Rtable/$DR);
	print "0.0 0.0 0.0 0.0 0.0 1.0 0.0\n";

for($i=1; $i<=$Ntable; $i++){
       $R=$i*$DR;

	if($R > 0.00){
	       $R1=0;
	       $R2=0;
	       $R3=0.5*(1+tanh($miu*($rc*0.4-$R))); # Pay attention to the unit !!! It's nm here !!!
	       $R4=0.5*(1-tanh($miu*($rc*0.4-$R))**2)*$miu; # Pay attention to the unit !!! It's nm here !!!
	print "$R 0.0 0.0 $R1 $R2 $R3 $R4\n";
	}
}
