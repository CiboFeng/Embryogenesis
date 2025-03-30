#!/usr/bin/perl

###############################################################################
# maketable4.pl will make a table of 1/R^10, 1/R^12 and the appropriate       #                           
# negative first derivatives for use with the CA model in gromacs: VERSION 4  #
# Values not provided for distances under 0.1 Angstroms                       #
# WARNING: Always look at the table to ensure no precision issues are present #
# Written by Paul Whitford, 12/6/08                                           #
###############################################################################


#what is the length of the table? (in nm. 100 in this example)
        $Rtable=180;
# what is the spacing of the table? (0.005 nm, here)
# NOTE: Gromacs internally uses a spacing of 0.002, so you may want to use 0.002 for consistency.
	$DR=0.05;
        $Ntable=int($Rtable/$DR);
	$pi=3.1415926;
	print "0.0 2.0 0.0\n";
for($i=1; $i<=$Ntable; $i++){
       $R=$i*$DR;
		$R1 = 1-cos(($R-180)/180*$pi);
#		$R2 = -sin(($R-180)/180*$pi);
		$R2 = -$pi/180*sin(($R-180)/180*$pi);
		print "$R $R1 $R2\n";
}
