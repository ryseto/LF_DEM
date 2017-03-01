#!/usr/bin/perl

# Usage:
# $ convertMaterialFunctions.pl st_[...].dat

use Math::Trig;
use IO::Handle;
use Getopt::Long;

my $st_data = $ARGV[0];

# Create output file name
$i = index($st_data, 'st_', 0)+3;
$j = index($st_data, '.dat', $i-1);
$name = substr($st_data, $i, $j-$i);

$ii = index($name, '_extension', 0);
if ($ii == -1) {
	$flow_type = "shear";
} else {
	$flow_type = "extension";
}

if ($flow_type eq "shear") {
	$rate = 0.5;
	@dtensor = (0, $rate, 0,
				$rate, 0, 0,
				0, 0, 0);
	@e0tensor = (-0.5*$rate, 0, 0,
				 0, -0.5*$rate, 0,
				 0,    0, 1*$rate);
	@e1tensor = (-$rate, 0, 0,
	0, $rate, 0,
	0, 0, 0);
	
	@e2tensor = (-$rate, 0, 0,
	0, $rate, 0,
	0, 0, 0);
	
	@e3tensor = (-$rate, 0, 0,
				  0, $rate, 0,
				  0, 0, 0);
} elsif ($flow_type eq "extension") {
	$rate = 1;
	@dtensor = ($rate,  0, 0,
				0, -$rate, 0,
				0,  0, 0);
	@e0tensor = (-0.5*$rate,    0, 0,
				    0, -0.5*$rate, 0,
				    0,   0 , $rate);
	@e3tensor = (0, $rate, 0,
				 $rate, 0, 0,
				 0, 0, 0);
}

$dd = 0;
for ($i=0; $i<9; $i++) {
	$dd += $dtensor[$i]*$dtensor[$i];
}
$e0e0 = 0;
for ($i=0; $i<9; $i++) {
	$e0e0 += $e0tensor[$i]*$e0tensor[$i];
}
$e3e3 = 0;
for ($i=0; $i<9; $i++) {
	$e3e3 += $e3tensor[$i]*$e3tensor[$i];
}

printf "$dd $e0e0 $e3e3\n";


$output = "mf_$name.dat";
printf "flow type is $flow_type\n";

open (OUT, "> ${output}");
open (IN_stress, "< ${st_data}");

$kappa_sum = 0;
$lambda0_sum = 0;
$lambda3_sum = 0;
$cnt_sum = 0;

&readHeader;
&InStress;

sub InStress {
	while (1) {
		$line = <IN_stress>;
		last unless defined $line;
		($time, $strain, $shear_rate, $sxx, $sxy, $sxz, $syz, $syy, $szz) = split(/\s+/, $line);
		@stresstensor = ($sxx, $sxz, $sxy,
		                 $sxz, $szz, $syz,
						 $sxy, $syz, $syy);
		#printf "$strain\n";
		$kappa = 0;
		for ($i=0; $i<9; $i++) {
			$kappa += $stresstensor[$i]*$dtensor[$i]/$dd ;
		}
		$lambda0 = 0;
		for ($i=0; $i<9; $i++) {
			$lambda0 += $stresstensor[$i]*$e0tensor[$i]/$e0e0;
		}
		$lambda3 = 0;
		for ($i=0; $i<9; $i++) {
			$lambda3 += $stresstensor[$i]*$e3tensor[$i]/$e3e3;
		}
		if ($strain > 1) {
			$kappa_sum += $kappa;
			$lambda0_sum += $lambda0;
			$lambda3_sum += $lambda3;
			$cnt_average += 1;
		}
		
		printf OUT "$strain $kappa $lambda0 $lambda3\n";
	}
	$kappa_average   = $kappa_sum/$cnt_average;
	$lambda0_average = $lambda0_sum/$cnt_average;
	$lambda3_average = $lambda3_sum/$cnt_average;
	
	printf  "$kappa_average $lambda0_average $lambda3_average\n";
}

sub readHeader {
	$line = <IN_stress>;
	$line = <IN_stress>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_stress>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_stress>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_stress>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_stress>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	$line = <IN_stress>; ($buf, $buf, $flow_type) = split(/\s+/, $line);
	if ($Ly == 0) {
		$number_of_header = 10;
	} else {
		$number_of_header = 10;
	}
	for ($i = 0; $i<$number_of_header; $i++) {
		$line = <IN_stress>;
		#printf "$line";
	}
	#	printf "=====\n";
}
