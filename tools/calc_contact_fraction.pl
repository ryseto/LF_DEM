#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;

$interaction_data = $ARGV[0];
$shearrate_data = $ARGV[1];
# Create output file name
$i = index($interaction_data, 'int_', 0)+4;
$j = index($interaction_data, '.dat', $i-1);
$name = substr($interaction_data, $i, $j-$i);

$i = index($name, 'N', 0)+1;
$j = index($name, 'VF', $i-1);
$number = substr($name, $i, $j-$i) ;
printf "num = $number \n";
printf "$name\n";

$output = "confrac_$name.dat";
open (IN_rheo, "< rheo_${name}.dat");

open (IN_stress, "< eigen_${name}.dat");

open (OUT, "> ${output}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
while (1){
	&InInteractions;
	last unless defined $line;
	$num ++;
}
close (OUT);

sub readHeader{
	$line = <IN_interaction>;
	$line = <IN_interaction>;
	$line = <IN_interaction>;
	$line = <IN_interaction>;
	$line = <IN_interaction>;
	$line = <IN_interaction>;
	#	printf "$line";
	$cnt = 0;
	while (1){
		$line = <IN_rheo>;
		($buf, $buf2) = split(/\s+/, $line);
		#printf "$buf\n";
		if ($cnt >= 52){
			return;
		}
		$cnt ++;
	}

}

sub readRheo{
	while (1) {
		$line = <IN_rheo>;
		($d1, $d2, $d3, $d4, $d5, $d6, $d7, $d8, $d9, $d10,
		$d11, $d12, $d13, $d14, $d15, $d16, $d17, $d18, $d19, $d20,
		$d21, $d22, $d23, $d24, $d25, $d26, $d27, $d28, $d29, $d30,
		$d31, $d32, $d33, $d34, $d35, $d36, $d37, $d38, $d39, $d40,
		$d41, $d42, $d43, $d44, $d45, $d46, $d47) = split(/\s+/, $line);
		#	last unless defined $line;
		$shear_strain_data = $d1;
		$viscosity_data = $d2;
		$pressure = $d38;

		
		if ($shear_strain_data >= $shear_strain){
			#			printf "$shearrate_data :::: $viscosity_data\n";
			#rintf "a = $shear_strain_data b =$viscosity_data\n";
			return;
		}
	}
	
}

sub readStress{
	
	while(1){
		$line = <IN_stress>;
		($shear_strain_data, $eigen1) = split(/\s+/, $line);
		printf "$shear_strain_data --- $shear_strain \n";
		if ($shear_strain_data >= $shear_strain){
			$stress_eigen_value = - $eigen1;
			#printf "$shear_strain_data :::: $shear_strain\n";
			#rintf "a = $shear_strain_data b =$viscosity_data\n";
			return;
		}
	}
	
	
}

sub InInteractions {
	$line = <IN_interaction>;
	($buf, $shear_strain, $num_interaction) = split(/\s+/, $line);
	#printf "$line\n";
	if ($buf != "#"){
		exit;
	}
	#	printf "shear strain ---> $shear_strain \n";
	
	
	&readRheo;
	&readStress;

	# 1, 2: numbers of the interacting particles
	# 3: 1=contact, 0=apart
	# 4, 5, 6: normal vector
	# 7: dimensionless gap = s - 2, s = 2r/(a1+a2)
	# 8: lubrication force
	# 9: Normal part of contact force
	# 10: Tangential part of contact force
	# 11: Colloidal force
	# 12: Viscosity contribution of contact xF
	# 13: N1 contribution of contact xF
	# 14: N2 contribution of contact xF
	
	$count_contact = 0;
	$count_frictional = 0;
	
	#	printf "$num_interaction\n";

	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction> ;
		#		/* 1, 2: numbers of the interacting particles
		#		* 3: 1=contact, 0=apart
		#		* 4, 5, 6: normal vector
		#		* 7: dimensionless gap = s - 2, s = 2r/(a1+a2)
		#		* 8: normal     of lubrication force
		#		* 9: tangential of lubrication force
		#		* 10: normal part     of contact force
		#		* 11: tangential part of contact force
		#		* 12: normal colloidal force
		#		* 13: Viscosity contribution of contact xF
		#		* 14: N1 contribution of contact xF
		#		* 15: N2 contribution of contact xF
		#		* 16: friction state
		#		*      0 = not frictional
		#		*      1 = non-sliding
		#		*      2 = sliding
		#		*/
		($i, $j, $contact, $nx, $ny, $nz, #1---6
		$gap, $f_lub_norm, $f_lub_tan , $fc_n, $fc_tan, $fcol, #7--12
		$sxz_cont_xF, $n1_cont_xF, $n2_cont_xF, $friction, # 13--16
		) = split(/\s+/, $line);
		
		
		### critical load model
		if ($contact == 1){
			$count_contact ++;
		}
		if ($friction != 0){
			$count_frictional ++;
		}
		### repulsion
		#		printf "$gap\n";

#		if ($gap < 0.02){
#			$count_contact ++;
#
#		}
#		if ($contact == 1){
#			$count_frictional ++;
#		}
		
	}

	if ( $count_contact != 0) {
		if ($shear_strain_data > 5) {
			$contact_fraction = $count_frictional / $count_contact ;
			$stress = $viscosity_data*$shearrate_data;
			$stress_eigen_value_nonscale = $stress_eigen_value*$shearrate_data;
			$pressure_nonscale = $pressure*$shearrate_data;
			$contact_number = 2*$count_contact / $number;
			$fric_contact_number = 2*$count_frictional  / $number;
			
			printf OUT "$contact_fraction $stress $stress_eigen_value_nonscale $pressure_nonscale $contact_number $fric_contact_number\n";
		}

	}
}

