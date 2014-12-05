#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;

$force_factor = 0.000001;

$y_section = 0;
$yap_radius = 1;

$particle_data = $ARGV[0];
if ($#ARGV >= 1){
	$force_factor = $ARGV[1];
}
if ($#ARGV == 2){
	$y_section = $ARGV[2];
	printf "section $y_section\n";
}

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);


$interaction_data = "int_${name}.dat";
printf "$interaction_data\n";
$output = "y_$name.yap";
$output2 = "nvec_$name.dat";
$out_gaps = "gaps_$name.dat";
$out_mathm_pos = "mpos_$name.dat";
$out_mathm_int = "mint_$name.dat";

open (IN_rheo, "< rheo_${name}.dat");

${sum_fmax} = 0;
${cnt} = 0;
while (1) {
	$line = <IN_rheo>;
	($d1, $d2, $d3, $d4, $d5, $d6, $d7, $d8, $d9, $d10,
	$d11, $d12, $d13, $d14, $d15, $d16, $d17, $d18, $d19, $d20,
	$d21, $d22, $d23, $d24, $d25, $d26, $d27, $d28, $d29, $d30,
	$d31, $d32, $d33, $d34, $d35, $d36, $d37, $d38, $d39, $d40) = split(/\s+/, $line);
	last unless defined $line;
	if ($d1 > 2){
		${sum_fmax} += $d29;
		${cnt} ++;
	}
}

#$fmax_ave = ${sum_fmax}/${cnt};
#printf "fmax = $fmax_ave $cnt \n";
#$force_factor = 0.3/$fmax_ave;
$force_factor = 0.0003;
#$force_factor = 0.003;
#printf  "$fmax_ave\n";
#exit;

open (OUT, "> ${output}");
open (OUT2, "> ${output2}");
open (OUTG, "> ${out_gaps}");
open (OUTMP, "> ${out_mathm_pos}");
open (OUTMI, "> ${out_mathm_int}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
$first=1;
$c_traj=0;
$num = 0;

while (1){
	
	&InParticles;
	
	&InInteractions;
	
	printf "$shear_rate $shear_stress\n";
	last unless defined $line;
	&OutYaplotData;
	$num ++;

}
close (OUT);

sub readHeader{
	
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	for ($i = 0; $i < 16; $i++) {
		$line = <IN_particle>;
	}

	
	for ($i = 0; $i < 24; $i++) {
		$line = <IN_interaction>;
	}
	printf "$np, $VF, $Lx, $Ly, $Lz\n";
}

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
    ($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
	
	# h_xzstress << sp << c_xzstressXF << sp << c_xzstressGU << sp << b_xzstress
	# 1: number of the particle
	# 2: radius
	# 3, 4, 5: position
	# 6, 7, 8: velocity
	# 9, 10, 11: angular velocity
	# 12: viscosity contribution of lubrication
	# 13: viscosity contributon of contact GU xz
	# 14: viscosity contributon of brownian xz
	# (15: angle for 2D simulation)
	
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_particle> ;
        ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
		$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
		
		if (true){
			#printf OUTMP "$line";
			printf OUTMP "$i $x $y $z $a\n";
		}
		
		$radius[$i] = $a;
        $posx[$i] = $x;
        $posy[$i] = $y;
        $posz[$i] = $z;
		$omegay[$i] = $oy;
		$ang[$i] = $angle;
		if ($radius_max < $a){
			$radius_max = $a;
		}
    }
	
	#	if ($c_traj == 0) {
	#		$min_dist_origin = 100;
	#		$min_dist_origin1 = 100;
	#		$min_dist_origin2 = 100;
	#		for ($i = 0; $i < $np; $i ++) {
	#			#$sqdistance = ($posx[$i]*$posx[$i]+$posy[$i]*$posy[$i]+$posz[$i]*$posz[$i]);
	#			$sqdistance = &calcsqdist($posx[$i],$posy[$i],$posz[$i], 0,0,0);
	#			if ($sqdistance < $min_dist_origin) {
	#				$min_dist_origin = $sqdistance;
	#				$center = $i;
	#			}
	#			$sqdistance = &calcsqdist($posx[$i],$posy[$i],$posz[$i], 0,0,5);
	#			if ($sqdistance < $min_dist_origin1) {
	#				$min_dist_origin1 = $sqdistance;
	#				$uppder = $i;
	#			}
	#			$sqdistance = &calcsqdist($posx[$i],$posy[$i],$posz[$i], 0,0,-5);
	#			if ($sqdistance < $min_dist_origin2) {
	#				$min_dist_origin2 = $sqdistance;
	#				$lower = $i;
	#			}
	#
	#		}
	#	}
	#
	#	$trajx[$c_traj] = $posx[$center];
	#	$trajy[$c_traj] = $posy[$center];
	#	$trajz[$c_traj] = $posz[$center];
	#	$trajx2[$c_traj] = $posx[$uppder];
	#	$trajy2[$c_traj] = $posy[$uppder];
	#	$trajz2[$c_traj] = $posz[$uppder];
	#	$trajx3[$c_traj] = $posx[$lower];
	#	$trajy3[$c_traj] = $posy[$lower];
	#	$trajz3[$c_traj] = $posz[$lower];
	$c_traj++;
}

sub calcsqdist {
    ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
	$dist = ($x1-$x2)*($x1-$x2);
	$dist += ($y1-$y2)*($y1-$y2);
	$dist += ($z1-$z2)*($z1-$z2);
	return $dist;
}

sub InInteractions {
	$line = <IN_interaction>;
	($buf, $shear_strain_i, $num_interaction) = split(/\s+/, $line);

#	if ($shear_strain_i != $shear_strain) {
#		printf "$shear_strain_i  !=  $shear_strain\n";
#		exit(1);
#	}

	if ($buf != "#"){
		exit(1);
	}
	printf OUTG "$shear_rate\n";
	
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
	
	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction> ;
		#		"#1: particle 1 label\n"
		#		"#2: particle 2 label\n"
		#		"#3: contact state (0 = no contact, 1 = frictionless contact, 1 = non-sliding frictional, 2 = sliding frictional)\n"
		#		"#4: normal vector, oriented from particle 1 to particle 2 x\n"
		#		"#5: normal vector, oriented from particle 1 to particle 2 y\n"
		#		"#6: normal vector, oriented from particle 1 to particle 2 z\n"
		#		"#7: dimensionless gap = s-2, s = 2r/(a1+a2)\n"
		#		"#8: norm of the normal part of the lubrication force\n"
		#		"#9: tangential part of the lubrication force x\n"
		#		"#10: tangential part of the lubrication force y\n"
		#		"#11: tangential part of the lubrication force z\n"
		#		"#12: norm of the normal part of the contact force\n"
		#		"#13: tangential part of the contact force, x\n"
		#		"#14: tangential part of the contact force, y\n"
		#		"#15: tangential part of the contact force, z\n"
		#		"#16: norm of the normal repulsive force\n"
		#		"#17: Viscosity contribution of contact xF\n";
		($i, $j, $contact, $nx, $ny, $nz, #1---6
		$gap, $f_lub_norm, # 7, 8
		$f_lub_tan_x, $f_lub_tan_y, $f_lub_tan_z, # 9, 10, 11
		$fc_norm, # 12
		$fc_tan_x, $fc_tan_y, $fc_tan_z, # 13, 14, 15
		$fr_norm, $s_xF
		) = split(/\s+/, $line);
		if ($num==$num_mathm){
			printf OUTMI "$line";
		}
		# $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
		$int0[$k] = $i;
		$int1[$k] = $j;
		$contactstate[$k] = $contact;
		$domega[$k] = $omegay[$i] - $omegay[$j];
		$F_lub[$k] = $f_lub_norm;
		$Fc_n[$k] = $fc_norm;
		$Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
		#$f_normal = $fc_norm + $fr_norm + $f_lub_norm;
		#$f_normal = $f_lub_norm;
		$f_normal = $fc_norm;
		$force[$k] = $f_normal;

		$S_bf[$k] =  $s_xF;
		$fricstate[$k] = $friction;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		$ft_x[$k] = $f_lub_tan_x + $fc_tan_x;
		$ft_y[$k] = $f_lub_tan_y + $fc_tan_y;
		$ft_z[$k] = $f_lub_tan_z + $fc_tan_z;
		
		$Gap[$k] = $gap;
		printf OUTG "$gap ";
	}
	printf OUTG "\n";
}

sub OutYaplotData{
	
	if ($first == 0){
		printf OUT "\n";
	} else {
		$first = 0;
	}
	$postext = $Lz/2+2;
	$postext2 = $Lz/2+3;
	$shear_rate_text = int($shear_rate*1e5);
	$shear_rate_text *= 1e-5;

	
	#printf OUT "y 10\n";
	#printf OUT "@ 2\n";
	#printf OUT "t -2 0 $postext shear rate = $shear_rate_text \n";
	#printf OUT "t -2 0 $postext2 shear stress = $shear_stress \n";
	
	#$meterx = $Lx/2+1;
	#$meterbottom = -$Lz/2;
	#$metertop = $meterbottom + $shear_rate*1000;
	printf OUT "@ 3\n";
	&OutMeter($shear_rate, 0.01, $Lx/2+2, -$Lz/2+2);
	printf OUT "@ 4\n";
	#&OutMeter($shear_stress, 2, -$Lx/2-3.5, -$Lz/2+2);
	&OutStress($shear_stress, 5);


	
	#
	#	printf OUT "y 7\n";
	#	printf OUT "r 0.1\n";
	#    printf OUT "@ 5\n";
	#	for ($i = 0; $i < $c_traj; $i++){
	#		$xs = $trajx[$i];
	#		$ys = $trajy[$i];
	#		$zs = $trajz[$i];
	#		$xe = $trajx[$i+1];
	#		$ye = $trajy[$i+1];
	#		$ze = $trajz[$i+1];
	#		if (abs($zs-$ze) < 1
	#			&& abs($xs-$xe) < 1
	#			) {
	#			printf OUT "l $xs $ys $zs $xe $ye $ze\n";
	#		}
	#		#		$xs = $trajx[$i];
	#		#$ys = $trajy[$i];
	#		#$zs = $trajz[$i];
	#		#		printf OUT "c $xs $ys $zs\n";
	#	}
	#    printf OUT "@ 3\n";
	#	for ($i = 0; $i < $c_traj; $i++){
	#		$xs = $trajx2[$i];
	#		$ys = $trajy2[$i];
	#		$zs = $trajz2[$i];
	#		$xe = $trajx2[$i+1];
	#		$ye = $trajy2[$i+1];
	#		$ze = $trajz2[$i+1];
	#		if (abs($zs-$ze) < 1
	#			&& abs($xs-$xe) < 1
	#			) {
	#				printf OUT "l $xs $ys $zs $xe $ye $ze\n";
	#			}
	#		#$xs = $trajx2[$i];
	#		#$ys = $trajy2[$i];
	#		#$zs = $trajz2[$i];
	#		#		printf OUT "c $xs $ys $zs\n";
	#	}
	#    printf OUT "@ 4\n";
	#	for ($i = 0; $i < $c_traj; $i++){
	#		$xs = $trajx3[$i];
	#		$ys = $trajy3[$i];
	#		$zs = $trajz3[$i];
	#		$xe = $trajx3[$i+1];
	#		$ye = $trajy3[$i+1];
	#		$ze = $trajz3[$i+1];
	#		if (abs($zs-$ze) < 1
	#			&& abs($xs-$xe) < 1
	#			) {
	#				printf OUT "l $xs $ys $zs $xe $ye $ze\n";
	#			}
	#		#$xs = $trajx3[$i];
	#		#$ys = $trajy3[$i];
	#		#$zs = $trajz3[$i];
	#		#printf OUT "c $xs $ys $zs\n";
	#	}
	
	printf OUT "y 1\n";
    printf OUT "@ 10\n";
	$r = $yap_radius*$radius[0];
	printf OUT "r $r\n";
    for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		if ($y_section == 0 ||
			abs($posy[$i]) < $y_section ){
				printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
			}
		
    }
#	printf OUT "y 3\n";
#	printf OUT "@ 4\n";
#	$r = 1.02*$yap_radius*$radius[0];
#	printf OUT "r $r\n";
#	for ($i = 0; $i < $np; $i ++){
#		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
#			$r = 1.02*$yap_radius*$radius[$i];
#			printf OUT "r $r\n";
#		}
#		#		if ($i % 100 == 0){
#		#			$col = $i/100 + 2;
#		#			printf OUT "@ $col\n";
#		#		}
#		if ($y_section == 0 ||
#			abs($posy[$i]) < $y_section ){
#				$yy = $posy[$i]+0.1;
#				printf OUT "c $posx[$i] $yy $posz[$i] \n";
#			}
#		
#	}
#

	
	printf OUT "y 2\n";
	printf OUT "r 0.2\n";

	printf OUT "@ 5\n"; # static
	for ($k = 0; $k < $num_interaction; $k ++){
		if ($contactstate[$k] > 1) {
			&OutString2($int0[$k],  $int1[$k]);
			#&OutContact($int0[$k], $int1[$k], $contactstate[$k]);
		}
	}
	#
	#
#		printf OUT "y 4\n";
#		printf OUT "r 0.4\n";
#		printf OUT "@ 0\n"; # static
#		for ($k = 0; $k < $num_interaction; $k ++){
#			$force = $Fc_n[$k];
#	        if ($F_lub[$k] < 0) {
#				$force += - $F_lub[$k];
#			}
#			if ($Gap[$k] < 0) {
#				if ($fricstate[$k] == 1 && abs($domega[$k]) < 100) {
#					if ( abs($posx[$int0[$k]]-$posx[$int1[$k]]) < 10
#						&& abs($posz[$int0[$k]]-$posz[$int1[$k]]) < 10){
#							$xx = 0.5*($posx[$int0[$k]] + $posx[$int1[$k]]);
#							$zz = 0.5*($posz[$int0[$k]] + $posz[$int1[$k]]);
#							$tmpoy = int $oy;
#							printf OUT "t $xx -0.2 $zz $tmpoy \n";
	#					}
	#
	#			}
	#		}
	#    }
	#		printf OUT "y 3\n";
	#		printf OUT "r 0.4\n";
	#		printf OUT "@ 7\n"; # static
	#		for ($k = 0; $k < $num_interaction; $k ++){
	#			$force = $Fc_n[$k];
	#	        if ($F_lub[$k] < 0) {
	#				$force += - $F_lub[$k];
	#			}
	#			if ($Gap[$k] < 0) {
	#				if ($fricstate[$k] == 0) {
	#					&OutString2($int0[$k],  $int1[$k]);
	#				}
	#			}
	#	    }
	#
	
	#	printf OUT "y 9\n";
	#	printf OUT "r 0.35\n";
	#	printf OUT "@ 6\n"; # dynamic
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		$force = $Fc_n[$k];
	#        if ($F_lub[$k] < 0) {
	#			$force += - $F_lub[$k];
	#		}
	#		if ($Gap[$k] < 0) {
	#			if ($fricstate[$k] == 2) {
	#				&OutString2($int0[$k],  $int1[$k]);
	#			}
	#		}
	#    }
	#	printf OUT "r 0.2\n";
	#	printf OUT "@ 0\n"; # dynamic
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		$force = $Fc_n[$k];
	#        if ($F_lub[$k] < 0) {
	#			$force += - $F_lub[$k];
	#		}
	#		if ($Gap[$k] < 0) {
	#			if ($fricstate[$k] == 0) {
	#				# static
	#				&OutString2($int0[$k],  $int1[$k]);
	#			}
	#		}
	#    }
	
	
	#printf OUT "r 0.2\n";
	$force_factor = 0.05;
	printf OUT "y 3\n";
	printf OUT "@ 3\n";
	for ($k=0; $k<$num_interaction; $k++){
		#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
		#$force = $Fcol[$k];
		#$force = $Fc_n[$k];
		if ($force[$k] < 0){
			$force = -$force[$k];
			$string_width = ${force_factor}*${force};
			#&OutString2($int0[$k], $int1[$k]);
			&OutString_width($int0[$k], $int1[$k]);
		}
	}
	printf OUT "y 4\n";
	printf OUT "@ 4\n";
	for ($k = 0; $k < $num_interaction; $k ++){
		if ($force[$k] > 0){
			$force = $force[$k];
			$string_width = ${force_factor}*${force};
			#&OutString2($int0[$k], $int1[$k]);
			&OutString_width($int0[$k], $int1[$k]);
		}
	}
	#	printf OUT "y 2\n";
#	printf OUT "@ 6\n";
#	printf OUT "r 0.2\n";
#	for ($k = 0; $k < $num_interaction; $k ++){
##		if ($Gap[$k] < 0) {
##			&OutString2($int0[$k], $int1[$k]);
##		}
#		if ($is_contact[$k] == 1) {
#			&OutContact($int0[$k], $int1[$k]);
#			#&OutString2($int0[$k], $int1[$k]);
#		}
#	}
	#	$stressfactor = 0.005;
	#	printf OUT "y 4\n";
	#	printf OUT "@ 3\n";
	#for ($k = 0; $k < $num_interaction; $k ++){
	#		#$force = $force[$k];
	#	$string_width = $stressfactor*abs($S_bf[$k]);
	#		#&OutString2($int0[$k], $int1[$k]);
	#		&OutString_width($int0[$k], $int1[$k]);
	#}
#	printf OUT "@ 5\n";
#	for ($k = 0; $k < $num_interaction; $k ++){
#		if ($S_bf[$k] > 0){
#			#$force = $force[$k];
#			$string_width = $stressfactor*abs($S_bf[$k]);
#			#&OutString2($int0[$k], $int1[$k]);
#			&OutString_width($int0[$k], $int1[$k]);
#		}
	#	}
	
	#
	#	printf OUT "y 10\n";
	#	printf OUT "@ 5\n";
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($Gap[$k] < 0) {
	#			$xx = $sp_x[$k] + 10*$xi_x[$k];
	#			$yy = $sp_y[$k] + 10*$xi_y[$k];
	#			$zz = $sp_z[$k] + 10*$xi_z[$k];
	#			printf OUT "l $sp_x[$k] $sp_y[$k] $sp_z[$k] $xx $yy $zz\n";
	#
	#		}
	#    }
	#	printf OUT "y 9\n";
	#	printf OUT "@ 4\n";
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($Gap[$k] < 0) {
	#			$xx = $sp_x[$k] - $nrvec_x[$k];
	#			$yy = $sp_y[$k] - $nrvec_y[$k];
	#			$zz = $sp_z[$k] - $nrvec_z[$k];
	#			printf OUT "l $sp_x[$k] $sp_y[$k] $sp_z[$k] $xx $yy $zz\n";
	#		}
	#    }
	#	$ff=1.0/1000;
	#	$stressfactor = 0.00001;

	#	$ff=1.0/1000;
	#	printf OUT "@ 3\n";
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($Gap[$k] < 0) {
	#			$xx = $sp_x[$k] + $nrvec_x[$k];
	#			$yy = $sp_y[$k] + $nrvec_y[$k];
	#			$zz = $sp_z[$k] + $nrvec_z[$k];
	#			printf OUT "l $sp_x[$k] $sp_y[$k] $sp_z[$k] $xx $yy $zz\n";
	#		}
	#    }
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($Gap[$k] < 0) {
	#			$xx = $sp_x[$k] - $ff*$ft_x[$k];
	#			$yy = $sp_y[$k] - $ff*$ft_y[$k];
	#			$zz = $sp_z[$k] - $ff*$ft_z[$k];
	#			printf OUT "l $sp_x[$k] $sp_y[$k] $sp_z[$k] $xx $yy $zz\n";
	#		}
	#    }
	

#	if ($Ly == 0){
#		printf OUT "y 6\n";
#		printf OUT "@ 1\n";
#		for ($i = 0; $i < $np; $i ++){
#			&OutCross($i);
#		}
#	}
	
	&OutBoundaryBox;
	
	#	$maxS=0;
	#for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($maxS < $Sxz_lub[$k]){
	#			$maxS = $Sxz_lub[$k];
	#		}
	#}
	
	#	printf OUT "y 4\n";
	#    printf OUT "@ 5\n";
	#    for ($k = 0; $k < $num_interaction; $k ++){
	#		$string_with = $Sxz_lub[$k]/$maxS;
	#		printf OUT "r ${string_with}\n";
	#		&OutCircle_middle($int0[$k],  $int1[$k]);
	#	}
	#	$zpos = $Lz / 2 + 1;
	#	printf OUT sprintf("t 0 0 %3.2f %2.4f\n", $zpos, $maxS);
	#
	#	printf OUT "y 5\n";
	#	printf OUT "@ 2\n";
	#printf OUT "r 0.3\n";
	#    for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($Gap[$k] < 0){
	#			&OutString($int0[$k],  $int1[$k]);
	#			&OutNvec($k);
	#		}
	#	}
	#	printf OUT2 "\n";
	#
	#
}

sub OutBoundaryBox{
	$x0 = -$Lx/2;
	$x1 = -$Lx/2 + $shear_disp / 2;
	$z1 = $Lz/2;
	$x2 = $Lx/2;
	$z2 = 0;
	$x3 = $Lx/2 - $shear_disp / 2;
	$z3 = -$Lz/2;
	
	printf OUT "y 7\n";
	printf OUT "@ 6\n";
#	printf OUT "r 0.2\n";
#	if ($shear_stress == 0.5) {
#		printf OUT "@ 6\n";
#	} else {
#		printf OUT "@ 3\n";
#	}
	#	printf OUT "l -$lx2 0 0 $lx2 0 0\n";
	#	printf OUT "l $x0 0.01 0 $x1 0.01 $z1\n";
	#	printf OUT "l $x2 0.01 $z2 $x3 0.01 $z3\n";
	
	
	#$yb = 0.1;
	if($Ly == 0){
		$lx2 = $Lx/2+1;
		$ly2 = $Ly/2+1;
		$lz2 = $Lz/2+1;
		#printf OUT "p 4 -$lx2 $yb $lz2 $lx2 $yb $lz2 $lx2 $yb -$lz2 -$lx2 $yb -$lz2\n";
		printf OUT "l -$lx2 0 $lz2    $lx2 0 $lz2\n";
		printf OUT "l -$lx2 0 -$lz2   $lx2 0 -$lz2\n";
		printf OUT "l -$lx2 0 -$lz2  -$lx2 0 $lz2\n";
		printf OUT "l $lx2 0 $lz2   $lx2 0 -$lz2\n";
	} else {
		$lx2 = $Lx/2;
		$ly2 = $Ly/2;
		$lz2 = $Lz/2;
		
		printf OUT "l -$lx2 -$ly2 $lz2    $lx2 -$ly2 $lz2\n";
		printf OUT "l -$lx2 -$ly2 -$lz2   $lx2 -$ly2 -$lz2\n";
		printf OUT "l -$lx2 -$ly2 -$lz2  -$lx2 -$ly2 $lz2\n";
		printf OUT "l  $lx2 -$ly2 $lz2    $lx2 -$ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 $lz2    $lx2 $ly2 $lz2\n";
		printf OUT "l -$lx2 $ly2 -$lz2   $lx2 $ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 -$lz2  -$lx2 $ly2 $lz2\n";
		printf OUT "l  $lx2 $ly2 $lz2    $lx2 $ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 $lz2    -$lx2 -$ly2 $lz2\n";
		printf OUT "l $lx2 $ly2 $lz2    $lx2 -$ly2 $lz2\n";
		printf OUT "l $lx2 $ly2 -$lz2    $lx2 -$ly2 -$lz2\n";
		printf OUT "l -$lx2 $ly2 -$lz2    -$lx2 -$ly2 -$lz2\n";

		
	}
	
}

sub OutNvec {
    ($k) = @_;
	$nx = $nrvec_x[$k];
	$ny = $nrvec_y[$k];
	$nz = $nrvec_z[$k];
	printf OUT2 "$nx $ny $nz\n";
}

sub OutString_width {
    ($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i] - 0.01;
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j] - 0.01;
    $zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2 ;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
		printf OUT "r ${string_width}\n";
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
	}
	
}

sub OutString2{
    ($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i] - 0.02;
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j] - 0.02;
    $zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
	}
}

sub OutContact{
	($i, $j, $cs) = @_;
	$xi = $posx[$i];
	$xj = $posx[$j];
	
	$yi = $posy[$i]-0.022;
	$yj = $posy[$j]-0.022;
	
	$zi = $posz[$i];
	$zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
	
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
		if ( sqrt($sq_dist) > $radius[$i] + $radius[$j]+0.001) {
			$tmp = sqrt($sq_dist);
			printf "$tmp $radius[$i]  $radius[$j] $contactstate[$k]\n";
			exit(1)
		}
		if ($cs == 2) {
			printf OUT "r 0.3\n";
		} elsif ($cs == 3) {
			printf OUT "r 0.1\n";
		}
		$xm = ($radius[$j]*$xi+$radius[$i]*$xj)/($radius[$i]+$radius[$j]);
		$zm = ($radius[$j]*$zi+$radius[$i]*$zj)/($radius[$i]+$radius[$j]);
		$norm = sqrt(($posx[$j]-$posx[$i])**2+($posz[$j]-$posz[$i])**2);
		$vecx = ($posx[$j] - $posx[$i])/$norm;
		$vecz = ($posz[$j] - $posz[$i])/$norm;
		$contact_xi = $xm + 0.3*$vecz;
		$contact_zi = $zm - 0.3*$vecx;
		$contact_xj = $xm - 0.3*$vecz;
		$contact_zj = $zm + 0.3*$vecx;
		printf OUT "s $contact_xi $yi $contact_zi $contact_xj $yj $contact_zj\n";
	}
}


sub OutString {
	($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i];
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j];
    $zj = $posz[$j];
	
	if (abs($xi-$xj) < $radius_max*5
		&&  abs($yi-$yj) < $radius_max*5
		&&  abs($zi-$zj) < $radius_max*5){
			if ( $y_section == 0
				|| abs($yi) < $y_section
				|| abs($yj) < $y_section){
					printf OUT "s $xi $yi $zi $xj $yj $zj\n";
				}
		}
}

sub OutMeter {
	($value, $maxvalue, $xx, $zz) = @_;
	if ( $value > 3*$maxvalue ){
		$value = 0;
	}
	if ( $value > $maxvalue ){
		$value = $maxvalue;
	}
	$metertop = $zz + ($Lz-4)*$value/$maxvalue;
	$xx2 = $xx + 1.5;
	printf OUT "p 4 $xx 0 $zz $xx 0 $metertop $xx2 0 $metertop $xx2 0 $zz\n";
}

sub OutStress {
	($value, $maxvalue) = @_;
	$xx = 0.5*$Lz*$value/$maxvalue;
	$xxTip = $xx + 2;
	$arrowhead = 1;
	$arrowwidth = 0.5;
	$zz0 = $Lz/2+2;

	$zzB1 = $zz0-$arrowwidth;
	$zzB2 = $zz0-$arrowhead;
	$zzT1 = $zz0+$arrowwidth;
	$zzT2 = $zz0+$arrowhead;

	$xxTip = $xx + 2*$arrowhead;
	printf OUT "p 7 -$xx 0 $zzB1 $xx 0 $zzB1 $xx 0 $zzB2 $xxTip 0 $zz0 $xx 0 $zzT2 $xx 0 $zzT1 -$xx 0 $zzT1\n";
	printf OUT "p 7 $xx 0 -$zzB1 -$xx 0 -$zzB1 -$xx 0 -$zzB2 -$xxTip 0 -$zz0 -$xx 0 -$zzT2 -$xx 0 -$zzT1 $xx 0 -$zzT1\n";
}


sub OutCircle_middle {
    ($i, $j) = @_;
	$xi = $posx[$i];
    $yi = $posy[$i];
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j];
    $zj = $posz[$j];
	$ai = $radius[$i];
	$aj = $radius[$j];
    $xc = $posx[$i]*$aj/($ai+$aj) + $posx[$j]*$ai/($ai+$aj);
    $yc = $posy[$i]*$aj/($ai+$aj) + $posy[$j]*$ai/($ai+$aj) -0.01;
    $zc = $posz[$i]*$aj/($ai+$aj) + $posz[$j]*$ai/($ai+$aj);
	
	if (abs($xi-$xj) < $radius_max*5
		&&  abs($yi-$yj) < $radius_max*5
		&&  abs($zi-$zj) < $radius_max*5){
			printf OUT "c $xc $yc $zc\n";
		}
}

sub OutCross {
    ($i) = @_;
	$a = $radius[$i];
    $xi = $posx[$i];
    $yi = $posy[$i] - 0.01;
    $zi = $posz[$i];
    $angle = $ang[$i];
	$ux =  $a*cos($angle);
	$uz = -$a*sin($angle);
	$xa = $xi - $ux;
	$ya = $yi - 0.01;
	$za = $zi - $uz;
	$xb = $xi + $ux;
	$yb = $yi - 0.01;
	$zb = $zi + $uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
	$ux =  $a*sin($angle);
	$uz =  $a*cos($angle);
	$xa = $xi - $ux;
	$ya = $yi - 0.01;
	$za = $zi - $uz;
	$xb = $xi + $ux;
	$yb = $yi - 0.01;
	$zb = $zi + $uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
	
}
