#!/usr/bin/perl

# Usage:
# $ ....
use Math::Trig;
use IO::Handle;
use Getopt::Long;
my $particle_data = $ARGV[0];
my $twodim = 0;
my $pi = atan2(1, 1) * 4;
$printinfo = 1;
if ($printinfo) {
	printf "$particle_data\n";
}
$i = index($particle_data, "par_", 0)+4;
$iD3 = index($particle_data, "D3N");
$j = index($particle_data, ".dat", $i-1);
$name = substr($particle_data, $i, $j-$i);
my $dim = 3;
if ($iD3 eq -1){
	$dim = 2;
}
if ($printinfo) {
	printf "$dim dimension\n";;
}

$interaction_data = "int_${name}.dat";
$rheology_data = "data_${name}.dat";
$filename_yap = "y_$name.yap";
$outputss = "ss_$name.dat";
$outputDataReconstruct = "r_$name.dat";

open (OUT, "> ${filename_yap}");
open (OUTSS, "> ${outputss}");
open (OUTDR, "> ${outputDataReconstruct}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
open (IN_rheo, "< ${rheology_data}");

${sum_fmax} = 0;
${cnt} = 0;

my $friction = 1;
if ($friction == 1) {
	$number_of_header_data = 45;
	$elm_jamming = 37;
} else {
	$number_of_header_data = 44;
	$elm_jamming = 36;
}
for ($i = 0; $i<$number_of_header_data; $i++) {
	$line = <IN_rheo>;
	printf "$line";
}
$i=0;
&readHeader;
&yaplotColor;


$first = 1;
$first_int = 1;
$checkpoint = 1;
$shear_strain_previous = 0;
$shearrate_positive = 1;
$shear_direction = 1;
$cntjamming = 0;
$output = 1;
$ii = 0;
$cnt_history = 0;

while (1) {
	$target_stress = $stress[$ii];
	$shear_strain2 = $strain[$ii];
	$ii ++;
	&InParticles;
	last unless defined $line;
	if ($shearrate_positive > 0) {
		if ($shear_strain < $shear_strain_previous) {
			$shearrate_positive = -1;
		}
		$checkpoint = 0;
	} else {
		if ($shear_strain > $shear_strain_previous) {
			$shearrate_positive = 1;
			$checkpoint = 1;
		} else {
			$checkpoint = 0;
		}
	}
	$shear_strain_previous = $shear_strain;
	&InInteractions;
	if ($reversibility_test) {
		if ($first || $checkpoint == 1) {
			&keepInitialConfig;
		}
	}
	&OutYaplotData;
	#printf "$jamming  $shear_strain\n";
	if (0 && $jamming > 0) {
		$evenodd = $cntjamming % 2;
		printf "even or odd $evenodd, $cntjamming  $jamming  $shear_strain $time\n";
		if ($evenodd == 0) {
			printf "output\n";
			#&OutYaplotData;
			&recordTrajectory;
			if ($cnt_history >= 2) {
				if ($cnt_history != 2) {
					printf OUT "\n";
				} else {
					for ($i = 0; $i < $np; $i ++){
						$xx = $traj_posx[$i][0];
						$zz = $traj_posz[$i][0];
						$rr = sqrt($xx*$xx + $zz*$zz);
						if ($rr > 13 && $rr < 18) {
							$mark[$i] = 1;
						} else {
							$mark[$i] = 0;
						}
					}
				}
				printf OUT "y 1 \n";
				printf OUT "@ 8 \n";
				printf OUT "r 0.1 \n";
				$vmax = 0;
				for ($i = 0; $i < $np; $i ++){
					for ($j = $cnt_history-1; $j < $cnt_history; $j++) {
						$xs[$i] = $traj_posx[$i][$j-1];
						$zs[$i] = $traj_posz[$i][$j-1];
						$xe = $traj_posx[$i][$j];
						$ze = $traj_posz[$i][$j];
						if (abs($xs[$i]-$xe) < 10
							&& abs($zs[$i]-$ze) < 10 ) {
								$vx[$i] = $xe - $xs[$i];
								$vz[$i] = $ze - $zs[$i];
								#  printf OUT "s $xs[$i] 0 $zs[$i]  $xe 0 $ze \n";
								$v = $vx[$i]*$vx[$i] + $vy[$i]*$vy[$i];
								if ($v > $vmax) {
									$vmax = $v;
								}
							}
					}
				}
				
				printf OUT "y 1 \n";
				for ($i = 0; $i < $np; $i ++){
					if ($mark[$i] == 1) {
						printf OUT "@ 12 \n";
					} else {
						printf OUT "@ 10 \n";
					}
					printf OUT "r $radius[$i]\n";
					printf OUT "c $xs[$i] 0.01 $zs[$i]\n";
				}
				if (1) {
					$vfactor = 2.0/sqrt($vmax);
					printf OUT "y 2 \n";
					printf OUT "@ 11 \n";
					printf OUT "r 0.2 \n";
					for ($i = 0; $i < $np; $i ++){
						$xe = $xs[$i] + $vfactor*$vx[$i];
						$ze = $zs[$i] + $vfactor*$vz[$i];
						printf OUT "s $xs[$i] 0 $zs[$i]  $xe 0 $ze \n";
					}
				}
				
			}
		}
		$cntjamming ++;
	}
}

close (OUT);

close (IN_particle);
close (IN_interaction);

sub recordTrajectory {
	for ($i = 0; $i < $np; $i ++){
		$traj_posx[$i][$cnt_history] = $posx[$i];
		$traj_posy[$i][$cnt_history] = $posy[$i];
		$traj_posz[$i][$cnt_history] = $posz[$i];
	}
	$cnt_history++;
}

##################################################################
sub keepInitialConfig {
	for ($i = 0; $i < $np; $i ++){
		$posx_init[$i] = $posx[$i];
		$posy_init[$i] = $posy[$i];
		$posz_init[$i] = $posz[$i];
		$ang_init[$i] = $ang[$i];
		$radius_init[$i] = $radius[$i];
	}
}

sub readHeader {
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $flwtyp) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $dataunit) = split(/\s+/, $line);
	
	if ($Ly == 0) {
		$number_of_header = 9;
	} else {
		$number_of_header = 7;
	}
	for ($i = 0; $i<$number_of_header; $i++) {
		$line = <IN_particle>;
		if ($printinfo) {
			printf "$line";
		}
	}
	printf "---\n";
	if ($Ly == 0) {
		$number_of_header_int = 20;
	} else {
		$number_of_header_int = 20;
	}
	for ($i = 0; $i<$number_of_header_int; $i++) {
		$line = <IN_interaction>;
	}
	#	$xo = $Lx/2;
	#	$yo = $Ly/2;
	#	$zo = $Lz/2;
	$xo = 0;
	$yo = 0;
	$zo = 0;
}

sub yaplotColor {
	#printf OUT "\@1 0 0 0 \n";
	#printf OUT "\@1 50 100 205 \n";
	printf OUT "\@0 25 50 102 \n"; ## bg
	#    printf OUT "\@1 255 255 255  \n"; #bg
	printf OUT "\@1 255 255 255 \n";
	printf OUT "\@2 200 200 200 \n";
	printf OUT "\@3 50 150 255 \n"; # blue
	printf OUT "\@4 50 200 50 \n"; # green
	#printf OUT "\@5 255 100 100 \n"; #red
	printf OUT "\@5 0 0 255 \n"; #blue
	printf OUT "\@6 50 200 50 \n"; # green
	printf OUT "\@7 255 255 0 \n"; # yellow
	#printf OUT "\@7 255 0 0 \n"; # red
	#printf OUT "\@8 255 255 255\n";
	printf OUT "\@8 0 0 0\n";
	printf OUT "\@9 50 50 50\n";
	#printf OUT "\@8 224 143 0 \n";
	#printf OUT "\@9 67 163 230 \n";
	#printf OUT "\@8 253 105 6 \n";
	#printf OUT "\@9 109 109 109 \n";
	printf OUT "\@10 224 240 253 \n";
	printf OUT "\@11 250 50 50 \n";
	printf OUT "\@12 0 0 255 \n";
	printf OUT "\@13 220 220 220 \n";
	printf OUT "\@14 210 210 210 \n";
	printf OUT "\@15 200 200 200 \n";
	printf OUT "\@16 190 190 190 \n";
	printf OUT "\@17 180 180 180 \n";
	printf OUT "\@18 170 170 170 \n";
	printf OUT "\@19 160 160 160 \n";
	printf OUT "\@20 150 150 150 \n";
	printf OUT "\@21 140 140 140 \n";
	printf OUT "\@22 130 130 130 \n";
	printf OUT "\@23 120 120 120 \n";
	printf OUT "\@24 110 110 110 \n";
	printf OUT "\@25 100 100 100 \n";
	printf OUT "\@26 90 90 90 \n";
	printf OUT "\@27 80 90 90 \n";
	printf OUT "\@28 70 70 70\n";
	printf OUT "\@29 60 60 60 \n";
	printf OUT "\@30 50 50 50 \n";
	printf OUT "\@31 40 40 40 \n";
	printf OUT "\@32 30 30 30 \n";
	printf OUT "\@33 20 20 20 \n";
	printf OUT "\@34 10 10 10 \n";
	printf OUT "\@35 0 0 0 \n";
}

sub InParticles {
	$radius_max = 0;
	##  Snapshot Header
	$j = 0;
	
	while (1) {
		$line = <IN_particle>;
		($buf, $val) = split(" : ", $line);
		$val =~ s/(\n|\r)//g;
		($buf1) = split(/\s+/, $buf);
		if ($buf1 ne '#') {
			last;
		} else {
			$ssHeader[$j++] = $val;
		}
		last unless defined $line;
	}
	$shear_strain = $ssHeader[0];
	$shear_disp = $ssHeader[1];
	$shear_rate = $ssHeader[2];
	$target_stress = $ssHeader[3];
	$time = $ssHeader[4];
	@a1 = ($Lx, 0);
	@a2 = ($shear_disp, $Lz);
	
	#$viscosity = $ssHeader[3];
	#$normalstressdiff1 = $ssHeader[4];
	$jamming = 0;
	
	while (1) {
		$line_rheo = <IN_rheo>;
		my @d = split(/\s+/, $line_rheo);
		$time_rheo = $d[0];
		if ($d[$elm_jamming] > 0) {
			$jamming = $d[$elm_jamming];
			$shear_direction *= -1;
		}
		if ($time_rheo >= $time) {
			#	printf "$time_rheo == $time jamming $jamming\n";
			$stressdata = $d[$elm_jamming];
			last;
		}
		last unless defined $line_rheo;;
		#$i++;
	}
	
	for ($i = 0; $i < $np; $i ++){
		if ($i > 0) {
			$line = <IN_particle>;
		}
		if ($output == 1) {
			# LF_DEM version v3.0-53-g9f8a2585-dirty
			# np 500
			# VF 0.78
			# Lx 51.6458
			# Ly 0
			# Lz 51.6458
			# flow_type shear
			# data in r units.
			#1: particle index
			#2: radius
			#3: position x
			#4: position z
			#5-7: velocity (x, y, z)
			#8-10: angular velocity (x, y, z)
			#11: angle
			#				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			#	$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			if ($dim eq 3) {
				($ip, $a, $x, $y, $z, $vx, $vz, $vy, $ox, $oz, $oy, $connum) = split(/\s+/, $line);
			} else {
				($ip, $a, $x, $z, $vx, $vz, $vy, $ox, $oz, $oy, $angle, $connum) = split(/\s+/, $line);
			}
			#
			#
			$ang[$i] = $angle;
			$radius[$i] = $a;
			if ($xz_shift) {
				$x += $Lx/2;
				$z += $Lz/2;
				if ($x > $Lx/2) {
					$x -= $Lx;
				}
				if ($z > $Lz/2) {
					$z -= $Lz;
				}
			}
			$posx[$i] = $x-$xo;
			$posy[$i] = $y-$yo;
			$posz[$i] = $z-$zo;
			$velx[$i] = $vx;
			$vely[$i] = $vy;
			$velz[$i] = $vz;
			$omegax[$i] = $ox;
			$omegay[$i] = $oy;
			$omegaz[$i] = $oz;
			$contactnumber[$i] = $connum;
			if ($radius_max < $a) {
				$radius_max = $a;
			}
		}
	}
}

sub InInteractions{
	#$line = <IN_interaction>;
	#($buf, $shear_strain_i, $num_interaction) = split(/\s+/, $line);
	#printf "int $buf $shear_strain_i $num_interaction\n";
	while (1) {
		$line = <IN_interaction>;
		($buf, $val) = split(" : ", $line);
		($buf1) = split(/\s+/, $buf);
		if ($buf1 ne '#') {
			last;
		} else {
			$val =~ s/(\n|\r)//g;
			$ssHeader[$j++] = $val;
		}
		last unless defined $line;
	}
	
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
	$k = 0;
	while (true) {
		if ($k > 0) {
			$line = <IN_interaction>;
		}
		#1: particle 1 label
		#2: particle 2 label
		#3: contact state (0 = no contact, 1 = frictionless contact, 2 = non-sliding frictional, 3 = sliding frictional)
		#4-6: normal vector, oriented from particle 1 to particle 2
		#7: dimensionless gap = s-2, s = 2r/(a1+a2)
		#8: normal part of the lubrication force (positive for compression)
		#9-11: tangential part of the lubrication force
		#12: norm of the normal part of the contact force
		#13-15: tangential part of the contact force
		#16: norm of the normal repulsive force
		#17: Viscosity contribution of contact xF
		
		($i, $j, $contact, $nx, $nz, $ny, #1---6
		$gap, $f_lub_norm, # 7, 8
		$f_lub_tan_x, $f_lub_tan_z, $f_lub_tan_y, # 9, 10, 11
		$fc_norm, # 12
		$fc_tan_x, $fc_tan_z, $fc_tan_y, # 13, 14, 15
		$fr_norm, $s_xF) = split(/\s+/, $line); # 16, 17
		
		if ($i eq '#' || $i eq NU) {
			($buf, $val) = split(" : ", $line);
			$val =~ s/(\n|\r)//g;
			$shear_strain = $val;
			last;
		}
		if (! defined $i) {
			last;
		}
		if ($output == 1) {
			$int0[$k] = $i;
			$int1[$k] = $j;
			$contactstate[$k] = $contact;
			#            if ($contact > 0) {
			#                #$force[$k] = $fc_norm + $f_lub_norm + $fr_norm;
			#                $force[$k] = $fc_norm + $f_lub_norm + $fr_norm;
			#                #$force[$k] = $fc_norm;
			#            } else {
			#                $force[$k] = $f_lub_norm + $fr_norm;
			#                #$force[$k] = 0;
			#            }
			$F_lub[$k] = $f_lub_norm;
			$Fc_n[$k] = $fc_norm;
			$Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
			$F_lub_t[$k] = sqrt($f_lub_tan_x**2+ $f_lub_tan_z**2 + $f_lub_tan_y**2);
			$S_bf[$k] =  $s_xF;
			#$force[$k] += $Fc_t[$k] + $F_lub_t[$k];
			$nrvec_x[$k] = $nx;
			$nrvec_y[$k] = $ny;
			$nrvec_z[$k] = $nz;
			$Gap[$k] = $gap;
			$distance[$k] = $radius[$i] + $radius[$j] + $gap;
			$k++;
		}
	}
	$num_interaction = $k;
}

sub OutYaplotData{
	if ($first == 0) {
		printf OUT "\n";
	} else {
		$first = 0;
		if (0) {
			for ($i = 0; $i < $np; $i ++){
				$xx = $posx[$i];
				$zz = $posz[$i];
				$rr = sqrt($xx*$xx + $zz*$zz);
				if ($rr > 13 && $rr < 18) {
					$mark[$i] = 1;
				} else {
					$mark[$i] = 0;
				}
			}
		}
	}
	#	printf OUT "@ 7\n";
	# &OutStress($shear_direction*$target_stress, 5);
	
	printf OUT "y 1\n";
	#	printf OUT "@ 2\n";
	printf OUT "@ 8\n";
	## visualize particles
	if ($monodisperse) {
		
		printf OUT "r $radius[0]\n";
		for ($i = 0; $i < $np; $i++) {
			printf OUT "c $posx[$i] $posy[$i] $posz[$i]  \n";
		}
	} else {
		for ($i = 0; $i < $np; $i++) {
			$rr = $yap_radius*$radius[$i];
			if (1) {
				#if ($mark[$i] == 1) {
				if ($contactnumber[$i] >= 2) {
					printf OUT "@ 8 \n";
				} else {
					printf OUT "@ 1 \n";
				}
			}
			printf OUT "r $rr\n";
			printf OUT "c $posx[$i] $posy[$i] $posz[$i]  \n";
		}
	}
	if (0) {
		printf OUT "y 3\n";
		printf OUT "@ 0\n";
		$npjamming = 0;
		$totalcontact = 0;
		for ($i = 0; $i < $np; $i++) {
			#		if ($contactnumber[$i] == 0) {
			#	$rr = $yap_radius*$radius[$i];
			#	printf OUT "r $rr\n";
			printf OUT "t $posx[$i] -0.1 $posz[$i] $contactnumber[$i] \n";
			if ($contactnumber[$i] >= 2) {
				$npjamming ++;
				$totalcontact += $contactnumber[$i];
			}
		}
		if ($npjamming > 0) {
			$coordinationnumber = $totalcontact / $npjamming;
			printf "$npjamming $totalcontact  : z = $coordinationnumber \n";
			$zt = 0.53*$Lz;
			printf OUT "@ 8\n";
			printf OUT "r 20\n";
			printf OUT "t 0 -0.1 $zt z = $coordinationnumber\n";
		}
	}
	
	#	printf OUT "y 2\n";
	#	printf OUT "@ 0\n";
	#	## visualize particles
	#	printf OUT "r 0.1\n";
	#	for ($i = 0; $i < $np; $i++) {
	#		printf OUT "c $posx[$i] -0.01 $posz[$i]  \n";
	#	}
	
	## visualize contact network
	if (0) {
		printf OUT "y 2\n";
		printf OUT "r 0.2\n";
		printf OUT "@ 6\n"; # static
		for ($k = 0; $k < $num_interaction; $k ++) {
			if ($contactstate[$k] >= 2) {
				$i = $int0[$k];
				$j = $int1[$k];
				if ($contactnumber[$i] >= 2 && $contactnumber[$j] >= 2) {
					&OutString2($int0[$k], $int1[$k]);
				}
			}
		}
		if (0) {
			printf OUT "y 2\n";
			printf OUT "r 0.2\n";
			printf OUT "@ 2\n"; # static
			for ($k = 0; $k < $num_interaction; $k ++) {
				if ($contactstate[$k] == 1) {
					$i = $int0[$k];
					$j = $int1[$k];
					if ($contactnumber[$i] >= 2 && $contactnumber[$j] >= 2) {
						&OutString2($int0[$k], $int1[$k]);
					}
				}
			}
		}
	}
	# visualize force chain network
	#
	if (1) {
		printf OUT "y 4\n";
		printf OUT "@ 2\n";
		for ($k = 0; $k < $num_interaction; $k ++) {
			$i = $int0[$k];
			$j = $int1[$k];
			if ($contactstate[$k] >= 2) {
				if ($contactnumber[$i] >= 2 && $contactnumber[$j] >= 2) {
					#$force = $F_lub[$k] + $Fc_n[$k];
					$force = sqrt($Fc_n[$k]**2 + $Fc_t[$k]**2);
					&OutString_width($int0[$k], $int1[$k], $force_factor*abs($force), -0.01);
				}
				#&OutString_width($int0[$k], $int1[$k], $force_factor*abs($force)/$target_stress, 0.01);
			}
		}
	}
	if (0) {
		printf OUT "@ 5\n";
		for ($k = 0; $k < $num_interaction; $k ++) {
			$force = $F_lub[$k] + $Fc_n[$k];
			if ($force < 0) {
				&OutString_width($int0[$k], $int1[$k], $force_factor*abs($force), 0.05);
			}
		}
	}
	#    if ($cnt eq $confout) {
	#        if (0) {
	#            for ($k = 0; $k < $num_interaction; $k ++) {
	#                $force = $F_lub[$k] + $Fc_n[$k];
	#                $i = $int0[$k];
	#                $j = $int1[$k];
	#                $xi = $posx[$i];
	#                $zi = $posz[$i];
	#                $xj = $posx[$j];
	#                $zj = $posz[$j];
	#                $sq_dist = ($xi-$xj)**2 + ($zi-$zj)**2;
	#                $min = ($radius[$i] + $radius[$j]+1)**2;
	#                if ($sq_dist < $min) {
	#                    printf "$xi $zi $xj $zj $force\n";
	#                }
	#            }
	#            exit;
	#        }
	#    }
	if (0){
		
		printf OUT "y 6\n";
		printf OUT "@ 7\n";
		$etaTotal = 0;
		for ($k = 0; $k < $num_interaction; $k ++) {
			#$force = $F_lub[$k] + $Fc_n[$k];
			$nx = $nrvec_x[$k];
			$ny = $nrvec_y[$k];
			$eta_lub = -$F_lub[$k]*$distance[$k]*($nx*$ny);
			$eta_con = -$Fc_n[$k]*$distance[$k]*($nx*$ny);
			$eta[$k] = $eta_lub + $eta_con;
			$etaTotal += $eta[$k];
			if ($eta[$k] > 0) {
				printf OUT "@ 7\n";
			} else {
				printf OUT "@ 5\n";
			}
			&OutString_width($int0[$k], $int1[$k],  $force_factor*abs($eta[$k]), 0.01);
			#        }
		}
		$visapprox = $etaTotal/($Lx*$Lz);
	}
	#    &calcContributions;
	#    for ($k = 0; $k < $num_interaction; $k ++) {
	#    $force = $F_lub[$k] + $Fc_n[$k];
	#    if ($force > 0) {
	#        printf OUT "@ 7\n";
	#    } else {
	#        printf OUT "@ 5\n";
	#    }
	#    &OutString_width($int0[$k], $int1[$k], $force_factor*abs($force[$k]), 0.02);
	#    #&OutString_width($int0[$k], $int1[$k], $force_factor*abs($n1), 0.01);
	# }
	
	if (0) {
		printf OUT "y 3\n";
		printf OUT "@ 5\n";
		for ($k = 0; $k < $num_interaction; $k ++) {
			#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
			if ($force[$k] < 0) {
				&OutString_width($int0[$k], $int1[$k], -$force_factor*$force[$k], 0.02);
			}
		}
	}
	## visualize rotation in 2D
	if ($Ly == 0) {
		if (0) {
			printf OUT "y 6\n";
			printf OUT "@ 8\n";
			for ($i = 0; $i < $np; $i++) {
				OutCross($i);
			}
		}
	}
	if ($reversibility_test) {
		printf OUT "y 5\n";
		printf OUT "@ 7\n";
		$r = $radius_init[0];
		printf OUT "r $r\n";
		for ($i = 0; $i < $np; $i++) {
			if ($i >= 1 && $radius_init[$i] != $radius_init[$i-1]) {
				$r = $yap_radius*$radius_init[$i];
				printf OUT "r $r\n";
			}
			if ($Ly == 0) {
				printf OUT "c $posx_init[$i] 0.01 $posz_init[$i] \n";
			} else {
				printf OUT "c $posx_init[$i] $posy_init[$i] $posz_init[$i] \n";
			}
		}
	}
	&OutBoundaryBox;
}

sub OutBoundaryBox {
	$x0 = -$Lx/2;
	$x1 = -$Lx/2 + $shear_disp / 2;
	$z1 = $Lz/2;
	$x2 = $Lx/2;
	$z2 = 0;
	$x3 = $Lx/2 - $shear_disp / 2;
	$z3 = -$Lz/2;
	
	printf OUT "y 7\n";
	printf OUT "@ 6\n";
	
	if ($Ly == 0) {
		$lx2 = $Lx/2;
		$ly2 = $Ly/2;
		$lz2 = $Lz/2;
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
	if ($axis) {
		printf OUT "l -$lx2 0 0 $lx2 0 0\n";
		printf OUT "l 0 0 -$lz2 0 0 $lz2\n";
	}
}

sub OutString_width {
	($i, $j, $w, $delta) = @_;
	$xi = $posx[$i];
	$yi = $posy[$i] - $delta;
	$zi = $posz[$i];
	$xj = $posx[$j];
	$yj = $posy[$j] - $delta;
	$zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2 ;
	printf OUT "r $w\n";
	my $max_sq_dist = ($radius[$i] + $radius[$j]+1)**2;
	if ($sq_dist < $max_sq_dist) {
		printf OUT "s $xi $yi $zi $xj $yj $zj \n";
	} else {
		for ($ix = -1; $ix <=1; $ix ++) {
			for ($iz = -1; $iz <=1; $iz ++) {
				$xip = $xi + $ix*$a1[0] + $iz*$a2[0];
				$zip = $zi + $ix*$a1[1] + $iz*$a2[1];
				$sq_dist = ($xip-$xj)**2 + ($zip-$zj)**2;
				if ($sq_dist < $max_sq_dist) {
					printf OUT "s $xip $yi $zip $xj $yj $zj \n";
					$i_done = 1;
					$xjp = $xj - $ix*$a1[0] - $iz*$a2[0];
					$zjp = $zj - $ix*$a1[1] - $iz*$a2[1];
					$sq_dist = ($xi-$xjp)**2 + ($zi-$zjp)**2;
					printf OUT "s $xi $yi $zi $xjp $yj $zjp \n";
					last;
				}
			}
		}
	}
}

sub OutString2{
	$offset2 = -0.012;
	($i, $j) = @_;
	$xi = $posx[$i];
	$yi = $posy[$i]+$offset2;
	$zi = $posz[$i];
	$xj = $posx[$j];
	$yj = $posy[$j]+$offset2;
	$zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1) {
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
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
		&&  abs($zi-$zj) < $radius_max*5) {
			printf OUT "s $xi $yi $zi $xj $yj $zj\n";
		}
}

sub OutCross {
	$offset1 = -0.01;
	($i) = @_;
	$a = $radius[$i];
	$xi = $posx[$i];
	$yi = $posy[$i]+$offset1;
	$zi = $posz[$i];
	$angle = $ang[$i];
	$ux =  $a*cos($angle);
	$uz = -$a*sin($angle);
	$xa = $xi-$ux;
	$ya = $yi+$offset1;
	$za = $zi-$uz;
	$xb = $xi+$ux;
	$yb = $yi+$offset1;
	$zb = $zi+$uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
	$ux = $a*sin($angle);
	$uz = $a*cos($angle);
	$xa = $xi-$ux;
	$ya = $yi+$offset1;
	$za = $zi-$uz;
	$xb = $xi+$ux;
	$yb = $yi+$offset1;
	$zb = $zi+$uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
}

sub calcContributions {
	$n1approx_pos_con = 0;
	$n1approx_neg_con = 0;
	$n1approx_pos_lub = 0;
	$n1approx_neg_lub = 0;
	
	#    $n1_factor = 0.015/$viscosity;
	$n1_factor = $force_factor;
	printf OUT "y 5\n";
	#    printf OUT "@ 7\n";
	for ($k = 0; $k < $num_interaction; $k ++) {
		#$force = $F_lub[$k] + $Fc_n[$k];
		$nx = $nrvec_x[$k];
		$ny = $nrvec_y[$k];
		$n1lub = -$F_lub[$k]*$distance[$k]*($nx*$nx - $ny*$ny);
		$n1con = -$Fc_n[$k]*$distance[$k]*($nx*$nx - $ny*$ny);
		#        $force = $Fc_n[$k];
		#if (1 || $force[$k] >= 0) {
		#    &OutString_width($int0[$k], $int1[$k], $force_factor*$force[$k], 0.01);
		#}
		#        if ($force > 0) {
		
		if ($n1lub > 0) {
			$n1approx_pos_lub += $n1lub;
		} else {
			$n1approx_neg_lub += $n1lub;
		}
		if ($n1con > 0) {
			$n1approx_pos_con += $n1con;
		} else {
			$n1approx_neg_con += $n1con;
		}
		$n1[$k] = $n1lub + $n1con;
		if ($n1[$k] > 0) {
			printf OUT "@ 7\n";
		} else {
			printf OUT "@ 5\n";
		}
		# &OutString_width($int0[$k], $int1[$k], $force_factor*$force[$k], 0.01);
		&OutString_width($int0[$k], $int1[$k], $n1_factor*abs($n1[$k]), 0.01);
		#        }
	}
	
	if ($cnt eq $confout) {
		for ($k = 0; $k < $num_interaction; $k ++) {
			$force = $F_lub[$k] + $Fc_n[$k];
			$i = $int0[$k];
			$j = $int1[$k];
			$xi = $posx[$i];
			$zi = $posz[$i];
			$xj = $posx[$j];
			$zj = $posz[$j];
			$sq_dist = ($xi-$xj)**2 + ($zi-$zj)**2;
			$min = ($radius[$i] + $radius[$j]+1)**2;
			if ($sq_dist < $min) {
				printf "$xi $zi $xj $zj $force $n1[$k] $eta[$k]\n";
			}
		}
		for ($i = 0; $i < $np; $i++) {
			printf OUTSS "$posx[$i]  $posz[$i] $radius[$i] \n";
		}
		exit;
	}
	
	#    $barsize = 0.1;
	$barsize = 10;
	$n1approxP = $n1approx_pos_con/($Lx*$Lz*$viscosity);
	$n1approxN = $n1approx_neg_con/($Lx*$Lz*$viscosity);
	$testposition1 = $Lx/2 + 1;
	$testposition2 = $Lx/2 + 3;
	
	printf OUT "@ 7\n";
	printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $n1approxP $testposition1 0 $n1approxP \n";
	printf OUT "@ 5\n";
	printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $n1approxN $testposition1 0 $n1approxN \n";
	
	$n1approxCon = $n1approxP + $n1approxN;
	$testposition1 = $Lx/2 + 3.2;
	$testposition2 = $Lx/2 + 5.2;
	
	if ($n1approx > 0) {
		printf OUT "@ 7\n";
	} else {
		printf OUT "@ 5\n";
	}
	printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $n1approxCon $testposition1 0 $n1approxCon \n";
	
	$testposition1 = -$Lx/2 - 3.2;
	$testposition2 = -$Lx/2 - 5.2;
	
	$n1approxP = $n1approx_pos_lub/($Lx*$Lz*$viscosity);
	$n1approxN = $n1approx_neg_lub/($Lx*$Lz*$viscosity);
	printf OUT "@ 7\n";
	printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $n1approxP $testposition1 0 $n1approxP \n";
	printf OUT "@ 5\n";
	printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $n1approxN $testposition1 0 $n1approxN \n";
	
	$n1approxLub = $n1approxP + $n1approxN;
	$testposition1 = -$Lx/2 - 1;
	$testposition2 = -$Lx/2 - 3;
	if ($n1approx > 0) {
		printf OUT "@ 7\n";
	} else {
		printf OUT "@ 5\n";
	}
	printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $n1approxLub $testposition1 0 $n1approxLub \n";
	
	#    $ratio = $normalstressdiff1/$n1approx;
	#    printf "ratio = $ratio\n";
	#    $n1approx = 20*$n1approx/$viscosity;
	#$totaln1 = 20*$normalstressdiff1/$viscosity;
	#    if ($normalstressdiff1 > 0 ) {
	#    printf OUT "@ 7\n";
	#} else {
	#    printf OUT "@ 5\n";
	#}
	#printf OUT "r 1\n";
	#printf OUT "p 4 $testposition1 0 0 $testposition2 0 0 $testposition2 0 $totaln1 $testposition1 0 $totaln1 \n";
	
	if ($shear_strain > 1) {
		$n1approx = ($n1approxCon + $n1approxLub);
		####
		$viscosityOUT = 6*$pi*$viscosity;
		$visapproxOUT = 6*$pi*$visapprox;
		$normalstressdiff1OUT = $normalstressdiff1/$viscosity;
		$n1approxOUT = $n1approx;
		printf OUTDR "$shear_strain $viscosityOUT $visapproxOUT $normalstressdiff1OUT $n1approxOUT\n";
	}
}
sub OutStress {
	($value, $maxvalue) = @_;
	
	$xx = 0.5*$Lz*abs($value)/$maxvalue;
	$xxTip = $xx + 2;
	$arrowhead = 2;
	$arrowwidth = 1;
	$zz0 = $Lz/2+2;
	
	$zzB1 = $zz0-$arrowwidth;
	$zzB2 = $zz0-$arrowhead;
	$zzT1 = $zz0+$arrowwidth;
	$zzT2 = $zz0+$arrowhead;
	
	$xxTip = $xx + 2*$arrowhead;
	if ($value > 0) {
		printf OUT "p 7 -$xx 0 $zzB1 $xx 0 $zzB1 $xx 0 $zzB2 $xxTip 0 $zz0 $xx 0 $zzT2 $xx 0 $zzT1 -$xx 0 $zzT1\n";
		printf OUT "p 7 $xx 0 -$zzB1 -$xx 0 -$zzB1 -$xx 0 -$zzB2 -$xxTip 0 -$zz0 -$xx 0 -$zzT2 -$xx 0 -$zzT1 $xx 0 -$zzT1\n";
	} else {
		printf OUT "p 7 $xx 0 $zzB1 -$xx 0 $zzB1 -$xx 0 $zzB2 -$xxTip 0 $zz0 -$xx 0 $zzT2 -$xx 0 $zzT1 $xx 0 $zzT1\n";
		printf OUT "p 7 -$xx 0 -$zzB1 $xx 0 -$zzB1 $xx 0 -$zzB2 $xxTip 0 -$zz0 $xx 0 -$zzT2 $xx 0 -$zzT1 -$xx 0 -$zzT1\n";
	}
	
}
