#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [output_interval] [xz_shift]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;
use Getopt::Long;

my $particle_data = $ARGV[0];
my $yap_radius = 1;
my $force_factor = 0.0001;
my $output_interval = 1;
my $xz_shift = 1;
my $axis = 0;
my $reversibility_test = 0;
my $monodisperse = 0;
my $epsilondot=0;
my $timeper=0;
my $scale=0.5;
my $axisrotation=1;
my $pd_color = 0;
my $draw_cross = 0;
my $flow_type = "shear";
my $draw_trajectory = 0;
my $phi6_data = 1;

my $ii_min = -1;
my $ii_max = 1;
my $jj_min = -2;
my $jj_max = 1;

GetOptions(
'forcefactor=f' => \$force_factor,
'interval=i' => \$output_interval,
'shift=f' => \$xz_shift,
'axis' => \$axis,
'reversibility' => \$reversibility_test,
'monodisperse' => \$monodisperse);

printf "force_factor = $force_factor\n";
printf "output_interval = $output_interval\n";
printf "xz_shift = $xz_shift\n";
printf "axis = $axis\n";
printf "reversibility = $reversibility_test\n";
printf "monodisperse = $monodisperse\n";

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$interaction_data = "int_${name}.dat";
$outputname = "y_$name.yap";

open (OUT, "> ${outputname}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
$cnt2 = 0;
&readHeader;
#&yaplotColor;
&yapColor;

$tcnt = 0;
$cnt_interval = 0;
$first = 1;
$first_int = 1;
$checkpoint = 1;
$shear_strain_previous = 0;
$shearrate_positive = 1;
$particlecolor = 3;
$i_trac = 3;
$strainIntFile = -1;

my $magicangle = atan(0.5*(sqrt(5)-1));
my $cos_ma = cos($magicangle);
my $sin_ma = sin($magicangle);

while (1) {
	if ($cnt_interval == 0 ||
		$cnt_interval % $output_interval == 0) {
			$output = 1;
		} else {
			$output = 0;
		}
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
	if (0) {
		&InInteractions;
	}
	
	if ($reversibility_test) {
		if ($first || $checkpoint == 1) {
			&keepInitialConfig;
		}
	}
	
	if ($output == 1) {
		&OutYaplotData;
	}
	$cnt_interval ++;
}

close (OUT);

close (IN_particle);
close (IN_interaction);

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
	$line = <IN_particle>; ($buf, $buf, $flow_type) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $dataunit) = split(/\s+/, $line);
	
	if ($Ly == 0) {
		$number_of_header = 9;
	} else {
		$number_of_header = 7;
	}
	printf "=====\n";
	for ($i = 0; $i<$number_of_header; $i++) {
		$line = <IN_particle>;
		printf "$line";
	}
	printf "=====\n";
	for ($i = 0; $i<20; $i++) {
		$line = <IN_interaction>;
		printf "$line";
	}
	printf "=====\n";

}

sub InParticles {
	
	$radius_max = 0;
	##  Snapshot Header
	$j = 0;
	while (1) {
		$line = <IN_particle>;
		($buf, $val) = split(" : ", $line);
		($buf1) = split(/\s+/, $buf);
		if ($buf1 ne '#') {
			last;
		} else {
			$ssHeader[$j++] = $val;
		}
		last unless defined $line;
	}
	$shear_strain = $ssHeader[0];
	printf "strain = $shear_strain";
	$shear_disp = $ssHeader[1];
	$epsilondot = $ssHeader[2];
	$timeper = $ssHeader[3]/2;
	$retrim = $ssHeader[4];
	for ($i = 0; $i < $np; $i ++){
		if ($i > 0) {
			$line = <IN_particle>;
		}
		if (1) {
			# 1: number of the particle
			# 2: radius
			# 3, 4, 5: position
			# 6, 7, 8: velocity
			# 9, 10, 11: angular velocity
			# 12: viscosity contribution of lubrication
			# 13: viscosity contributon of contact GU xz
			# 14: viscosity contributon of brownian xz
			# (15: angle for 2D simulation)
			#                ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			#    $h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			
			# 3D
			# ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz) = split(/\s+/, $line);
			# 2D
			if ($dim eq 3) {
				($ip, $a, $x, $y, $z, $vx, $vz, $vy, $ox, $oz, $oy, $connum) = split(/\s+/, $line);
			} else {
				($ip, $a, $x, $z, $vx, $vz, $vy, $ox, $oz, $oy, $angle, $connum) = split(/\s+/, $line);
			}
			#
			$ang[$i] = $angle;
			$radius[$i] = $a;
			if ($xz_shift) {
				$x += $Lx/2;
				$z += $Lz/2;
				#if ($x > $Lx/2) {
				#    $x -= $Lx;
				# }
				#if ($z > $Lz/2) {
				#    $z -= $Lz;
				#}
			}
			$posx[$i] = $x;
			$posy[$i] = $y;
			$posz[$i] = $z;
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
	#	$line = <IN_interaction>;
	#    ($buf, $shear_strain_i, $num_interaction) = split(/\s+/, $line);
	# printf "int $buf $shear_strain_i $num_interaction\n";
	
	if ($strainIntFile < 0) {
		$lineInt = <IN_interaction>;
		if ( $strainIntFile == -1) {
			($buf, $strainIntFile) = split(" : ", $lineInt);
			printf "ssint = $strainIntFile\n" ;
		} else {
			($buf, $val) = split(" : ", $lineInt);
		}
		
		while (1) {
			$lineInt = <IN_interaction>;
			($buf, $val) = split(" : ", $lineInt);
			($buf1) = split(/\s+/, $buf);
			if ($buf1 ne '#') {
				last;
			} else {
				$ssHeader[$j++] = $val;
			}
			last unless defined $lineInt;
		}
	}
	if ($strainIntFile <= $shear_strain || $strainIntFile == -2) {
		$strainIntFile = -2;

		$k = 0;
		while (true) {
			if ($k > 0) {
				$lineInt = <IN_interaction>;
			}
			($i, $j, $contact, $nx, $ny, $nz, #1---6
			$gap, $f_lub_norm, # 7, 8
			$f_lub_tan_x, $f_lub_tan_y, $f_lub_tan_z, # 9, 10, 11
			$fc_norm, # 12
			$fc_tan_x, $fc_tan_y, $fc_tan_z, # 13, 14, 15
			$fr_norm, $s_xF) = split(/\s+/, $lineInt);
			#			printf "k = $k \n" ;
			if ($i eq '#' || $i eq NU) {
				#   ($buf, $shear_strain) = split(" : ", $lineInt);
				$val =~ s/(\n|\r)//g;
				#				$shear_strain = $val;
				printf "A $shear_strain \n" ;
				last;
			}
			if (! defined $i) {
				printf "$lineInt";
				printf "B $k $i \n" ;
				last;
			}
			#		last unless defined $line;
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
			if ($output == 1) {
				$int0[$k] = $i;
				$int1[$k] = $j;
				$contactstate[$k] = $contact;
				$F_lub[$k] = $f_lub_norm;
				#$F_lub_tan = sqrt($f_lub_tan_x**2 + $f_lub_tan_y**2 + $f_lub_tan_z**2);
				
				$Fc_n[$k] = $fc_norm;
				#$Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
				
				if ($contact > 0) {
					#$force[$k] = $fc_norm + $f_lub_norm + $fr_norm;
					$fn = $fc_norm + $f_lub_norm + $fr_norm;
					#$force[$k] = $fc_norm;
				} else {
					$fn = $fc_norm + $f_lub_norm + $fr_norm;
					#$force[$k] = $f_lub_norm + $fr_norm;
					#$force[$k] = 0;
				}
				$ftx = $f_lub_tan_x + $fc_tan_x;
				$fty = $f_lub_tan_y + $fc_tan_y;
				$ftz = $f_lub_tan_z + $fc_tan_z;
				$ft = sqrt($ftx*$ftx + $fty*$fty + $ftz*$ftz);
				$force[$k] = sqrt($fn*$fn + $ft*$ft);
				$S_bf[$k] =  $s_xF;
				$nrvec_x[$k] = $nx;
				$nrvec_y[$k] = $ny;
				$nrvec_z[$k] = $nz;
				$Gap[$k] = $gap;
			}
			$k++;
		}
		printf "k=$k\n" ;
		$num_interaction = $k;
	}
}

sub OutYaplotData{
	if ($first == 0) {
		printf OUT "\n";
	} else {
		$first = 0;
	}
	printf OUT "y 1\n";
	printf OUT "@ 8\n";
	## visualize particles
	$epsilondot = 1;
	$ax = exp( $epsilondot*($timeper))*$cos_ma*$cos_ma+exp(-$epsilondot*($timeper))*$sin_ma*$sin_ma;
	$az = exp(-$epsilondot*($timeper))*$cos_ma*$sin_ma-exp( $epsilondot*($timeper))*$sin_ma*$cos_ma;
	$bx = exp(-$epsilondot*($timeper))*$cos_ma*$sin_ma-exp( $epsilondot*($timeper))*$sin_ma*$cos_ma;
	$bz = exp(-$epsilondot*($timeper))*$cos_ma*$cos_ma+exp( $epsilondot*($timeper))*$sin_ma*$sin_ma;
	#        if ($cnt2 == 10) {
	#            exit;
	#        }
	#        $cnt2 += 1;
	$lx2 = $Lx/2;
	$ly2 = $Ly/2;
	$lz2 = $Lz/2;
	if ($flow_type == "extension") {
		&OutParticleExtension;
	} else {
		&OutParticleShear;
	}
	if ($draw_trajectory) {
		for ($ii = $ii_min; $ii <= $ii_max; $ii++) {
			for ($jj = $jj_min; $jj <= $jj_max; $jj++) {
				$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
				$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
				$xx = $posx[$i_trac] + $pd_xshift;
				$zz = $posz[$i_trac] + $pd_zshift;
				if ($axisrotation) {
					$xx2 = $cos_ma*$xx - $sin_ma*$zz;
					$zz2 = $sin_ma*$xx + $cos_ma*$zz;
				} else {
					$xx2 = $xx;
					$zz2 = $zz;
				}
				if (abs($xx2)<$scale*$Lx && abs($zz2) < $scale*$Lz) {
					$tx[$tcnt] = $xx2;
					$tz[$tcnt] = $zz2;
					$tcnt ++
				}
			}
		}
		printf OUT "y 7\n";
		printf OUT "r 0.3\n";
		printf OUT "@ 0\n";
		for ($k = 0; $k < $tcnt; $k ++) {
			printf OUT "c $tx[$k] -0.01 $tz[$k] \n"
		}
	}
	## visualize contact network
	#	if (1) {
	#		printf OUT "y 2\n";
	#		printf OUT "r 0.2\n";
	#		printf OUT "@ 6\n"; # static
	#		for ($k = 0; $k < $num_interaction; $k ++) {
	#			if ($contactstate[$k] >= 2) {
	#				$i = $int0[$k];
	#				$j = $int1[$k];
	#				if ($contactnumber[$i] >= 2 && $contactnumber[$j] >= 2) {
	#					&OutString2($int0[$k], $int1[$k]);
	#				}
	#			}
	#		}
	#		printf OUT "y 2\n";
	#		printf OUT "r 0.2\n";
	#		printf OUT "@ 2\n"; # static
	#		for ($k = 0; $k < $num_interaction; $k ++) {
	#			if ($contactstate[$k] == 1) {
	#				$i = $int0[$k];
	#				$j = $int1[$k];
	#				if ($contactnumber[$i] >= 2 && $contactnumber[$j] >= 2) {
	#					&OutString2($int0[$k], $int1[$k]);
	#				}
	#			}
	#		}
	#	}
	
	#	printf OUT "y 2\n";
	#	printf OUT "r 0.2\n";
	#	printf OUT "@ 0\n"; # static
	#	for ($k = 0; $k < $num_interaction; $k ++) {
	#		if ($Gap[$k] < 0) {
	#			$w = $force_factor*$forcetmp;
	#			for ($ii = -2; $ii <= 2; $ii++) {
	#				for ($jj = -3; $jj <= 2; $jj++) {
	#					$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
	#					$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
	#					$bulk = 1;
	#					&OutString_width($int0[$k], $int1[$k], 0.2, 0.01, $pd_xshift, $pd_zshift);
	#					if ($bulk == 0) {
	#						&OutString_width_nvec($i, $j, $nrvec_x[$k], $nrvec_y[$k], $nrvec_z[$k],	0.2, 0.01, $pd_xshift, $pd_zshift);
	#					}
	#				}
	#			}
	#		}
	#	}
	if (0) {
		printf OUT "y 4\n";

		$cont_bond = 0;
		for ($k = 0; $k < $num_interaction; $k ++) {
			#$forcetmp = $force[$k];
			#$forcetmp = $F_lub[$k];
			if ($contactstate[$k] == 0) {
				# $w = $force_factor*$forcetmp;
				$cont_bond++ ;
				printf OUT "@ 7\n";
			} else {
				printf OUT "@ 5\n";
			}
			for ($ii = $ii_min; $ii <= $ii_max; $ii++) {
				for ($jj = $jj_min; $jj <= $jj_max; $jj++) {
					$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
					$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
					$bulk = 1;
					&OutString_width($int0[$k], $int1[$k], 0.2, 0.01, $pd_xshift, $pd_zshift);
					#						if ($bulk == 0) {
					#	&OutString_width_nvec($i, $j, $nrvec_x[$k], $nrvec_y[$k], $nrvec_z[$k], $Gap[$k], $w, 0.01, $pd_xshift, $pd_zshift);
					#}
				}
			}

		}
		printf "$cont_bond $num_interaction\n";
	}
	
	## visualize force chain network
	if (0) {
		printf OUT "y 4\n";
		printf OUT "@ 7\n";
		for ($k = 0; $k < $num_interaction; $k ++) {
			$forcetmp = $force[$k];
			#$forcetmp = $F_lub[$k];
			if ($forcetmp > 0.1) {
				$w = $force_factor*$forcetmp;
				for ($ii = $ii_min; $ii <= $ii_max; $ii++) {
					for ($jj = $jj_min; $jj <= $jj_max; $jj++) {
						$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
						$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
						$bulk = 1;
						&OutString_width($int0[$k], $int1[$k], $w, 0.01, $pd_xshift, $pd_zshift);
						if ($bulk == 0) {
							&OutString_width_nvec($i, $j, $nrvec_x[$k], $nrvec_y[$k], $nrvec_z[$k], $Gap[$k], $w, 0.01, $pd_xshift, $pd_zshift);
						}
					}
				}
			}
		}
	}
	
	#		printf OUT "y 3\n";
	#		printf OUT "@ 5\n";
	#		for ($k = 0; $k < $num_interaction; $k ++) {
	#			$forcetmp = $F_lub[$k];
	#			if ($forcetmp < 0) {
	#				$w = -$force_factor*$forcetmp;
	#				for ($ii = -2; $ii <= 2; $ii++) {
	#					for ($jj = -3; $jj <= 2; $jj++) {
	#						$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
	#						$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
	#						$bulk = 1;
	#						&OutString_width($int0[$k], $int1[$k], $w, 0.01, $pd_xshift, $pd_zshift);
	#						if ($bulk == 0) {
	#							&OutString_width_nvec($i, $j, $nrvec_x[$k], $nrvec_y[$k], $nrvec_z[$k], $w, 0.01, $pd_xshift, $pd_zshift);
	#						}
	#					}
	#				}
	#			}
	#		}
	
	#		printf OUT "y 5\n";
	#		printf OUT "@ 5\n";
	#
	#		for ($i = 0; $i < $np; $i++) {
	#			$rr = $yap_radius*$radius[$i];
	#			printf OUT "r $rr\n";
	#			#			if ($i == $np-1) {
	#			#				printf OUT "@ 4\n";
	#			#			}
	#
	#			$xi = $posx[$i];
	#			$yi = $posy[$i];
	#			$zi = $posz[$i];
	#			$xj = $posx[$i] + 0.1*$velx[$i];
	#			$yj = $posy[$i];
	#			$zj = $posz[$i] + 0.1*$velz[$i];
	#
	#			printf OUT "l $xi $yi $zi $xj $yj $zj\n";
	#
	#		}
	
	
	
	## visualize rotation in 2D
	
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

sub OutParticleShear {
	for ($i = 0; $i < $np; $i++) {
		$rr = $yap_radius*$radius[$i];
		printf OUT "r $rr\n";
		#			if ($i == $np-1) {
		#				printf OUT "@ 4\n";
		#			}
	LOOP: for ($ii = -1; $ii <= 1; $ii++) {
		for ($jj = -1; $jj <= 1; $jj++) {
			$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
			$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
			$xx =  $posx[i] + $pd_xshift;
			$zz =  $posz[i] + $pd_zshift;
			if ($xx> 0 && $xx < $Lx &&  $zz > 0 && $zz < $Lz){
				printf OUT "c $xx 0 $zz \n";
			}
			last LOOP if ($xx<$Lx && $xx > 0 && $zz > 0 && $zz < $Lz);
		}
	}
	}
}


sub OutParticleExtension {
	$particlecolor = 8;
	printf OUT "@ 8\n";
	$j = 0;
	for ($ii = $ii_min; $ii <= $ii_max; $ii++) {
		for ($jj = $jj_min; $jj <= $jj_max; $jj++) {
			if ($pd_color) {
				if ($particlecolor == 10){
					$particlecolor = 8;
				}
				printf OUT "@ $particlecolor\n";
				$particlecolor ++;
			}
			for ($i = 0; $i < $np; $i++) {
				$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
				$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
				$xx = $posx[$i] + $pd_xshift;
				$yy = $posy[$i] - $Ly/2;
				$zz = $posz[$i] + $pd_zshift;
				if ($axisrotation) {
					$xx2 = cos($magicangle)*$xx - sin($magicangle)*$zz;
					$zz2 = sin($magicangle)*$xx + cos($magicangle)*$zz;
				} else {
					$xx2 = $xx;
					$zz2 = $zz;
					
				}
				#if ($xx2> 0 && $xx2 < $Lx &&  $zz2 > 0 && $zz2 < $Lz){
				if (abs($xx2)<$scale*$Lx && abs($zz2) < $scale*$Lz) {
					#                    if (abs($ii) >= 2 || abs($jj) >= 2) {
					#   printf "i,j = $ii, $jj\n"
					#}
					if ($phi6_data == 1) {
						$op = 119*(1-$phi6abs[$i])+4;
						$color = int $op;
						## printf OUT "@ $color\n";
					} else {
						if ($draw_trajectory) {
							if ($i == $i_trac) {
								printf OUT "@ 5\n";
							} else {
								printf OUT "@ 8\n";
							}
						}
					}
					$rr = $yap_radius*$radius[$i];
					printf OUT "r $rr\n";
					$posx2[$j] = $xx2;
					$posy2[$j] = $posy[$i];
					$posz2[$j] = $zz2;
					$radius2[$j] = $radius[$i];
					$ang2[$j] = $ang[$i];
					$j++;
					printf OUT "c $xx2 $yy $zz2 \n";
				}
				# last LOOP2 if ($xx2<$Lx && $xx2 > 0 && $zz2 > 0 && $zz2 < $Lz);
			}
		}
	}
	if ($draw_cross) {
		printf OUT "y 6\n";
		printf OUT "@ 1\n";
		for ($i = 0; $i < $j; $i++) {
			OutCross($i);
		}
	}
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
		
		#		printf OUT "l 0 0 $Lz    $Lx 0 $Lz\n";
		#printf OUT "l 0 0 0   $Lx 0 0\n";
		#printf OUT "l 0 0 0  0  0 $Lz\n";
		#printf OUT "l $Lx 0 $Lz   $Lx 0 0\n";
		
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
	($i, $j, $w, $delta, $dx, $dz) = @_;
	$xi = $posx[$i] + $dx;
	$yi = $posy[$i] - $delta - $ly2;
	$zi = $posz[$i] + $dz;
	$xj = $posx[$j] + $dx;
	$yj = $posy[$j] - $delta - $ly2;
	$zj = $posz[$j] + $dz;
	
	if ($axisrotation) {
		$xi2 = cos($magicangle)*$xi - sin($magicangle)*$zi;
		$zi2 = sin($magicangle)*$xi + cos($magicangle)*$zi;
		$xj2 = cos($magicangle)*$xj - sin($magicangle)*$zj;
		$zj2 = sin($magicangle)*$xj + cos($magicangle)*$zj;
		$xi = $xi2;
		$zi = $zi2;
		$xj = $xj2;
		$zj = $zj2;
	}
	
	if (abs($xi)<$scale*$Lx &&  abs($zi)<$scale*$Lz){
		$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2 ;
		if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1) {
			$bulk = 1;
			printf OUT "r $w\n";
			printf OUT "s $xi $yi $zi $xj $yj $zj\n";
			#			if ($Gap[$k] < -0.1) {
			#				$xx = 0.5*($xi +$xj);
			#				$zz = 0.5*($zi +$zj);
			#				printf OUT "t $xx 0 $zz $Gap[$k]\n";
			#			}
		} else {
			$bulk = 0;
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

sub OutString_width_nvec {
	($i, $j, $nx, $ny, $nz, $gap, $w, $delta, $dx, $dz) = @_;
	$rr =  $radius[$i] + $radius[$j] + $gap;
	printf OUT "r $w\n";
	$xi = $posx[$i] + $dx;
	$yi = $posy[$i] - $delta- $ly2;
	$zi = $posz[$i] + $dz;
	$xj = $posx[$i] + $rr*$nx + $dx;
	$yj = $posy[$i] - $delta- $ly2;
	$zj = $posz[$i] + $rr*$nz + $dz;
	if ($axisrotation) {
		$xi2 = cos($magicangle)*$xi - sin($magicangle)*$zi;
		$zi2 = sin($magicangle)*$xi + cos($magicangle)*$zi;
		$xj2 = cos($magicangle)*$xj - sin($magicangle)*$zj;
		$zj2 = sin($magicangle)*$xj + cos($magicangle)*$zj;
		$xi = $xi2;
		$zi = $zi2;
		$xj = $xj2;
		$zj = $zj2;
	}
	if (abs($xi)<$scale*$Lx &&  abs($zi)<$scale*$Lz){
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
	}
	$xi = $posx[$j] + $dx;
	$yi = $posy[$j] - $delta- $ly2;
	$zi = $posz[$j] + $dz;
	$xj = $posx[$j] - $rr*$nx + $dx;
	$yj = $posy[$j] - $delta- $ly2;
	$zj = $posz[$j] - $rr*$nz + $dz;
	if ($axisrotation) {
		$xi2 = cos($magicangle)*$xi - sin($magicangle)*$zi;
		$zi2 = sin($magicangle)*$xi + cos($magicangle)*$zi;
		$xj2 = cos($magicangle)*$xj - sin($magicangle)*$zj;
		$zj2 = sin($magicangle)*$xj + cos($magicangle)*$zj;
		$xi = $xi2;
		$zi = $zi2;
		$xj = $xj2;
		$zj = $zj2;
	}
	if (abs($xi)<$scale*$Lx &&  abs($zi)<$scale*$Lz){
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
	$a = $radius2[$i];
	$xi = $posx2[$i];
	$yi = $posy2[$i]+$offset1;
	$zi = $posz2[$i];
	$angle = $ang2[$i];
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

sub yapColor {
	printf OUT "\@0 0 0 0 \n";
	#printf OUT "\@1 50 100 205 \n";
	printf OUT "\@1 25 50 102 \n";
	#printf OUT "\@1 255 255 255 \n";
	printf OUT "\@2 200 200 200 \n";
	printf OUT "\@3 50 150 255 \n";
	printf OUT "\@4 50 200 50 \n";
	printf OUT "\@5 255 100 100 \n";
	printf OUT "\@6 50 200 50 \n";
	printf OUT "\@7 255 255 0 \n";
	printf OUT "\@8 255 255 255\n";
	printf OUT "\@9 150 150 150\n";
	printf OUT "\@10 224 143 0 \n";
	printf OUT "\@11 67 50 230 \n";
	printf OUT "\@12 253 105 6 \n";
	printf OUT "\@13 109 109 109 \n";
}

#
#sub yapColor {
#	printf OUT "@ 0 100 100 100\n";
#	printf OUT "@ 1 0 0 0\n";
#	printf OUT "@ 2 255 255 255\n";
#	printf OUT "@ 125 50 50 50\n";
#	printf OUT "@ 126 255 255 0\n";
#	printf OUT"@ 3 120 27 134\n";
#	printf OUT"@ 4 114 28 139\n";
#	printf OUT"@ 5 109 28 143\n";
#	printf OUT"@ 6 103 28 148\n";
#	printf OUT"@ 7 98 28 153\n";
#	printf OUT"@ 8 92 29 157\n";
#	printf OUT"@ 9 87 29 162\n";
#	printf OUT"@ 10 82 29 167\n";
#	printf OUT"@ 11 78 31 171\n";
#	printf OUT"@ 12 76 35 174\n";
#	printf OUT"@ 13 74 39 178\n";
#	printf OUT"@ 14 72 42 181\n";
#	printf OUT"@ 15 70 46 185\n";
#	printf OUT"@ 16 68 50 189\n";
#	printf OUT"@ 17 65 53 192\n";
#	printf OUT"@ 18 63 57 196\n";
#	printf OUT"@ 19 63 62 197\n";
#	printf OUT"@ 20 63 66 199\n";
#	printf OUT"@ 21 63 71 200\n";
#	printf OUT"@ 22 63 75 202\n";
#	printf OUT"@ 23 62 80 204\n";
#	printf OUT"@ 24 62 85 205\n";
#	printf OUT"@ 25 62 89 207\n";
#	printf OUT"@ 26 62 94 207\n";
#	printf OUT"@ 27 63 98 207\n";
#	printf OUT"@ 28 64 102 206\n";
#	printf OUT"@ 29 64 107 206\n";
#	printf OUT"@ 30 65 111 206\n";
#	printf OUT"@ 31 66 115 205\n";
#	printf OUT"@ 32 67 119 205\n";
#	printf OUT"@ 33 67 124 204\n";
#	printf OUT"@ 34 69 127 202\n";
#	printf OUT"@ 35 70 130 200\n";
#	printf OUT"@ 36 71 134 198\n";
#	printf OUT"@ 37 73 137 196\n";
#	printf OUT"@ 38 74 140 193\n";
#	printf OUT"@ 39 75 144 191\n";
#	printf OUT"@ 40 77 147 189\n";
#	printf OUT"@ 41 78 150 187\n";
#	printf OUT"@ 42 80 152 183\n";
#	printf OUT"@ 43 82 155 180\n";
#	printf OUT"@ 44 84 157 177\n";
#	printf OUT"@ 45 86 160 174\n";
#	printf OUT"@ 46 88 162 171\n";
#	printf OUT"@ 47 90 164 167\n";
#	printf OUT"@ 48 91 167 164\n";
#	printf OUT"@ 49 94 168 161\n";
#	printf OUT"@ 50 96 170 157\n";
#	printf OUT"@ 51 98 171 153\n";
#	printf OUT"@ 52 101 173 150\n";
#	printf OUT"@ 53 103 175 146\n";
#	printf OUT"@ 54 106 176 143\n";
#	printf OUT"@ 55 108 178 139\n";
#	printf OUT"@ 56 111 179 136\n";
#	printf OUT"@ 57 113 180 132\n";
#	printf OUT"@ 58 116 181 129\n";
#	printf OUT"@ 59 119 182 125\n";
#	printf OUT"@ 60 122 183 122\n";
#	printf OUT"@ 61 125 184 119\n";
#	printf OUT"@ 62 128 185 115\n";
#	printf OUT"@ 63 130 186 112\n";
#	printf OUT"@ 64 134 186 109\n";
#	printf OUT"@ 65 137 187 106\n";
#	printf OUT"@ 66 140 187 104\n";
#	printf OUT"@ 67 143 188 101\n";
#	printf OUT"@ 68 146 188 98\n";
#	printf OUT"@ 69 150 188 95\n";
#	printf OUT"@ 70 153 189 92\n";
#	printf OUT"@ 71 156 189 90\n";
#	printf OUT"@ 72 159 189 88\n";
#	printf OUT"@ 73 163 189 86\n";
#	printf OUT"@ 74 166 189 84\n";
#	printf OUT"@ 75 169 189 82\n";
#	printf OUT"@ 76 173 189 80\n";
#	printf OUT"@ 77 176 189 78\n";
#	printf OUT"@ 78 179 189 76\n";
#	printf OUT"@ 79 182 188 74\n";
#	printf OUT"@ 80 185 187 73\n";
#	printf OUT"@ 81 188 187 72\n";
#	printf OUT"@ 82 191 186 71\n";
#	printf OUT"@ 83 195 185 69\n";
#	printf OUT"@ 84 198 184 68\n";
#	printf OUT"@ 85 201 184 67\n";
#	printf OUT"@ 86 203 183 66\n";
#	printf OUT"@ 87 206 181 65\n";
#	printf OUT"@ 88 208 179 64\n";
#	printf OUT"@ 89 210 177 63\n";
#	printf OUT"@ 90 213 176 62\n";
#	printf OUT"@ 91 215 174 61\n";
#	printf OUT"@ 92 217 172 61\n";
#	printf OUT"@ 93 220 171 60\n";
#	printf OUT"@ 94 221 168 59\n";
#	printf OUT"@ 95 222 165 58\n";
#	printf OUT"@ 96 224 162 58\n";
#	printf OUT"@ 97 225 159 57\n";
#	printf OUT"@ 98 226 156 56\n";
#	printf OUT"@ 99 227 153 56\n";
#	printf OUT"@ 100 229 150 55\n";
#	printf OUT"@ 101 229 146 54\n";
#	printf OUT"@ 102 229 142 53\n";
#	printf OUT"@ 103 229 137 53\n";
#	printf OUT"@ 104 229 133 52\n";
#	printf OUT"@ 105 230 128 51\n";
#	printf OUT"@ 106 230 124 50\n";
#	printf OUT"@ 107 230 120 49\n";
#	printf OUT"@ 108 230 115 48\n";
#	printf OUT"@ 109 229 110 47\n";
#	printf OUT"@ 110 228 104 46\n";
#	printf OUT"@ 111 227 99 45\n";
#	printf OUT"@ 112 226 93 44\n";
#	printf OUT"@ 113 226 88 43\n";
#	printf OUT"@ 114 225 82 42\n";
#	printf OUT"@ 115 224 77 41\n";
#	printf OUT"@ 116 223 72 40\n";
#	printf OUT"@ 117 222 66 39\n";
#	printf OUT"@ 118 222 60 38\n";
#	printf OUT"@ 119 221 55 37\n";
#	printf OUT"@ 120 220 49 36\n";
#	printf OUT"@ 121 220 44 35\n";
#	printf OUT"@ 122 219 38 34\n";
#
#}
