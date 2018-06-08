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

my $output_interval = 1;
my $xz_shift = 0;
my $axis = 0;
my $reversibility_test = 0;
my $monodisperse = 0;
my $twodim = 0;
my $pi = atan2(1, 1) * 4;
my $shear_strain_min = 1;
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
printf "$particle_data\n";
$i = index($particle_data, "par_", 0)+4;
$iD3 = index($particle_data, "D3N");
$j = index($particle_data, ".dat", $i-1);
$name = substr($particle_data, $i, $j-$i);
my $dim = 3;
if ($iD3 eq -1){
    $dim = 2;
}
printf "$dim dimension\n";;

$interaction_data = "int_${name}.dat";

$outputDataReconstruct = "r_$name.dat";

open (OUTDR, "> ${outputDataReconstruct}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");

&readHeader;
&yaplotColor;

$cnt_interval = 0;
$first = 1;
$first_int = 1;
$checkpoint = 1;
$shear_strain_previous = 0;
$shearrate_positive = 1;

$jbin = 100;
$delta_angle = $pi/$jbin;

for ($j = 0; $j < $jbin; $j ++){
    $n1bin[$j] = 0;
}

for ($j = 0; $j < $jbin; $j ++){
    $etabin[$j] = 0;
}


$n1_count = 0;
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
	&InInteractions;
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

$nangle = $pi + atan2($ny, $nx);
$iangle = int($nangle/$delta_angle);

$jbin = 100;
for ($j = 0; $j < $jbin; $j ++){
    $angle = $delta_angle*$j + $delta_angle/2;
    $angle_diff = $angle - 3*$pi/4;
    if ($angle_diff < -$pi/2) {
        $angle_diff += $pi;
    }
    $n1 = $jbin*$n1bin[$j]/($n1_count);
    $eta = 6*$pi*$jbin*$etabin[$j]/($n1_count);
    $force = $forcebin[$j]/($n1_count);
    printf "$angle_diff $n1 $eta $force\n";
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
	$line = <IN_particle>; ($buf, $buf, $flwtyp) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $dataunit) = split(/\s+/, $line);
	
	if ($Ly==0) {
		$number_of_header = 8;
	} else {
		$number_of_header = 7;
	}
	for ($i = 0; $i<$number_of_header; $i++) {
		$line = <IN_particle>;
		printf "$line";
	}
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
    #	$target_stress = $ssHeader[3];
    $viscosity = $ssHeader[3];
    $normalstressdiff1 = $ssHeader[4];
    
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
                ($ip, $a, $x, $y, $z, $vx, $vz, $vy, $ox, $oz, $oy) = split(/\s+/, $line);
            } else {
                ($ip, $a, $x, $z, $vx, $vz, $vy, $ox, $oz, $oy, $angle) = split(/\s+/, $line);
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
	}

	
    if ($shear_strain > $shear_strain_min) {
        $n1_count ++;
    }
	## visualize force chain network
    printf OUT "y 6\n";
    printf OUT "@ 7\n";
    $etaTotal = 0;
    
    $eta_denominator = ($Lx*$Lz);
    for ($k = 0; $k < $num_interaction; $k ++) {
        $force = $F_lub[$k] + $Fc_n[$k];
        $nx = $nrvec_x[$k];
        $ny = $nrvec_y[$k];
        $eta_lub = -$F_lub[$k]*$distance[$k]*($nx*$ny);
        $eta_con = -$Fc_n[$k]*$distance[$k]*($nx*$ny);
        $eta = $eta_lub + $eta_con;
        $etaTotal += $eta;
        
        if ($shear_strain > $shear_strain_min) {
            $nangle = $pi + atan2($ny, $nx);
            $iangle = int($nangle/$delta_angle);
            #        printf "$iangle\n";
            $etabin[$iangle] += $eta/$eta_denominator;
            $forcebin[$iangle] += $force;
        }
    }
    $visapprox = $etaTotal/($Lx*$Lz);
    
    $n1approx_pos_con = 0;
    $n1approx_neg_con = 0;
    $n1approx_pos_lub = 0;
    $n1approx_neg_lub = 0;

    $force_factor = 0.015/$viscosity;
    printf OUT "y 5\n";
    #    printf OUT "@ 7\n";
    $n1_denominator = ($Lx*$Lz*$viscosity);
    for ($k = 0; $k < $num_interaction; $k ++) {
        $nx = $nrvec_x[$k];
        $ny = $nrvec_y[$k];
        $n1lub = -$F_lub[$k]*$distance[$k]*($nx*$nx - $ny*$ny);
        $n1con =  -$Fc_n[$k]*$distance[$k]*($nx*$nx - $ny*$ny);

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
        $n1 = $n1lub + $n1con;
        if ($shear_strain > $shear_strain_min) {
            $nangle = $pi + atan2($ny, $nx);
            $iangle = int($nangle/$delta_angle);
            $n1bin[$iangle] += $n1/$n1_denominator;
        }
    }
    $barsize = 10;
    $n1approxP = $n1approx_pos_con/($Lx*$Lz*$viscosity);
    $n1approxN = $n1approx_neg_con/($Lx*$Lz*$viscosity);
    $n1approxCon = $n1approxP + $n1approxN;

    $testposition1 = -$Lx/2 - 3.2;
    $testposition2 = -$Lx/2 - 5.2;

    $n1approxP = $n1approx_pos_lub/($Lx*$Lz*$viscosity);
    $n1approxN = $n1approx_neg_lub/($Lx*$Lz*$viscosity);

    $n1approxLub = $n1approxP + $n1approxN;
    $testposition1 = -$Lx/2 - 1;
    $testposition2 = -$Lx/2 - 3;

    
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


