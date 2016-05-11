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
my $force_factor = 0.001;
my $output_interval = 1;
my $xz_shift = 0;
my $axis = 0;
my $reversibility_test = 0;
my $monodisperse = 0;
my $rotatingobserver = 0;
#my $np_movable = 500;
#my $rout = 37.5726;
#my $rout = 31.8812;
#my $rout = 40.4141;

my $rin = $rout/2;
my $calcrheology = 0;

if ($calcrheology == 1) {
	$rotatingobserver = 0;
}

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

$j = index($name, '_fric', 1);## not flexible...
$initconfig = substr($name, 0, $j);
printf "$initconfig\n";

open (IN_CONFIG, "< ${initconfig}.dat");

$line = <IN_CONFIG>;
$line = <IN_CONFIG>;
($buf, $np1, $np2, $vf, $lx, $ly, $lz, $np_in, $np_out, $rin, $rout) = split(/\s+/, $line);

my $np_movable = $np1+$np2;

$interaction_data = "int_${name}.dat";

if ($calcrheology == 0) {
	$output = "y_$name.yap";
	open(OUT, "> ${output}");
} else {
	$output = "CGrheo_$name.dat";
	open(OUT, "> ${output}");
}
open(IN_particle, "< ${particle_data}");
open(IN_interaction, "< ${interaction_data}");


&readHeader;
if ($calcrheology == 0) {
	#### &yaplotColor;
	&yaplotColorGradient;
}

$cnt_interval = 0;
$first = 1;
$checkpoint = 1;
$shear_strain_previous = 0;
$shearrate_positive = 1;

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
	$force_in = 0;
	$force_out = 0;
	# &InInteractions;
	if ($reversibility_test) {
		if ($first || $checkpoint == 1) {
			&keepInitialConfig;
		}
	}
	if ($calcrheology == 0) {
		if ($output == 1) {
			&OutYaplotData;
		}
	} else {
		printf "$shear_strain_i $shear_strain\n";
		$viscosity_in =  3*$force_in/(2*$rin);
		$viscosity = 3*$force_out/(2*$rout);
		printf OUT "$shear_strain $shear_rate $viscosity $viscosity_in\n";
	}
	$cnt_interval ++;
}

close (OUT);

close (IN_particle);
close (IN_interaction);

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
	for ($i = 0; $i<11; $i++) {
		$line = <IN_particle>;
	}	
	for ($i = 0; $i<24; $i++) {
		$line = <IN_interaction>;
	}
}

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
	if (defined $line) {
		# 1 sys.get_shear_strain()
		# 2 sys.shear_disp
		# 3 getRate()
		# 4 target_stress_input
		# 5 sys.get_time()
		# 6 sys.angle_external_magnetic_field
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
		$max_stress = 0;
		$min_stress = 0;
		for ($i = 0; $i < $np; $i ++){
			$line = <IN_particle>;
			if ($output == 1) {
				# 1: number of the particle
				# 2: radius
				# 3, 4, 5: position
				# 6, 7, 8: velocity
				# 9, 10, 11: angular velocity
				# 12: viscosity contribution of lubrication
				# 13: viscosity contributon of contact GU xz
				# 14: viscosity contributon of brownian xz
				# (15: angle for 2D simulation)
				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz, $s_rr, $s_tt, $s_rt, $angle) = split(/\s+/, $line);
				#
				$ang[$i] = $angle;
				$radius[$i] = $a;
				if ($xz_shift) {
					$x += $xz_shift;
					$z += $xz_shift;
					if ($x > $Lx/2) {
						$x -= $Lx;
					}
					if ($z > $Lz/2) {
						$z -= $Lz;
					}
				}
				if ($rotatingobserver) {
					$posxtmp[$i] = $x;
					$posy[$i] = $y;
					$posztmp[$i] = $z;
				} else {
					$posx[$i] = $x;
					$posy[$i] = $y;
					$posz[$i] = $z;
				}
				$velx[$i] = $vx;
				$vely[$i] = $vy;
				$velz[$i] = $vz;
				$omegax[$i] = $ox;
				$omegay[$i] = $oy;
				$omegaz[$i] = $oz;
				$omegay[$i] = $oy;
				$stress_output[$i] = $s_tt + $s_rr; # particle pressure
				if ($max_stress < $s_tt) {
					$max_stress = $s_tt;
				}
				if ($min_stress > $s_tt) {
					$min_stress = $s_tt;
				}
				if ($radius_max < $a) {
					$radius_max = $a;
				}
			}
		}
		printf "$min_stress $max_stress \n";
		if ($rotatingobserver) {
			if ($first) {
				$rout = sqrt($posxtmp[$np-1]*$posxtmp[$np-1] + $posztmp[$np-1]*$posztmp[$np-1]);
				#printf "$rout \n";
			}
			$theta = -0.25*$shear_strain;
			#printf "$theta, $shear_strain, $rout \n";
			for ($i = 0; $i < $np; $i ++){
				$posx[$i] = $posxtmp[$i]*cos($theta)-$posztmp[$i]*sin($theta);
				$posz[$i] = $posxtmp[$i]*sin($theta)+$posztmp[$i]*cos($theta);
			}
		}
	}
}

sub InInteractions {
	$line = <IN_interaction>;
	($buf, $shear_strain_i, $num_interaction) = split(/\s+/, $line);
	if ($buf != "#") {
		exit(1);
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
	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction>;
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
		#if ($output == 1) {
		($i, $j, $contact, $nx, $ny, $nz, #1---6
		$gap, $f_lub_norm, # 7, 8
		$f_lub_tan_x, $f_lub_tan_y, $f_lub_tan_z, # 9, 10, 11
		$fc_norm, # 12
		$fc_tan_x, $fc_tan_y, $fc_tan_z, # 13, 14, 15
		$fr_norm, $s_xF) = split(/\s+/, $line);
		$int0[$k] = $i;
		$int1[$k] = $j;
		$contactstate[$k] = $contact;
        $force[$k] = $fc_norm + $f_lub_norm + $fr_norm;
        #$force[$k] = $f_lub_norm + $fr_norm;
		$F_lub[$k] = $f_lub_norm;
		$Fc_n[$k] = $fc_norm;
		$Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
		$S_bf[$k] = $s_xF;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		$Gap[$k] = $gap;
		if ($i < $np_movable && $j >= $np_movable) {
			#			$fx = $f_lub_tan_x + $fc_tan_x + $force[$k]*$nx;
			#$fz = $f_lub_tan_z + $fc_tan_z + $force[$k]*$nz;
			$fx = $force[$k]*$nx;
			$fz = $force[$k]*$nz;
			$posr = sqrt($posx[$j]*$posx[$j] + $posz[$j]*$posz[$j]);
			$nrx = $posx[$j]/$posr;
			$nrz = $posz[$j]/$posr;
			$f_tan = -$fx*$nrz + $fz*$nrx;
			if ($posr/$rout > 0.9) {
				$force_out += $f_tan;
			} else {
				$force_in += $f_tan;
			}
		}
	}
}

sub OutYaplotData{
	if ($first == 0) {
		printf OUT "\n";
	} else {
		$first = 0;
		$theta0 = atan2($posz[$np_movable], $posx[$np_movable]);
		$rwheel = sqrt($posz[$np_movable]*$posz[$np_movable] + $posx[$np_movable]*$posx[$np_movable]);
	}
	printf OUT "y 1\n";
	printf OUT "@ 8\n";
	$r = $yap_radius*$radius[0];
	printf OUT "r $r\n";
	$switch = 0;
	## visualize particles
	if ($monodisperse) {
		printf OUT "r $radius[0]\n";
		for ($i = 0; $i < $np; $i++) {
			printf OUT "c $posx[$i] $posy[$i] $posz[$i]\n";
		}
	} else {
		printf "@ 8\n";
		for ($i = 0; $i < $np_movable; $i++) {
			printf OUT "r $radius[$i]\n";
			$normalized_stress = -$stress_output[$i]/40000;
			if ($normalized_stress > 1) {
				$normalized_stress = 1;
			} elsif ($normalized_stress < 0) {
				$normalized_stress = 0;
			}
			#			if ($normalized_stress > 0) {
			$istress =  6 + int(120*$normalized_stress);
			if ($istress > 125) {
				$istress = 125;
			}
			$color = $istress;
			printf OUT "@ $color \n";
			printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
		}
		printf OUT "y 2\n";
		printf OUT "@ 2\n";
		printf OUT "r $radius[$np_movable]\n";
		for ($i = $np_movable; $i < $np; $i++) {
			printf OUT "c $posx[$i] $posy[$i] $posz[$i]\n";
		}
		printf OUT "@ 2\n";
		printf OUT "r 0.5\n";
		$xo = $Lx/2;
		$zo = $Lz/2;
		$np_in_end = $np_in + $np_movable;
		$theta = atan2($posz[$np_movable], $posx[$np_movable]) - $theta0;
		for ($j = 0; $j < 4; $j++) {
			$xcross = $rwheel*cos($theta+$j*pi/2);
			$zcross = $rwheel*sin($theta+$j*pi/2);
			printf OUT "s 0 -0.01 0 $xcross -0.01 $zcross\n";
		}
	}
	## visualize contact network
	#	printf OUT "y 2\n";
	#	printf OUT "r 0.2\n";
	#	printf OUT "@ 2\n"; # static
	#	for ($k = 0; $k < $num_interaction; $k ++) {
	#		if ($contactstate[$k] == 2) {
	#			&OutString2($int0[$k], $int1[$k]);
	#		}
	#	}
	## visualize force chain network
#	printf OUT "y 4\n";
#	printf OUT "@ 7\n";
#	for ($k = 0; $k < $num_interaction; $k ++) {
#		#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
#		if ($int0[$k] < $np_movable || $int1[$k] < $np_movable) {
#			if ($force[$k] > 0) {
#				&OutString_width($int0[$k], $int1[$k], $force_factor*$force[$k], 0.01);
#			}
#		}
#	}
#	printf OUT "y 3\n";
#	printf OUT "@ 5\n";
#	for ($k = 0; $k < $num_interaction; $k ++) {
#		#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
#		if ($force[$k] < 0) {
#			&OutString_width($int0[$k], $int1[$k], -$force_factor*$force[$k], -0.02);
#		}
#	}
	#
	## visualize rotation in 2D
#	printf OUT "y 6\n";
#	printf OUT "@ 0\n";
#	for ($i = 0; $i < $np_movable; $i++) {
#		OutCross($i);
#	}
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
				printf OUT "c $posx_init[$i] 0.01 $posz_init[$i]\n";
			} else {
				printf OUT "c $posx_init[$i] $posy_init[$i] $posz_init[$i]\n";
			}
		}
	}
	#	&OutBoundaryBox;
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
		#printf OUT "p 4 -$lx2 $yb $lz2 $lx2 $yb $lz2 $lx2 $yb -$lz2 -$lx2 $yb -$lz2 \n";
		printf OUT "l -$lx2 0  $lz2  $lx2 0  $lz2\n";
		printf OUT "l -$lx2 0 -$lz2  $lx2 0 -$lz2\n";
		printf OUT "l -$lx2 0 -$lz2 -$lx2 0  $lz2\n";
		printf OUT "l  $lx2 0  $lz2  $lx2 0 -$lz2\n";
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
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1) {
		printf OUT "r $w\n";
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
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
	$offset1 = 0.01;
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



sub yaplotColor {
	printf OUT "\@0 0 0 0\n";
	#printf OUT "\@1 50 100 205\n";
	printf OUT "\@1 25 50 102\n";
	#printf OUT "\@1 255 255 255\n";
	printf OUT "\@2 200 200 200\n";
	printf OUT "\@3 50 150 255\n";
	printf OUT "\@4 50 200 50\n";
	printf OUT "\@5 255 100 100\n";
	printf OUT "\@6 50 200 50\n";
	printf OUT "\@7 255 255 0\n";
	printf OUT "\@8 255 255 255\n";
	printf OUT "\@9 150 150 150\n";
	#printf OUT "\@8 224 143 0\n";
	#printf OUT "\@9 67 163 230\n";
	#printf OUT "\@8 253 105 6\n";
	#printf OUT "\@9 109 109 109\n";
	printf OUT "\@10 250 250 250\n";
	printf OUT "\@11 240 240 240\n";
	printf OUT "\@12 230 230 230\n";
	printf OUT "\@13 220 220 220\n";
	printf OUT "\@14 210 210 210\n";
	printf OUT "\@15 200 200 200\n";
	printf OUT "\@16 190 190 190\n";
	printf OUT "\@17 180 180 180\n";
	printf OUT "\@18 170 170 170\n";
	printf OUT "\@19 160 160 160\n";
	printf OUT "\@20 150 150 150\n";
	printf OUT "\@21 140 140 140\n";
	printf OUT "\@22 130 130 130\n";
	printf OUT "\@23 120 120 120\n";
	printf OUT "\@24 110 110 110\n";
	printf OUT "\@25 100 100 100\n";
	printf OUT "\@26 90 90 90\n";
	printf OUT "\@27 80 90 90\n";
	printf OUT "\@28 70 70 70\n";
	printf OUT "\@29 60 60 60\n";
	printf OUT "\@30 50 50 50\n";
	printf OUT "\@31 40 40 40\n";
	printf OUT "\@32 30 30 30\n";
	printf OUT "\@33 20 20 20\n";
	printf OUT "\@34 10 10 10\n";
	printf OUT "\@35 0 0 0\n";
}

sub yaplotColorGradient_1 {
	printf OUT "\@0 0 0 0\n";
	#printf OUT "\@1 50 100 205\n";
	printf OUT "\@1 25 50 102\n";
	#printf OUT "\@1 255 255 255\n";
	printf OUT "\@2 200 200 200\n";
	printf OUT "\@6 42 71 238\n";
	printf OUT "\@7 45 77 238\n";
	printf OUT "\@8 49 83 239\n";
	printf OUT "\@9 53 88 239\n";
	printf OUT "\@10 56 94 239\n";
	printf OUT "\@11 60 99 239\n";
	printf OUT "\@12 64 105 240\n";
	printf OUT "\@13 67 111 240\n";
	printf OUT "\@14 71 116 240\n";
	printf OUT "\@15 75 122 240\n";
	printf OUT "\@16 79 127 241\n";
	printf OUT "\@17 82 133 241\n";
	printf OUT "\@18 86 138 241\n";
	printf OUT "\@19 90 144 241\n";
	printf OUT "\@20 93 149 242\n";
	printf OUT "\@21 97 153 242\n";
	printf OUT "\@22 100 157 242\n";
	printf OUT "\@23 104 161 243\n";
	printf OUT "\@24 108 165 243\n";
	printf OUT "\@25 111 169 243\n";
	printf OUT "\@26 115 173 244\n";
	printf OUT "\@27 118 177 244\n";
	printf OUT "\@28 122 181 244\n";
	printf OUT "\@29 126 185 245\n";
	printf OUT "\@30 129 189 245\n";
	printf OUT "\@31 133 193 245\n";
	printf OUT "\@32 136 197 246\n";
	printf OUT "\@33 140 200 246\n";
	printf OUT "\@34 143 203 246\n";
	printf OUT "\@35 146 206 247\n";
	printf OUT "\@36 149 208 247\n";
	printf OUT "\@37 153 211 247\n";
	printf OUT "\@38 156 214 248\n";
	printf OUT "\@39 159 216 248\n";
	printf OUT "\@40 162 219 248\n";
	printf OUT "\@41 166 222 248\n";
	printf OUT "\@42 169 225 249\n";
	printf OUT "\@43 172 227 249\n";
	printf OUT "\@44 175 230 249\n";
	printf OUT "\@45 178 233 250\n";
	printf OUT "\@46 182 235 250\n";
	printf OUT "\@47 184 237 250\n";
	printf OUT "\@48 187 238 250\n";
	printf OUT "\@49 190 239 250\n";
	printf OUT "\@50 193 241 251\n";
	printf OUT "\@51 196 242 251\n";
	printf OUT "\@52 199 243 251\n";
	printf OUT "\@53 202 245 251\n";
	printf OUT "\@54 205 246 251\n";
	printf OUT "\@55 208 247 252\n";
	printf OUT "\@56 211 249 252\n";
	printf OUT "\@57 214 250 252\n";
	printf OUT "\@58 216 252 252\n";
	printf OUT "\@59 219 253 252\n";
	printf OUT "\@60 222 253 250\n";
	printf OUT "\@61 225 253 247\n";
	printf OUT "\@62 227 253 244\n";
	printf OUT "\@63 230 253 241\n";
	printf OUT "\@64 232 253 238\n";
	printf OUT "\@65 235 253 234\n";
	printf OUT "\@66 237 253 231\n";
	printf OUT "\@67 240 253 228\n";
	printf OUT "\@68 242 253 225\n";
	printf OUT "\@69 245 253 222\n";
	printf OUT "\@70 247 253 219\n";
	printf OUT "\@71 250 253 215\n";
	printf OUT "\@72 252 253 212\n";
	printf OUT "\@73 254 253 209\n";
	printf OUT "\@74 253 252 206\n";
	printf OUT "\@75 253 251 202\n";
	printf OUT "\@76 253 250 199\n";
	printf OUT "\@77 253 249 196\n";
	printf OUT "\@78 252 248 193\n";
	printf OUT "\@79 252 247 189\n";
	printf OUT "\@80 252 246 186\n";
	printf OUT "\@81 251 245 183\n";
	printf OUT "\@82 251 244 180\n";
	printf OUT "\@83 251 243 176\n";
	printf OUT "\@84 250 243 173\n";
	printf OUT "\@85 250 242 170\n";
	printf OUT "\@86 250 241 167\n";
	printf OUT "\@87 249 238 163\n";
	printf OUT "\@88 248 236 160\n";
	printf OUT "\@89 248 233 156\n";
	printf OUT "\@90 247 231 153\n";
	printf OUT "\@91 246 229 150\n";
	printf OUT "\@92 245 226 146\n";
	printf OUT "\@93 245 224 143\n";
	printf OUT "\@94 244 221 140\n";
	printf OUT "\@95 243 219 136\n";
	printf OUT "\@96 242 216 133\n";
	printf OUT "\@97 242 214 130\n";
	printf OUT "\@98 241 212 126\n";
	printf OUT "\@99 240 209 123\n";
	printf OUT "\@100 239 206 120\n";
	printf OUT "\@101 238 202 116\n";
	printf OUT "\@102 238 198 113\n";
	printf OUT "\@103 237 195 110\n";
	printf OUT "\@104 236 191 107\n";
	printf OUT "\@105 235 187 103\n";
	printf OUT "\@106 234 183 100\n";
	printf OUT "\@107 233 180 97\n";
	printf OUT "\@108 232 176 94\n";
	printf OUT "\@109 231 172 91\n";
	printf OUT "\@110 230 168 87\n";
	printf OUT "\@111 229 165 84\n";
	printf OUT "\@112 228 161 81\n";
	printf OUT "\@113 227 157 78\n";
	printf OUT "\@114 226 152 76\n";
	printf OUT "\@115 225 147 73\n";
	printf OUT "\@116 225 142 71\n";
	printf OUT "\@117 224 137 69\n";
	printf OUT "\@118 223 132 67\n";
	printf OUT "\@119 222 128 64\n";
	printf OUT "\@120 221 123 62\n";
	printf OUT "\@121 220 118 60\n";
	printf OUT "\@122 219 113 57\n";
	printf OUT "\@123 218 108 55\n";
	printf OUT "\@124 217 103 53\n";
	printf OUT "\@125 216 99 51\n";
}

sub yaplotColorGradient {
	
	printf OUT "\@0 25 50 102\n";
	printf OUT "\@1 0 0 0\n";
	#printf OUT "\@1 50 100 205\n";

	#printf OUT "\@1 255 255 255\n";
	printf OUT "\@2 200 200 200\n";
	
	printf OUT "\@ 6 22 37 51\n";
	printf OUT "\@ 7 23 39 53\n";
	printf OUT "\@ 8 25 41 56\n";
	printf OUT "\@ 9 26 44 58\n";
	printf OUT "\@ 10 28 46 60\n";
	printf OUT "\@ 11 29 48 62\n";
	printf OUT "\@ 12 31 51 64\n";
	printf OUT "\@ 13 33 53 66\n";
	printf OUT "\@ 14 34 56 68\n";
	printf OUT "\@ 15 36 58 70\n";
	printf OUT "\@ 16 37 60 72\n";
	printf OUT "\@ 17 39 63 74\n";
	printf OUT "\@ 18 41 65 76\n";
	printf OUT "\@ 19 42 67 79\n";
	printf OUT "\@ 20 44 70 81\n";
	printf OUT "\@ 21 45 72 83\n";
	printf OUT "\@ 22 47 74 85\n";
	printf OUT "\@ 23 48 76 87\n";
	printf OUT "\@ 24 49 78 89\n";
	printf OUT "\@ 25 51 80 91\n";
	printf OUT "\@ 26 52 83 93\n";
	printf OUT "\@ 27 54 85 95\n";
	printf OUT "\@ 28 55 87 97\n";
	printf OUT "\@ 29 56 89 99\n";
	printf OUT "\@ 30 58 91 101\n";
	printf OUT "\@ 31 59 93 103\n";
	printf OUT "\@ 32 60 95 105\n";
	printf OUT "\@ 33 62 97 107\n";
	printf OUT "\@ 34 63 99 109\n";
	printf OUT "\@ 35 65 101 111\n";
	printf OUT "\@ 36 66 103 113\n";
	printf OUT "\@ 37 68 106 115\n";
	printf OUT "\@ 38 71 108 117\n";
	printf OUT "\@ 39 73 111 119\n";
	printf OUT "\@ 40 76 113 121\n";
	printf OUT "\@ 41 78 116 123\n";
	printf OUT "\@ 42 81 119 125\n";
	printf OUT "\@ 43 83 121 127\n";
	printf OUT "\@ 44 86 124 129\n";
	printf OUT "\@ 45 88 126 131\n";
	printf OUT "\@ 46 91 129 133\n";
	printf OUT "\@ 47 93 131 135\n";
	printf OUT "\@ 48 96 134 137\n";
	printf OUT "\@ 49 98 136 139\n";
	printf OUT "\@ 50 101 139 141\n";
	printf OUT "\@ 51 103 142 143\n";
	printf OUT "\@ 52 106 144 144\n";
	printf OUT "\@ 53 109 146 146\n";
	printf OUT "\@ 54 112 149 147\n";
	printf OUT "\@ 55 115 151 148\n";
	printf OUT "\@ 56 118 153 149\n";
	printf OUT "\@ 57 120 156 150\n";
	printf OUT "\@ 58 123 158 152\n";
	printf OUT "\@ 59 126 160 153\n";
	printf OUT "\@ 60 129 163 154\n";
	printf OUT "\@ 61 132 165 155\n";
	printf OUT "\@ 62 135 168 157\n";
	printf OUT "\@ 63 138 170 158\n";
	printf OUT "\@ 64 141 172 159\n";
	printf OUT "\@ 65 143 175 160\n";
	printf OUT "\@ 66 146 177 162\n";
	printf OUT "\@ 67 149 179 162\n";
	printf OUT "\@ 68 152 181 162\n";
	printf OUT "\@ 69 155 183 162\n";
	printf OUT "\@ 70 158 185 163\n";
	printf OUT "\@ 71 161 187 163\n";
	printf OUT "\@ 72 164 189 163\n";
	printf OUT "\@ 73 167 191 163\n";
	printf OUT "\@ 74 169 193 164\n";
	printf OUT "\@ 75 172 195 164\n";
	printf OUT "\@ 76 175 197 164\n";
	printf OUT "\@ 77 178 199 164\n";
	printf OUT "\@ 78 181 201 165\n";
	printf OUT "\@ 79 184 202 165\n";
	printf OUT "\@ 80 187 204 165\n";
	printf OUT "\@ 81 190 206 165\n";
	printf OUT "\@ 82 192 207 165\n";
	printf OUT "\@ 83 194 208 164\n";
	printf OUT "\@ 84 196 209 163\n";
	printf OUT "\@ 85 198 210 162\n";
	printf OUT "\@ 86 200 211 161\n";
	printf OUT "\@ 87 202 212 160\n";
	printf OUT "\@ 88 204 212 159\n";
	printf OUT "\@ 89 206 213 158\n";
	printf OUT "\@ 90 208 214 157\n";
	printf OUT "\@ 91 210 215 156\n";
	printf OUT "\@ 92 212 216 156\n";
	printf OUT "\@ 93 214 217 155\n";
	printf OUT "\@ 94 216 218 154\n";
	printf OUT "\@ 95 218 219 153\n";
	printf OUT "\@ 96 220 219 152\n";
	printf OUT "\@ 97 222 219 150\n";
	printf OUT "\@ 98 223 219 148\n";
	printf OUT "\@ 99 224 219 147\n";
	printf OUT "\@ 100 226 219 145\n";
	printf OUT "\@ 101 227 219 143\n";
	printf OUT "\@ 102 228 219 141\n";
	printf OUT "\@ 103 229 219 140\n";
	printf OUT "\@ 104 231 218 138\n";
	printf OUT "\@ 105 232 218 136\n";
	printf OUT "\@ 106 233 218 135\n";
	printf OUT "\@ 107 235 218 133\n";
	printf OUT "\@ 108 236 218 131\n";
	printf OUT "\@ 109 237 218 129\n";
	printf OUT "\@ 110 239 218 128\n";
	printf OUT "\@ 111 240 218 126\n";
	printf OUT "\@ 112 240 217 124\n";
	printf OUT "\@ 113 240 216 122\n";
	printf OUT "\@ 114 241 215 119\n";
	printf OUT "\@ 115 241 215 117\n";
	printf OUT "\@ 116 241 214 115\n";
	printf OUT "\@ 117 241 213 113\n";
	printf OUT "\@ 118 242 212 111\n";
	printf OUT "\@ 119 242 211 109\n";
	printf OUT "\@ 120 242 211 107\n";
	printf OUT "\@ 121 242 210 104\n";
	printf OUT "\@ 122 243 209 102\n";
	printf OUT "\@ 123 243 208 100\n";
	printf OUT "\@ 124 243 208 98\n";
	printf OUT "\@ 125 243 207 96\n";
}
