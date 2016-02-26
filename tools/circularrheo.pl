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
my $np_movable = 9000;
#my $rout = 81.3781; # 0.80
#my $rout = 82.4148; # 0.78
my $rout = 142.741;
#my $np_movable = 6000;
#my $rout = 118.562;

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
	&yaplotColor;
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
	&InInteractions;
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
		printf "$shear_strain_i $shear_strain \n";
		$viscosity_in =  3*$force_in/(2*$rin);
		$viscosity = 3*$force_out/(2*$rout);
		printf OUT "$shear_strain $shear_rate $viscosity $viscosity_in\n";
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
	for ($i = 0; $i<16; $i++) {
		$line = <IN_particle>;
	}
	for ($i = 0; $i<24; $i++) {
		$line = <IN_interaction>;
	}
}

sub yaplotColor {
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
	#printf OUT "\@8 224 143 0 \n";
	#printf OUT "\@9 67 163 230 \n";
	#printf OUT "\@8 253 105 6 \n";
	#printf OUT "\@9 109 109 109 \n";
	printf OUT "\@10 250 250 250 \n";
	printf OUT "\@11 240 240 240 \n";
	printf OUT "\@12 230 230 230 \n";
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
	$line = <IN_particle>;
	if (defined $line) {
		
		# 1 sys.get_shear_strain()
		# 2 sys.shear_disp
		# 3 getRate()
		# 4 target_stress_input
		# 5 sys.get_time()
		# 6 sys.angle_external_magnetic_field
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
		
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
				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
				$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
				
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
				if ($radius_max < $a) {
					$radius_max = $a;
				}
			}
		}
		if ($rotatingobserver) {
			if ($first) {
				$rout = sqrt($posxtmp[$np-1]*$posxtmp[$np-1] + $posztmp[$np-1]*$posztmp[$np-1]);
				#printf "$rout\n";
			}
			$theta = -0.25*$shear_strain;
			#printf "$theta, $shear_strain, $rout\n";
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
			printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
		}
	} else {
		printf OUT "@ 8\n";
		for ($i = 0; $i < $np; $i++) {
			printf OUT "r $radius[$i]\n";
			printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
		}
		printf OUT "r 1\n";
		printf OUT "@ 0\n";
		for ($i = $np_movable; $i < $np; $i+=10) {
			printf OUT "c $posx[$i] -0.01 $posz[$i] \n";
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
				printf OUT "c $posx_init[$i] 0.01 $posz_init[$i] \n";
			} else {
				printf OUT "c $posx_init[$i] $posy_init[$i] $posz_init[$i] \n";
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
