#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;
$y_section = 0;
$yap_radius = 1;
$avenumber = 1;
$mag = 1; ### magnetic simulation

$magnetoffset_y = -0.01;
$particle_data = $ARGV[0];
$output_interval = 1;
if ($#ARGV >= 1){
	$output_interval = $ARGV[1];
}
if ($#ARGV == 2){
	$xz_shift = $ARGV[2];
}

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);
$interaction_data = "int_${name}.dat";

printf "${name}\n";
if ($avenumber == 1){
	$yaplotfilename = "y_${name}.yap";
} else {
	$yaplotfilename = "y_${name}_av.yap";
}
printf "$yaplotfilename\n";

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
	if ($d1 > 2) {
		${sum_fmax} += $d29;
		${cnt} ++;
	}
}

open (OUT, "> ${yaplotfilename}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
$first=1;
$c_traj=0;
$num = 0;

&outputColorTable;

$shear_strain_next = 1;
$shear_strain_previous = 0;
$cnt_interval = 0;

$avecount = 0;
$output = 0;
while (1){
	&InParticles;
	last unless defined $line;
	if ($mag == 0) {
		&InInteractions;
	}
	if ($avenumber >= 2) {
		for ($i=0; $i < $np; $i++) {
			$aposx[$avecount][$i] = $posx[$i];
			$aposz[$avecount][$i] = $posz[$i];
		}
		for ($i=0; $i < $np; $i++) {
			$averagex[$i] = $posx[$i];
			$averagez[$i] = $posz[$i];
		}
		
		for ($j=1; $j<$avenumber; $j++){
			$jj = $avecount - $j;
			if ($jj < 0) {
				$jj += $avenumber;
			}
			
			for ($i=0; $i < $np; $i++) {
				$xpd[$i] = $aposx[$jj][$i];
				$zpd[$i] = $aposz[$jj][$i];
				if (abs($xpd[$i] - $posx[$i]) > 0.5*$Lx) {
					if ($xpd[$i] > $posx[$i]){
						$xpd[$i] -= $Lx;
					} else {
						$xpd[$i] += $Lx;
					}
				}
				if (abs($zpd[$i] - $posz[$i]) > 0.5*$Lz) {
					if ($zpd[$i] > $posz[$i]){
						$zpd[$i] -= $Lz;
					} else {
						$zpd[$i] += $Lz;
					}
				}
				$averagex[$i] += $xpd[$i];
				$averagez[$i] += $zpd[$i];
			}
		}
		for ($i=0; $i < $np; $i++) {
			$posx[$i] = $averagex[$i]/$avenumber;
			$posy[$i] = 0;
			$posz[$i] = $averagez[$i]/$avenumber;
		}
		#		for ($i=0; $i < $np; $i++) {
		#			$stdev[$i] = 0;
		#			for ($j=0; $j<$avenumber; $j++){
		#				$jj = $avecount - $j;
		#				if ($jj < 0) {
		#					$jj += $avenumber;
		#				}
		#				$xpd[$i] = $aposx[$jj][$i];
		#				$zpd[$i] = $aposz[$jj][$i];
		#				if (abs($xpd[$i] - $posx[$i]) > 0.5*$Lx) {
		#					if ($xpd[$i] > $posx[$i]){
		#						$xpd[$i] -= $Lx;
		#					} else {
		#						$xpd[$i] += $Lx;
		#					}
		#				}
		#				if (abs($zpd[$i] - $posz[$i]) > 0.5*$Lz) {
		#					if ($zpd[$i] > $posz[$i]){
		#						$zpd[$i] -= $Lz;
		#					} else {
		#						$zpd[$i] += $Lz;
		#					}
		#				}
		#
		#				$stdev[$i] += sqrt(($xpd[$i]-$posx[$i])**2+($zpd[$i]-$posz[$i])**2);
		#			}
		#			$stdev[$i] *= 1/$avenumber;
		#		}
		if ($avecount == $avenumber-1){
			$output = 1;
		}
	} else {
		$output = 1;
	}
	if ($output == 1) {
		if ($time > 0) {
			&OutYaplotData;
		}
		#		&OutputTxtFile;
	}
	$cnt_interval ++;
	$avecount ++;
	
	if ($avecount == $avenumber){
		$avecount = 0;
	}
}

open (OUTLAST, "> LastConfig.dat");
printf OUTLAST "$shear_disp\n";

for ($i = 0; $i < $np; $i++){
	$xx = 0.5*$Lx + $posx[$i];
	$yy = 0.5*$Ly + $posy[$i];
	$zz = 0.5*$Lz + $posz[$i];
	#	$mx = $magmom_x[$i];
	#	$my = $magmom_y[$i];
	#	$mz = $magmom_z[$i];
	$ms = $magsusceptibility[$i];
	printf OUTLAST "$xx $yy $zz $radius[$i] 0 0 0 $ms \n";
}
close (OUTLAST);

if ($mag) {
	close (OUTMAGPROF);
}

close (OUT);
close (OUTE);
close (IN_particle);
close (IN_interaction);

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
	#printf "$np, $VF, $Lx, $Ly, $Lz\n";
}

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
	if (defined $line){
		# 1 sys.get_shear_strain()
		# 2 sys.shear_disp
		# 3 getRate()
		# 4 target_stress_input
		# 5 sys.get_time()
		# 6 sys.angle_external_magnetic_field
		if ($mag) {
			($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress, $time, $fieldangle) = split(/\s+/, $line);
		} else {
			($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
		}
		# shear_rate/shear_rate0
		# shear_rate0 = Fr(0)/(6 pi eta0 a) = 1/6pi
		printf "$time\n";
		$shear_rate = $shear_rate;
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
			#			($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			#			$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			#
			if ($mag) {
				# ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz, $ms, $brownian_pressure, $contact_pressure) = split(/\s+/, $line);
				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ms) = split(/\s+/, $line);
				#					$magmom_x[$i] = $mx;
				#$magmom_y[$i] = $my;
				#$magmom_z[$i] = $mz;
				$magsusceptibility[$i] = $ms;
				#					$bpressure[$i] = $brownian_pressure;
				#				$con_pressure[$i] = $contact_pressure;
				#					$mm[$i] = sqrt($mx*$mx+$my*$my+$mz*$mz);
			} else {
				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
				$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
				$ang[$i] = $angle;
			}
			#		if (true){
			#			#printf OUTMP "$line";
			#			printf OUTMP "$i $x $y $z $a\n";
			#		}
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
			$posx[$i] = $x;
			$posy[$i] = $y;
			$posz[$i] = $z;
			$velx[$i] = $vx;
			$vely[$i] = $vy;
			$velz[$i] = $vz;
			$omegax[$i] = $ox;
			$omegay[$i] = $oy;
			$omegaz[$i] = $oz;
			$omegay[$i] = $oy;
			if ($radius_max < $a){
				$radius_max = $a;
			}
			
		}
		$c_traj++;
	}
}

sub outputColorTable {
	printf OUT "\@0 0 0 0 \n";
	printf OUT "\@1 50 100 205 \n";
	#printf OUT "\@1 255 255 255 \n";
	printf OUT "\@2 200 200 200 \n";
	printf OUT "\@3 255 127 0\n";
	printf OUT "\@4 50 200 50\n";
	printf OUT "\@5 255 100 100\n";
	printf OUT "\@6 50 150 255\n";
	printf OUT "\@7 255 255 0\n";
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
	
	$velocity0 = $shear_rate;
	$force0 = 6*pi*($shear_rate);
	
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
		if ($output == 1) {
			($i, $j, $contact, $nx, $ny, $nz, #1---6
			$gap, $f_lub_norm, # 7, 8
			$f_lub_tan_x, $f_lub_tan_y, $f_lub_tan_z, # 9, 10, 11
			$fc_norm, # 12
			$fc_tan_x, $fc_tan_y, $fc_tan_z, # 13, 14, 15
			$fr_norm, $s_xF
			) = split(/\s+/, $line);
			# $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
			$int0[$k] = $i;
			$int1[$k] = $j;
			$contactstate[$k] = $contact;
			#			$vx1 = $velocity0*($velx[$i] + $radius[$i]*($ny*$omegaz[$i]-$nz*$omegay[$i]));
			#			$vy1 = $velocity0*($vely[$i] + $radius[$i]*($nz*$omegax[$i]-$nx*$omegaz[$i]));
			#			$vz1 = $velocity0*($velz[$i] + $radius[$i]*($nx*$omegay[$i]-$ny*$omegax[$i]));
			#			$vx2 = $velocity0*($velx[$j] - $radius[$j]*($ny*$omegaz[$j]-$nz*$omegay[$j]));
			#			$vy2 = $velocity0*($vely[$j] - $radius[$j]*($nz*$omegax[$j]-$nx*$omegaz[$j]));
			#			$vz2 = $velocity0*($velz[$j] - $radius[$j]*($nx*$omegay[$j]-$ny*$omegax[$j]));
			#			$rvelx[$k] = ($vx2 - $vx1);
			#			$rvely[$k] = ($vy2 - $vy1);
			#			$rvelz[$k] = ($vz2 - $vz1);
			#			$rvel_dot_norm = $rvelx[$k]*$nx + $rvely[$k]*$ny + $rvelz[$k]*$nz;
			#			$rvelx[$k] -= $rvel_dot_norm*$nx;
			#			$rvely[$k] -= $rvel_dot_norm*$ny;
			#			$rvelz[$k] -= $rvel_dot_norm*$nz;
			#			$domega[$k] = $omegay[$i] - $omegay[$j];
			$f_normal = $fc_norm +$f_lub_norm + $fr_norm;
			#		$fx[$k] = $force0*($fc_tan_x + $f_lub_tan_x + $f_normal*$nx);
			#		$fy[$k] = $force0*($fc_tan_y + $f_lub_tan_y + $f_normal*$ny);
			#		$fz[$k] = $force0*($fc_tan_z + $f_lub_tan_z + $f_normal*$nz);
			#		$fx[$k] = $force0*($fc_tan_x + $f_lub_tan_x);
			#$fy[$k] = $force0*($fc_tan_y + $f_lub_tan_y);
			#$fz[$k] = $force0*($fc_tan_z + $f_lub_tan_z);
			#			$fx[$k] = $force0*($f_lub_tan_x);
			#			$fy[$k] = $force0*($f_lub_tan_y);
			#			$fz[$k] = $force0*($f_lub_tan_z);
			#$fx[$k] = $f_lub_tan_x; #+ $f_normal*$nx;
			#$fy[$k] = $f_lub_tan_y;# + $f_normal*$ny;
			#$fz[$k] = $f_lub_tan_z;#+ $f_normal*$nz;
			$F_lub[$k] = $f_lub_norm;
			$Fc_n[$k] = $fc_norm;
			$Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
			#			$f_normal = $fc_norm + $fr_norm + $f_lub_norm;
			#$f_normal = $f_lub_norm;
			$force[$k] = $f_normal;
			
			$S_bf[$k] =  $s_xF;
			#			$fricstate[$k] = $friction;
			$nrvec_x[$k] = $nx;
			$nrvec_y[$k] = $ny;
			$nrvec_z[$k] = $nz;
			#			$ft_x[$k] = $f_lub_tan_x + $fc_tan_x;
			#			$ft_y[$k] = $f_lub_tan_y + $fc_tan_y;
			#			$ft_z[$k] = $f_lub_tan_z + $fc_tan_z;
			
			$Gap[$k] = $gap;
			#	printf OUTG "$gap ";
		}
	}
	
	#printf OUTG "\n";
}

sub OutputTxtFile{
	if ($firstMag == 0){
		printf OUTMAG "\n";
	} else {
		$firstMag = 0;
	}
	
	for ($i = 0; $i < $np; $i ++){
		#		printf OUTMAG "$i $posx[$i] $posy[$i] $posz[$i] $magmom_x[$i] $magmom_y[$i] $magmom_z[$i]\n";
		printf OUTMAG ("%3d %.5f %.5f %.5f %.5f %.5f %.5f \n", $i, $posx[$i], $posy[$i], $posz[$i], $magmom_x[$i], $magmom_y[$i], $magmom_z[$i]);
		
	}
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
	#	printf OUT "@ 3\n";
	#&OutMeter($shear_rate, 0.3, $Lx/2+2, -$Lz/2+2);
	#printf OUT "@ 4\n";
	#&OutMeter($shear_stress, 2, -$Lx/2-3.5, -$Lz/2+2);
	#&OutStress($shear_stress, 25);
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
	printf OUT "@ 3\n";
	$r = $yap_radius*$radius[0];
	printf OUT "r $r\n";
	$switch = 0;
	for ($i = 0; $i < $np; $i ++){
		#		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
		#			$r = $yap_radius*$radius[$i];
		#			printf OUT "r $r\n";
		#		}
		if ($mag) {
			if ($switch == 0 &&
				($magsusceptibility[$i] == -1
				|| $magsusceptibility[$i] == 99)) {
					printf OUT "@ 8\n";
					$switch = 1;
				}
		}
		
		#		if ($stdev[$i] > 0.2) {
		#			printf OUT "@ 5\n";
		#		} else {
		#			printf OUT "@ 8\n";
		#		}
		
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		
		printf OUT "r $radius[$i]\n";
 		printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
	}
	
	if ($mag) {
		printf OUT "y 3\n";
		printf OUT "r 0.5\n";
		printf OUT "@ 0\n";
		for ($i = 0; $i < $np; $i ++){
			if ($magsusceptibility[$i] == 99
				|| $magsusceptibility[$i] == 101) {
					printf OUT "c $posx[$i] -0.02 $posz[$i] \n";
				}
		}
	}
	
	if (0) {
		printf OUT "y 2\n";
		printf OUT "r 0.3\n";
		printf OUT "@ 5\n"; # static
		for ($k = 0; $k < $num_interaction; $k ++){
			if ($contactstate[$k] == 2 && 0.035* $force[$k] > 0.05) {
				&OutString2($int0[$k],  $int1[$k]);
				#&OutContact($int0[$k], $int1[$k], $contactstate[$k]);
			}
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
	#						}
	#
	#				}
	#			}
	#	    }
	#	printf OUT "y 3\n";
	#	printf OUT "r 0.4\n";
	#	printf OUT "@ 7\n"; # static
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		$force = $Fc_n[$k];
	#		if ($F_lub[$k] < 0) {
	#			$force += - $F_lub[$k];
	#		}
	#		if ($fricstate[$k] == 0) {
	#			&OutString2($int0[$k],  $int1[$k]);
	#		}
	#	}
	
	
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
	#	$force_factor = 0.035;
	#	$force_factor = 0.02;
	#		printf OUT "y 3\n";
	#		printf OUT "@ 6\n";
	#		for ($k=0; $k<$num_interaction; $k++){
	#			#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
	#			#$force = $Fcol[$k];
	#			#$force = $Fc_n[$k];
	#			if ($force[$k] < 0){
	#
	#				$forceA = -$force[$k];
	#			$string_width = ${force_factor}*$forceA;
	#				#&OutString2($int0[$k], $int1[$k]);
	#				&OutString_width($int0[$k], $int1[$k]);
	#			}
	#		}
	if (0) {
		
		printf OUT "y 4\n";
		printf OUT "@ 7\n";
		for ($k = 0; $k < $num_interaction; $k ++){
			#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
			if ($force[$k] > 0){
				if ( $contactstate[$k] == 2 ){
					printf OUT "@ 7\n"; # static
				}
				if ( $contactstate[$k] == 3 ){
					printf OUT "@ 7\n"; # dynamic
				}
				$forceA = $force[$k];
				$string_width = ${force_factor}*$forceA;
				#&OutString2($int0[$k], $int1[$k]);
				&OutString_width($int0[$k], $int1[$k]);
			}
		}
	}
	#	printf OUT "y 4\n";
	#	printf OUT "@ 4\n";
	#	printf OUT "r 0.5\n";
	#	for ($k = 0; $k < $num_interaction; $k ++){
	#		&OutEnergyDissipation($int0[$k], $int1[$k], $k);
	#	}
	
	
	#		printf OUT "y 2\n";
	#		printf OUT "@ 6\n";
	#		printf OUT "r 0.2\n";
	#		for ($k = 0; $k < $num_interaction; $k ++){
	#	#		if ($Gap[$k] < 0) {
	#	#			&OutString2($int0[$k], $int1[$k]);
	#	#		}
	#			if ($is_contact[$k] == 1) {
	#				&OutContact($int0[$k], $int1[$k]);
	#				#&OutString2($int0[$k], $int1[$k]);
	#			}
	#		}
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
	
	if (0) {
		if ($mag) {
			printf OUT "y 6\n";
			printf OUT "@ 0\n";
			printf OUT "r 0.5\n";
			printf OUT "a 1\n";
			for ($i = 0; $i < $np; $i ++){
				&OutMagMoment($i);
			}
		} else {
			if ($Ly == 0){
				printf OUT "y 6\n";
				printf OUT "@ 1\n";
				for ($i = 0; $i < $np; $i ++){
					OutCross($i);
				}
			}
		}
	}
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
	
	if (0) {
		$fieldx = 3*sin($fieldangle);
		$fieldz = 3*cos($fieldangle);
		
		$xs = 23;
		$xe = $xs + $fieldx;
		$zs = 0;
		$ze = $zs + $fieldz;
		$xx = $xs + 3.2;
		$zz = $zs + 3.2;
		
		printf OUT "r 0.1 \n";
		printf OUT "s  $xs 0 $zs  $xe 0 $ze \n";
		printf OUT "t  $xx 0 $zs x \n";
		printf OUT "t  $xs 0 $zz z \n";
		
	}
	
	
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
	
}

sub OutNvec {
	($k) = @_;
	$nx = $nrvec_x[$k];
	$ny = $nrvec_y[$k];
	$nz = $nrvec_z[$k];
	printf OUT2 "$nx $ny $nz\n";
}

sub OutEnergyDissipation {
	($i, $j, $k) = @_;
	$xi = $posx[$i];
	$yi = $posy[$i] - 0.01;
	$zi = $posz[$i];
	$xj = $posx[$j];
	$yj = $posy[$j] - 0.01;
	$zj = $posz[$j];
	$velfactor = 0.05;
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2 ;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
		$xs = ($radius[$j]*$xi + $radius[$i]*$xj)/($radius[$j]+$radius[$i]);
		$ys = ($radius[$j]*$yi + $radius[$i]*$yj)/($radius[$j]+$radius[$i]);
		$zs = ($radius[$j]*$zi + $radius[$i]*$zj)/($radius[$j]+$radius[$i]);
		$xe = $xs + $velfactor*$rvelx[$k];
		$ye = $ys + $velfactor*$rvely[$k];
		$ze = $zs + $velfactor*$rvelz[$k];
		$xe2 = $xi + $radius[$i]*$nrvec_x[$k];
		$ye2 = $yi + $radius[$i]*$nrvec_y[$k];
		$ze2 = $zi + $radius[$i]*$nrvec_z[$k];
		#$energydis = abs($rvelx[$k]*$fx[$k] + $rvely[$k]*$fy[$k] + $rvelz[$k]*$fz[$k]);
		#		if ($shear_strain > 0.5) {
		#			printf "$energydis\n";
		#		}
		$energydis = 5000*$energy_dis[$k];
		if ($energydis > 1){
			$energydisINT = int ($energydis);
			#printf "$energydisINT\n";
			$energyColor = 35 - $energydisINT ;
			if ($energyColor < 10){
				#printf "$energydisINT\n";
				$energyColor = 10;
			}
		} else {
			$energyColor = 35;
		}
		#		printf OUT "r 0.1\n";
		#		printf OUT "s $xs $ys $zs $xe $ye $ze\n";
		#		printf OUT "s $xi $yi $zi $xe2 $ye2 $ze2\n";
		if ($contactstate[$k] == 2){
			#printf OUT "@ 35\n";
			#printf OUT "c $xs $ys $zs \n";
		} else {
			if ($energydisINT != 0) {
				#printf OUT "r $energydis \n";
				printf OUT "@ $energyColor\n";
				printf OUT "c $xs $ys $zs \n";
			}
		}
	}
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
	#//	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
	#if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
	printf OUT "s $xi $yi $zi $xj $yj $zj\n";
	#}
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
	$arrowhead = 2;
	$arrowwidth = 1;
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

sub OutMagMoment {
	($i) = @_;
	$a = $radius[$i];
	$xi = $posx[$i];
	$yi = $posy[$i]+$magnetoffset_y;
	$zi = $posz[$i];
	#	printf "$mm\n";
	$m0 = 1;
	if (abs($mm[$i]) > 0) {
		#		$xa = $xi - $magmom_x[$i]/$mm[$i] ;
		#		$ya = $yi - $magmom_y[$i]/$mm[$i] ;
		#		$za = $zi - $magmom_z[$i]/$mm[$i] ;
		#		$xb = $xi + $magmom_x[$i]/$mm[$i];
		#		$yb = $yi + $magmom_y[$i]/$mm[$i];
		#		$zb = $zi + $magmom_z[$i]/$mm[$i];
		$xa = $xi - $magmom_x[$i]/$m0;
		$ya = $yi - $magmom_y[$i]/$m0;
		$za = $zi - $magmom_z[$i]/$m0;
		$xb = $xi + $magmom_x[$i]/$m0;
		$yb = $yi + $magmom_y[$i]/$m0;
		$zb = $zi + $magmom_z[$i]/$m0;
		printf OUT "s $xa $ya $za $xb $yb $zb\n";
	}
}



