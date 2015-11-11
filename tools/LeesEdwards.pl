#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;

$particle_data = $ARGV[0];

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

#printf "$interaction_data\n";
$output = "le_$name.yap";

my $pos;
$pos = index($name, "_m");
if ($pos != -1) {
	$mag = 1;
} else {
	$mag = 0;
}
if ($mag) {
	$outputmp = "magprofile_$name.dat";
	open (OUTMAGPROF, "> $outputmp");
}
open (IN_particle, "< ${particle_data}");
open (OUT, "> ${output}");

&readHeader;
$first=1;
$c_traj=0;
$num = 0;
printf OUT "\@0 0 0 0 \n";
printf OUT "\@1 50 100 205 \n";
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

$shear_strain_next = 1;
$shear_strain_previous = 0;
while (1){
	&InParticles;
	last unless defined $line;
	&OutYaplotData;
	$num ++;
}

close (OUT);
close (IN_particle);

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
		
		if ($mag) {
			($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			$h_xzstress, $c_xzstressGU, $b_xzstress, $mx, $my, $mz, $ms) = split(/\s+/, $line);
			$magmom_x[$i] = $mx;
			$magmom_y[$i] = $my;
			$magmom_z[$i] = $mz;
			$magsusceptibility[$i] = $ms;
			$mm[$i] = sqrt($mx*$mx+$my*$my+$mz*$mz);
		} else {
			($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			$ang[$i] = $angle;
		}
		$radius[$i] = $a;
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
		$ang[$i] = $angle;
		if ($radius_max < $a){
			$radius_max = $a;
		}
	}
	$c_traj++;
}

sub calcsqdist {
	($x1, $y1, $z1, $x2, $y2, $z2) = @_;
	$dist = ($x1-$x2)*($x1-$x2);
	$dist += ($y1-$y2)*($y1-$y2);
	$dist += ($z1-$z2)*($z1-$z2);
	return $dist;
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
	
	printf OUT "y 1\n";
	printf OUT "@ 8\n";
	$r = $radius[0];
	printf OUT "r $r\n";
	$switch = 0;
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
	}
	
	


		
	
	$switch = 0;
	printf OUT "@ 8\n";
	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}

		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] + $Lx;
		printf OUT "c $x_image $posy[$i] $posz[$i] \n";
	}
	$switch = 0;
	printf OUT "@ 8\n";

	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}

		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] - $Lx;
		printf OUT "c $x_image $posy[$i] $posz[$i] \n";
	}
	$switch = 0;
	printf OUT "@ 8\n";
	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}

		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] + $shear_disp;
		$z_image = $posz[$i] + $Lz;
		printf OUT "c $x_image $posy[$i] $z_image \n";
	}
	$switch = 0;
	printf OUT "@ 8\n";

	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] + $shear_disp -$Lx;
		$z_image = $posz[$i] + $Lz;
		printf OUT "c $x_image $posy[$i] $z_image \n";
	}
	$switch = 0;
	printf OUT "@ 8\n";
	
	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] + $shear_disp +$Lx;
		$z_image = $posz[$i] + $Lz;
		printf OUT "c $x_image $posy[$i] $z_image \n";
	}
	
	$switch = 0;
	printf OUT "@ 8\n";
	
	
	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] - $shear_disp;
		$z_image = $posz[$i] - $Lz;
		printf OUT "c $x_image $posy[$i] $z_image \n";
	}
	
	$switch = 0;
	printf OUT "@ 8\n";

	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] - $shear_disp + $Lx;
		$z_image = $posz[$i] - $Lz;
		printf OUT "c $x_image $posy[$i] $z_image \n";
	}
	
	$switch = 0;
	printf OUT "@ 8\n";

	printf OUT "r $r\n";
	for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $yap_radius*$radius[$i];
			printf OUT "r $r\n";
		}
		if ($mag) {
			if ($switch == 0 &&
				$magsusceptibility[$i] <= 0){
					printf OUT "@ 9\n";
					$switch = 1;
				}
		}
		#		if ($i % 100 == 0){
		#			$col = $i/100 + 2;
		#			printf OUT "@ $col\n";
		#		}
		$x_image = $posx[$i] - $shear_disp - $Lx;
		$z_image = $posz[$i] - $Lz;
		printf OUT "c $x_image $posy[$i] $z_image \n";
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
	
	

	
#	if ($Ly == 0){
#		printf OUT "y 6\n";
#		printf OUT "@ 1\n";
#		for ($i = 0; $i < $np; $i ++){
#			&OutCross($i);
#		}
#	}
#	
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
	printf OUT "y 8\n";
	printf OUT "@ 4\n";
	printf OUT "r 0.5\n";
	if ($mag) {
		$fieldx = 20*cos($fieldangle);
		$fieldz = 20*sin($fieldangle);
		
		$xs = 0;
		$xe = $xs + $fieldx;
		$zs = 0;
		$ze = $zs + $fieldz;
		$xx = $xs + 3.2;
		$zz = $zs + 3.2;
		
		
		printf OUT "s  $xs -0.1 $zs  $xe 0 $ze \n";
		
	}
	
	
	printf OUT "y 7\n";
	printf OUT "@ 0\n";
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
	$ybox = -0.1;
	printf OUT "r 0.2\n";
	if($Ly == 0){
		$lx2 = $Lx/2+1;
		$ly2 = $Ly/2+1;
		$lz2 = $Lz/2+1;
		#printf OUT "p 4 -$lx2 $yb $lz2 $lx2 $yb $lz2 $lx2 $yb -$lz2 -$lx2 $yb -$lz2\n";
		printf OUT "l -$lx2 $ybox $lz2    $lx2 $ybox $lz2\n";
		printf OUT "l -$lx2 $ybox -$lz2   $lx2 $ybox -$lz2\n";
		printf OUT "l -$lx2 $ybox -$lz2  -$lx2 $ybox $lz2\n";
		printf OUT "l $lx2 $ybox $lz2   $lx2 $ybox -$lz2\n";
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
