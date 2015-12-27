#!/usr/bin/perl
use Math::Trig;
use IO::Handle;

$particle_data = $ARGV[0];
$output_strain = 9999;
if ($#ARGV >= 1){
	$output_strain = $ARGV[1];
}

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);
$interaction_data = "int_${name}.dat";
$output = "pov_$name.pov";
open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;

while (1){
	&InParticles;
	last unless defined $line;
	&InInteractions;
	if ($shear_strain > $output_strain) {
		last;
	}
}
&OutPovrayFile;

close (OUT);
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
	printf "$np, $VF, $Lx, $Ly, $Lz\n";
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
		
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
		printf "shear_disp\n";
		# shear_rate/shear_rate0
		# shear_rate0 = Fr(0)/(6 pi eta0 a) = 1/6pi
		#$shear_rate = $shear_rate;
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
			($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			$ang[$i] = $angle;
			
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
			
		}
		$c_traj++;
	}
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
	#printf OUTG "\n";
}

sub OutPovrayFile{
	$postext = $Lz/2+2;
	$postext2 = $Lz/2+3;
	$shear_rate_text = int($shear_rate*1e5);
	$shear_rate_text *= 1e-5;
	
	$r = $yap_radius*$radius[0];
	$switch = 0;
	&printHead;
	$xshift = 0;
	$yshift = -0.5;
	for ($i = 0; $i < $np; $i ++) {
		$posy[$i] += $xshift;
		if ($posx[$i] > $Ly/2) {
			$posx[$i] -= $Lx;
		} elsif ($posx[$i] < -$Lx/2) {
			$posx[$i] += $Lx;
		}

		$posy[$i] += $yshift;
		if ($posy[$i] > $Ly/2) {
			$posy[$i] -= $Ly;
		} elsif ($posy[$i] < -$Ly/2) {
			$posy[$i] += $Ly;
		}

		
	}

	
	for ($i = 0; $i < $np; $i ++) {
		$x = $posx[$i];
		$y = $posy[$i];
		$z = $posz[$i];
		$rr = $radius[$i];
		printf OUT "S(<$x, $y, $z>, $rr)\n";
	}
	
	$max_force = 0;
	for ($k = 0; $k < $num_interaction; $k ++){
		if ($max_force < $force[$k]){
			$max_force = $force[$k];
		}
	}
	
	for ($k = 0; $k < $num_interaction; $k ++){
		$r = 0.3*$force[$k]/$max_force;

		if ($contactstate[$k] != 0 && $r > 0) {
			&OutStringPovray($int0[$k], $int1[$k], $r);
		}
	}
}

sub OutBondPov{
	($xs, $ys, $zs, $xe, $ye, $ze, $f, $rs) = @_;
	
	$dist = sqrt(($xs-$xe)**2 + ($ys-$ye)**2 + ($zs-$ze)**2);
	
	if ($dist < 4 && $dist > 0 ){
		$dx = ($xe - $xs)*$rs/$dist;
		$dy = ($ye - $ys)*$rs/$dist;
		$dz = ($ze - $zs)*$rs/$dist;
		
		$xe = $xs + $dx;
		$ye = $ys + $dy;
		$ze = $zs + $dz;
		printf "$dx $dy $dz\n";
		if ($xs != $xe) {
			printf OUT  "B(<$xs, $ys, $zs>, <$xe, $ye, $ze>, $f)\n";
		}
	}
}

sub OutStringPovray{
	($i, $j, $f) = @_;
	$xi = $posx[$i];
	$yi = $posy[$i];
	$zi = $posz[$i];
	$xj = $posx[$j];
	$yj = $posy[$j];
	$zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
		printf OUT  "B(<$xi, $yi, $zi>, <$xj, $yj, $zj>, $f)\n";
	} else {
		if ($xi > $xj+6) {
			
			&OutBondPov($xj, $yj, $zj, $xi-$Lx, $yi, $zi, $f, $radius[$j]);
			&OutBondPov($xi, $yi, $zi, $xj+$Lx, $yj, $zj, $f, $radius[$i]);
		}
		if ($xj > $xi+6) {
			&OutBondPov($xj, $yj, $zj, $xi+$Lx, $yi, $zi, $f, $radius[$j]);
			&OutBondPov($xi, $yi, $zi, $xj-$Lx, $yj, $zj, $f, $radius[$i]);
		}
		
		if ($yi > $yj+6 && abs($xi-$xj) < 6 && abs($zi-$zj) < 6) {
			&OutBondPov($xj, $yj, $zj, $xi, $yi-$Ly, $zi, $f, $radius[$j]);
			&OutBondPov($xi, $yi, $zi, $xj, $yj+$Ly, $zj, $f, $radius[$i]);
		}
		if ($yj > $yi+6 && abs($xi-$xj) < 6 && abs($zi-$zj) < 6) {
			&OutBondPov($xj, $yj, $zj, $xi, $yi+$Ly, $zi, $f, $radius[$j]);
			&OutBondPov($xi, $yi, $zi, $xj, $yj-$Ly, $zj, $f, $radius[$i]);
		}
		if ($zi > $zj+6 && abs($xi-$xj) < 6 && abs($yi-$yj) < 6) {
			&OutBondPov($xj, $yj, $zj, $xi, $yi, $zi-$Lz, $f, $radius[$j]);
			&OutBondPov($xi, $yi, $zi, $xj, $yj, $zj+$Lz, $f, $radius[$i]);
		}
		if ($zj > $zi+6 && abs($xi-$xj) < 6 && abs($yi-$yj) < 6) {
			&OutBondPov($xj, $yj, $zj, $xi, $yi, $zi+$Lz, $f, $radius[$j]);
			&OutBondPov($xi, $yi, $zi, $xj, $yj, $zj-$Lz, $f, $radius[$i]);
		}
		#		if ($yj > $yi+6) {
		#			printf OUT "B(<$xi, $yi+$Ly, $zi>, <$xj, $yj, $zj>, $f)\n";
		#			printf OUT "B(<$xi, $yi, $zi>, <$xj, $yj-$Ly, $zj>, $f)\n";
		#		}
		#
		#		if ($zi > $zj+6) {
		#			printf OUT "B(<$xi, $yi, $zi-$Lz>, <$xj, $yj, $zj>, $f)\n";
		#			printf OUT "B(<$xi, $yi, $zi>, <$xj, $yj, $zj+$Lz>, $f)\n";
		#		}
		#		if ($zj > $zi+6) {
		#			printf OUT "B(<$xi, $yi, $zi+$Lz>, <$xj, $yj, $zj>, $f)\n";
		#			printf OUT "B(<$xi, $yi, $zi>, <$xj, $yj, $zj-$Lz>, $f)\n";
		#		}
	}
}

sub printHead
{
	printf OUT "
	global_settings { max_trace_level 100 }
	#include \"colors.inc\"
	#include \"shapes.inc\"
	#include \"skies.inc\"
	#include \"glass.inc\"
	#include \"woods.inc\"
	#include \"stones.inc\"
	#include \"metals.inc\"\n";
	

#	camera {
#		location <15, 30, 10>
#		right <21/22.5, 0, 0>
#		sky      <0, 0, 1>
#		look_at  <-1, 0, -1.5>
#		angle 45
#	}
	printf OUT "
	camera {
		location <0, 20, 0>
		right <21/22.5, 0, 0>
		sky      <0, 0, 1>
		look_at  <0, 0, 0>
		angle 45
	}
	light_source {
		100
		color White
	}
	background{color rgb <5/255,8/255,42/255>}
	";
	#	background{color rgb <0/255,0/255,0/255>}

	printf OUT "
	#macro S(p, r)
	sphere{ 0, r
		texture {
			pigment{ color rgbft <1, 1, 1, 0.3, 0.5> }
			finish { phong 1 reflection {0.3 metallic 0} }}
		interior {ior 1.4}
		translate p
	}
	#end
	";
	
	printf OUT "
	#macro B(p1, p2, r)
	object{
		cylinder{p1, p2, r}
		pigment{ color rgbf <1,1,0,0.1> }
		finish {
			phong 1
			ambient 0.3
		}}
	#end
	";
	
	#	printf OUT "
	#	#macro B(p1, p2, r)
	#	object{
	#		Round_Cylinder(p1, p2, 0.08, 0.3, 1)
	#
	#	}
	#	#end
	#	";
	
	#	sphere { 0, 0.1
	#		texture{ pigment{ color rgb<1, 1, 0>}}
	#		translate p1
	#	}
	#	sphere { 0, 0.1
	#		texture{ pigment{ color rgb<1, 1, 0>}}
	#		translate p2
	#	}
	#	printf OUT "
	#	plane { <0,0,-1>
	#		, 30
	#  pigment {
	#	  color White
	#  }
	#	}
	#	";
}


