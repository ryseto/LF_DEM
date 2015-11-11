#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;

$force_factor = 0.2;

$y_section = 0;
$yap_radius = 1;

$particle_data = $ARGV[0];
$dim = $ARGV[1];


# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

#$interaction_data = "int_${name}.dat";
#printf "$interaction_data\n";
$output = "p_$name.dat";

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

open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");

&readHeader;
$first=1;
$c_traj=0;
$num = 0;

printf OUT "# $np $VF, $Lx $Ly $Lz\n";


while (1){
	&InParticles;
	last unless defined $line;
	&OutData;
	$num ++;
	printf "$shear_rate\n";
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
	for ($i = 0; $i < 23; $i++) {
		$line = <IN_interaction>;
	}
	printf "--- $np, $VF, $Lx, $Ly, $Lz\n";

}

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
    ($buf, $shear_strain, $shear_disp) = split(/\s+/, $line);
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

sub OutData{
#	$npp = $np;
#	for ($i = 0; $i < $np; $i ++){
#		if ($posz[$i] < -$Lz/2+$radius[$i]){
#			$npp += 1;
#		}
#	}
	printf OUT "# $shear_disp\n";
	if ($dim == 2) {
		for ($i = 0; $i < $np; $i ++){
			printf OUT "$i $posx[$i] $posz[$i] $ang[$i] $radius[$i] $ang[$i]\n";
#			if ($posz[$i] < -$Lz/2+$radius[$i]){
#				$x = $posx[$i] + $shear_disp;
#				while ($x > $Lx/2){
#					$x -= $Lx
#				}
#				$z = $posz[$i] + $Lz;
#				printf OUT "$i $x $z $ang[$i] $radius[$i] \n";
#			}
		}
		
	} else {
		for ($i = 0; $i < $np; $i ++){
			printf OUT "$i $posx[$i] $posy[$i] $posz[$i] $radius[$i] $ang[$i]\n";
		}
	}
	

}






