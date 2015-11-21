#!/usr/bin/perl

# Usage:

use Math::Trig;

$particle_data = $ARGV[0];
$initconfig = $ARGV[1];

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '_preprocess', $i-1);
$name = substr($particle_data, $i, $j-$i);

$output_name= "${name}_.dat";
printf "$output_name\n";
open (OUT, "> ${output_name}");
open (IN_particle, "< ${particle_data}");
printf "$initconfig\n";
open (IN, "< ${initconfig}") or die "$!";
$line = <IN>;
printf "aaa $line\n";
printf OUT "$line";
$line = <IN>;
printf "bbb $line\n";
printf OUT "$line";


&readHeader;
$first=1;
$c_traj=0;
$num = 0;

while (1){
	&InParticles;
	
	last unless defined $line;
	$num ++;
}
printf "$np\n";
printf "|$initconfig|\n";

for ($i = 0; $i < $np; $i ++){
	if ($posx[$i] < 0 ){
		$posx[$i] += $Lx;
	}
	if ($posy[$i] < 0 ){
		$posy[$i] += $Ly;
	}
	printf OUT "$posx[$i] $posy[$i] $posz[$i] $radius[$i]\n";
}
close (OUT);

sub readHeader{
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	
	printf "$np, $VF, $Lx, $Ly, $Lz\n";
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
	if ($buf) {
		for ($i = 0; $i < $np; $i ++){
			$line = <IN_particle> ;
			($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
			$h_xzstress, $c_xzstressGU, $b_xzstress, $angle) = split(/\s+/, $line);
			
			if ($num==$num_mathm){
				printf OUTMP "$line";
			}
			
			$radius[$i] = $a;
			$posx[$i] = $x;
			$posy[$i] = $y;
			$posz[$i] = $z;
			
			$ang[$i] = $angle;
			if ($radius_max < $a){
				$radius_max = $a;
			}
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
