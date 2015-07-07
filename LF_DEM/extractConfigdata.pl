#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;

$particle_data = $ARGV[0];
if ($#ARGV == 1){
	$xz_shift = $ARGV[1];
}


# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);
#printf "$interaction_data\n";
$outputfilename = "config_$name.dat";

open (IN_data, "< data_${name}.dat");

$angle_magneticfield_old = -1;
$cnt = 0;



while (1) {
	$line = <IN_data>;
	($d1, $d2, $d3, $d4, $d5, $d6, $d7, $d8, $d9, $d10,
	$d11, $d12, $d13, $d14, $d15, $d16, $d17, $d18, $d19, $d20,
	$d21, $d22, $d23, $d24, $d25, $d26, $d27, $d28, $d29, $d30,
	$d31, $d32, $d33, $d34, $d35, $d36, $d37, $d38, $d39, $d40) = split(/\s+/, $line);
	last unless defined $line;
	if ($d1 != "#") {
		$angle_magneticfield = $d38;
		#printf "ang = $angle_magneticfield_old\n";
		$time = $d1;
		if ($angle_magneticfield_old != -1
			&& $angle_magneticfield != $angle_magneticfield_old) {
				$time[$cnt] = $time_old;
				$angle[$cnt] = $angle_magneticfield_old;
				$cnt++;
			}
		$angle_magneticfield_old =	$angle_magneticfield;
		#printf "$angle_magneticfield_old\n" ;
		$time_old = $time;
	}
}

open (OUT, "> ${outputfilename}");
open (IN_particle, "< ${particle_data}");
&readHeader;
$first=1;
$c_traj=0;
$num = 0;

$cnt = 0;
while (1){
	&InParticles;
	last unless defined $line;
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
	if (defined $line){
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress, $time) = split(/\s+/, $line);
		if ($time == $time[$cnt]){
			$output = 1;
			$ang =$angle[$cnt];
			$cnt++;
			printf OUT "# $ang\n";
		} else {
			$output = 0;
		}
		if ($output == 1) {
			printf "$ang\n";
		}
		for ($i = 0; $i < $np; $i ++){
			$line = <IN_particle> ;
			if ($output == 1) {
				($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
				$h_xzstress, $c_xzstressGU, $b_xzstress, $mx, $my, $mz, $ms) = split(/\s+/, $line);
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
				#printf OUT "$line";
				printf OUT "$ip $a $x $y $z $vx $vy $vz $ox $oy $oz $h_xzstress $c_xzstressGU $b_xzstress $mx $my $mz $ms\n";
			}
		}
	}
}




