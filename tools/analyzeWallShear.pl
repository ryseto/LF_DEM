#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [output_interval] [xz_shift]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;
use Getopt::Long;
use POSIX;

my $particle_data = $ARGV[0];

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$j = index($name, 'alls_1', 1);
$initconfig = substr($name, 0, $j+6);

printf "$initconfig\n";


open (IN_CONFIG, "< ${initconfig}.dat");
$line = <IN_CONFIG>;
$line = <IN_CONFIG>;
($buf, $np1, $np2, $vf, $lx, $ly, $lz, $np_in, $np_out, $z_bot, $z_top) = split(/\s+/, $line);

$np_mov = $np1+$np2;
printf "$z_bot $z_top\n";

# np1 np2 vf lx ly lz np_in np_out radius_in radius_out
close(IN_CONFIG);
#exit;
# Create output file name

$output = "WRdata_$name.dat";
printf "$output\n";

open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");

&readHeader;

$first = 1;
$output = 1;
$cnt_data = 0;
$shear_strain_steady_state = 5;

$kmax = 15;
$zdiff = ${z_top}-${z_bot};
$v_out = $zdiff*1;
$dz = $zdiff/$kmax;

for ($k = 0; $k < $kmax; $k++) {
	$average[$k] = 0;
	$zposition[$k] = 0;
	$cnt[$k] = 0;
	$particlearea[$k] = 0;
}

while (1) {
	&InParticles;
	last unless defined $line;
	if ($shear_strain > 200) {
		last;
	}
}

for ($k = 0; $k < $kmax; $k++) {
	if ($cnt[$k] != 0) {
		$ave_v_tan[$k] = $average_v[$k]/$cnt[$k];
	}
}


$total_area = 0;
$total_particle_area = 0;
$velo_wall = $zdiff;
for ($k = 0; $k < $kmax; $k++) {
	if ($cnt[$k] != 0) {
		$ave_v_tan = $average_v[$k]/$cnt[$k];
		if ($k < $kmax -1) {
			$gradient_v_tan = ($ave_v_tan[$k+1]-$ave_v_tan[$k])/$dz;
		} else {
			$gradient_v_tan = ($velo_wall-$ave_v_tan[$k])/$dz;
		}

		$z = $z_bot + $dz*$k;
		
		$area = $dz*$lx;
		$density = ($particlearea[$k]/$cnt_data)/$area;
		$total_particle_area += $particlearea[$k]/$cnt_data;
		$zz = ($z - $z_bot)/$zdiff;
		printf OUT "$z  $zz $ave_v_tan[$k] $gradient_v_tan $density $k $particlearea[$k] $cnt_data $area\n";
		$total_area += $area;
	}
}
$phi_check = $total_particle_area/$total_area;
printf "$total_area $phi_check\n";

$calctotal_area=$lx*$zdiff;
printf "$calctotal_area $lx $zdiff \n";



close (OUT);
close (IN_particle);

##################################################################

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

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
	printf "$line";
	if (defined $line) {
		# 1 sys.get_shear_strain()
		# 2 sys.shear_disp
		# 3 getRate()
		# 4 target_stress_input
		# 5 sys.get_time()
		# 6 sys.angle_external_magnetic_field
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
		printf "$shear_strain\n";
		$count_particle = 0;
		for ($i = 0; $i < $np; $i ++){
			$line = <IN_particle>;
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
			$ang[$i] = $angle;
			$radius[$i] = $a;

			if ($shear_strain > $shear_strain_steady_state && $i < $np_mov) {
				$z += $lz/2;
				$v_tan = $vx;
				$f_zpos = ($z - $z_bot)/$dz;
				
				$i_zpos = floor($f_zpos);
				if ($i_zpos >= 0 && $i_zpos < $kmax) {
					$average_v[$i_zpos] += $v_tan;
					$cnt[$i_zpos] ++;
					$particlearea[$i_zpos] += pi*$radius[$i]*$radius[$i];
					$count_particle ++;
				}
			}
			#$posx[$i] = $x;
			#$posy[$i] = $y;
			#$posz[$i] = $z;
			#$velx[$i] = $vx;
			#$vely[$i] = $vy;
			#$velz[$i] = $vz;
			#$omegax[$i] = $ox;
			#$omegay[$i] = $oy;
			#$omegaz[$i] = $oz;
			#$omegay[$i] = $oy;
		}
		printf "particle = $count_particle\n";
	
		if ($shear_strain > $shear_strain_steady_state) {
			$cnt_data ++;
			#		exit;
		}

	}
	
	
}
