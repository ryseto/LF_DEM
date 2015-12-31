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

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

# Create output file name

$output = "CWGdata_$name.dat";
printf "$output\n";

open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");

&readHeader;

$first = 1;
$output = 1;
$cnt_data = 0;
$shear_strain_steady_state = 5;

$kmax = 10;
$r_in = 22.8447 + 1;
$r_out = 45.6894 - 1;
$rdiff = ${r_out}-${r_in};
$dr = $rdiff/$kmax;

for ($k = 0; $k <= $kmax; $k++) {
	$average[$k] = 0;
	$radialposition[$k] = 0;
	$cnt[$k] = 0;
}

while (1) {
	&InParticles;
	last unless defined $line;
}


for ($k = 0; $k <= $kmax; $k++) {
	if ($cnt[$k] != 0) {
		$ave_v_tan = $average_v[$k]/$cnt[$k];
		$r = $r_in + $dr*$k;
		$rn = $r + $dr;
		$area = pi*($rn*$rn - $r*$r);
		$density = (pi*$cnt[$k]/$cnt_data)/$area;
		
		printf OUT "$r $ave_v_tan $density\n";
	}
}


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
	if (defined $line) {
		# 1 sys.get_shear_strain()
		# 2 sys.shear_disp
		# 3 getRate()
		# 4 target_stress_input
		# 5 sys.get_time()
		# 6 sys.angle_external_magnetic_field
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress) = split(/\s+/, $line);
		printf "$shear_strain\n";
		
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
				if ($shear_strain > $shear_strain_steady_state) {
					$pos_r[$i] = sqrt($x*$x + $z*$z);
					$v_dot_rvec = $vx*$x + $vz*$z;
					$v_tan_x = $vx - ${v_dot_rvec} * $x /$pos_r[$i];
					$v_tan_z = $vz - ${v_dot_rvec} * $z /$pos_r[$i];
					$v_tan[$i] = sqrt($v_tan_x*$v_tan_x+ $v_tan_z*$v_tan_z);

					$i_rpos = int(($pos_r[$i] - $r_in)/$dr);
					if ($i_rpos > 0 && $i_rpos <= $kmax) {
						$average_v[$i_rpos] += $v_tan[$i];
						$cnt[$i_rpos] ++;
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
		}
		if ($shear_strain > $shear_strain_steady_state) {
			$cnt_data ++;
		}
		
	}
}
