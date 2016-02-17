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

$j = index($name, 'cylinders0.5_', 1);
$initconfig = substr($name, 0, $j+14);

printf "$initconfig\n";


open (IN_CONFIG, "< ${initconfig}.dat");
$line = <IN_CONFIG>;
$line = <IN_CONFIG>;
($buf, $np1, $np2, $vf, $lx, $ly, $lz, $np_in, $np_out, $radius_in, $radius_out) = split(/\s+/, $line);

$np_mov = $np1+$np2;

printf "$radius_in $radius_out\n";
# np1 np2 vf lx ly lz np_in np_out radius_in radius_out
close(IN_CONFIG);
#exit;
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

$kmax = 8;
$r_in = $radius_in;
$r_out = $radius_out;
$rdiff = ${r_out}-${r_in};

$v_out = ($radius_out-$radius_in)*1;


$dr = $rdiff/$kmax;

printf "$radius_in $radius_out $rdiff $dr \n";
#exit;
for ($k = 0; $k < $kmax; $k++) {
	$average[$k] = 0;
	$radialposition[$k] = 0;
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


for ($k = 0; $k < $kmax; $k++) {
	if ($cnt[$k] != 0) {
		$ave_v_tan = $average_v[$k]/$cnt[$k];
		if ($k < $kmax -1) {
			$gradient_v_tan = ($ave_v_tan[$k+1]-$ave_v_tan[$k])/$dr;
		} else {
			$gradient_v_tan = 0;
		}
		$r = $r_in + $dr*$k + 0.5*$dr;
		$rnorm = ($r - $r_in)/($r_out - $r_in);
		$rn = $r + $dr;
		$area = pi*($rn*$rn - $r*$r);
		$density = ($particlearea[$k]/$cnt_data)/$area;
		
		printf OUT "$r $ave_v_tan[$k] $gradient_v_tan $density $rnorm\n";
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
				$pos_r2 = $x*$x + $z*$z;
				$pos_r = sqrt($pos_r2);
				$v_tan = ((-$vx*$z + $vz*$x)/$pos_r)/$v_out;
				$f_rpos = ($pos_r - $r_in)/$dr;
				$i_rpos = floor($f_rpos);
				if ($i_rpos >= 0 && $i_rpos < $kmax) {
					$average_v[$i_rpos] += $v_tan;
					$cnt[$i_rpos] ++;
					$particlearea[$i_rpos] += pi*$a*$a;
				} else {
					printf "@ $i $i_rpos   $pos_r\n";
					#exit;
				}
			} 			#$posx[$i] = $x;
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
		
		if ($shear_strain > $shear_strain_steady_state) {
			$cnt_data ++;
			#		exit;
		}

	}
	
}
