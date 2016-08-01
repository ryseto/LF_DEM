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
my $kmax = $ARGV[1];
my $d_strain = $ARGV[2];
my $strain_end = $ARGV[3];

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$j = index($name, '__', 1);
$initconfig = substr($name, 0, $j+1);

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


#$output_pos = "pos_$name.dat";
printf "$output\n";

#open (OUTpos, "> ${output_pos}");
open (IN_particle, "< ${particle_data}");

&readHeader;

$first = 1;
$output = 1;


$r_in = $radius_in;
$r_out = $radius_out;
$rmargin = $margin*($r_out - $r_in);
$r_range_in = $r_in + $rmargin/2;
$r_range_out = $r_out - $rmargin/2;
$rdiff = ${r_range_out}-${r_range_in};
$dr = $rdiff/$kmax;
$v_in = $radius_in;

printf "$r_in $r_out\n";
printf "$r_range_in $r_range_out\n";
printf "dr = $dr \n";

## $stress_in = $average_stress_rt[0]/$cnt[0];


&resetAverage;

$shear_strain_next = $d_strain;
$cnt = 0;
while (1) {
	&InParticles;
	last unless defined $line;
	if ($shear_strain > $shear_strain_next) {
		$output = "CWGdata_${name}_${cnt}.dat";
		open (OUT, "> ${output}");
		&outputAverage;
		close (OUT);
		&resetAverage;
		$shear_strain_next += $d_strain;
		${cnt} += 1;
	}
	if ($shear_strain > ${strain_end}) {
		last;
	}

}
#
#for ($k = 0; $k < $kmax; $k++) {
#	if ($cnt[$k] != 0) {
#		$ave_v_tan[$k] = $average_v[$k]/$cnt[$k];
#		$ave_stress_rr[$k] = $average_stress_rr[$k]/$cnt[$k];
#	}
#}
#close (OUTpos);
close (IN_particle);
##################################################################

sub resetAverage {
	$cnt_data = 0;

	for ($k = 0; $k < $kmax; $k++) {
		$average[$k] = 0;
		$radialposition[$k] = 0;
		$cnt[$k] = 0;
		$particlearea[$k] = 0;
		
		$average_stress_rr[$k] = 0;
		$average_stress_tt[$k] = 0;
		$average_stress_rt[$k] = 0;
		$average_stress_N1[$k] = 0;
		$average_stress_pp[$k] = 0;
		
	}
}

sub outputAverage {
	printf OUT "# strain = $shear_strain\n";
	printf OUT "# R_in = ${r_in}\n";
	printf OUT "# R_out = ${r_out}\n";
	printf OUT "# Rrange0 = $r_range_in\n";
	printf OUT "# Rrange1 = $r_range_out\n";
	printf OUT "# kmax = $kmax\n";
	printf OUT "#1 radial position\n";
	printf OUT "#2 velocity_theta \n";
	printf OUT "#3 area fraction\n";
	printf OUT "#4 stress tensor r-r\n";
	printf OUT "#5 stress tensor theta-theta\n";
	printf OUT "#6 stress tensor r-theta\n";
	printf OUT "#7 N1 = s_tt - s_rr\n";
	printf OUT "#8 P = -(1/2)*(s_tt + s_rr)\n";

	for ($k = 0; $k < $kmax; $k++) {
		if ($cnt[$k] != 0) {
			$ave_v_tan = $average_v[$k]/$cnt[$k];
			$r = $r_range_in + $dr*$k;
			$rmid = $r + 0.5*$dr;
			$rn = $r + $dr;
			$area = pi*($rn*$rn - $r*$r);
			$areafraction = ($particlearea[$k]/$cnt_data)/$area;
			$ave_stress_rr = $areafraction*$average_stress_rr[$k]/$cnt[$k];
			$ave_stress_tt = $areafraction*$average_stress_tt[$k]/$cnt[$k];
			$ave_stress_rt = $areafraction*$average_stress_rt[$k]/$cnt[$k];
			$ave_stress_N1 = $areafraction*$average_stress_N1[$k]/$cnt[$k];
			$ave_stress_PP = $areafraction*$average_stress_pp[$k]/$cnt[$k];
			printf OUT "$rmid $ave_v_tan $areafraction $ave_stress_rr $ave_stress_tt $ave_stress_rt $ave_stress_N1 $ave_stress_PP \n";
		}
	}
}



sub readHeader {
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	for ($i = 0; $i<11; $i++) {
		$line = <IN_particle>;
	}
	for ($i = 0; $i<16; $i++) {
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
			$stress_rr, $stress_tt, $stress_rt, $angle) = split(/\s+/, $line);
			$ang[$i] = $angle;
			$radius[$i] = $a;
			if ($shear_strain > $shear_strain_steady_state && $i < $np_mov) {
				$pos_r2 = $x*$x + $z*$z;
				$pos_r = sqrt($pos_r2);
				if ( $pos_r > $r_range_in && $pos_r < $r_range_out ) {
					$v_tan = ((-$vx*$z + $vz*$x)/$pos_r);
					$f_rpos = ($pos_r - $r_range_in)/$dr;
					$i_rpos = floor($f_rpos);
					if ($i_rpos >= 0 && $i_rpos < $kmax) {
						$average_v[$i_rpos] += $v_tan;
						$average_stress_rr[$i_rpos] += $stress_rr;
						$average_stress_tt[$i_rpos] += $stress_tt;
						$average_stress_rt[$i_rpos] += $stress_rt;
						$average_stress_pp[$i_rpos] += 0.5*(-$stress_rr-$stress_tt);
						$average_stress_N1[$i_rpos] += ($stress_tt-$stress_rr);
						$particlearea[$i_rpos] += pi*$a*$a;
						$cnt[$i_rpos] ++;
					} else {
						printf "@ $i $i_rpos   $pos_r\n";
						exit;
					}
				}
			} 			#$posx[$i] = $x;
			$posx[$i] = $x;
			$posz[$i] = $z;
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
