#!/usr/bin/perl

# Usage:
# $ averageConfigdata.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;
$particle_data = $ARGV[0];
$num_average = 1;
$magnetic_simulation = 0;

if (@ARGV >= 2) {
	$num_average = $ARGV[1];
}
$output_skip = 1;
if (@ARGV >= 3) {
	$output_skip = $ARGV[2];
}
$output_cnt = $output_skip-1;

if (@ARGV >= 4) {
	# unit of time
	$timeunit = $ARGV[3];
}

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

if ($num_average == 1) {
	$outputfilename = "config_$name.dat";
} else {
	$outputfilename = "config_${name}_av${num_average}.dat";
}

&importDataFile;

open (OUT, "> ${outputfilename}");
open (IN_particle, "< ${particle_data}");

&readHeader;

$lxhalf = $Lx/2;
$lyhalf = $Ly/2;
$lzhalf = $Lz/2;
$cnt = 0;
$frame = 0;
$avecount = 0;

my $magicangle = atan(0.5*(sqrt(5)-1));
my $cos_ma = cos($magicangle);
my $sin_ma = sin($magicangle);


for ($i=0; $i < $np; $i++) {
	$averagex[$i] = 0;
	$averagez[$i] = 0;
}

$periodicboundary = 1;
printf OUT "Lx $Lx Ly $Lz Periodic $periodicboundary\n";

while (1){
	&InParticles;
	if ($time > 0) {
		if ($num_average == 1) {
			&directOutput;
		} else {
			&averageOutput;
		}
	}
	if ($time > 300) {
		
		printf "exit 300 $time \n";
		exit;
	}
	printf "$line\n";
	last unless defined $line;
}

close (OUT);
close (IN_particle);

sub readHeader {
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $flowtype) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $datatype) = split(/\s+/, $line);
	#	for ($i = 0; $i < 16; $i++) {
	#		$line = <IN_particle>;
	#	}
	
	if ($magnetic_simulation == 1) {
		#	$iskip = 9;
		$iskip = 16;
		
	} else {
		$iskip = 7;
		
	}
	for ($i = 0; $i < $iskip; $i++) {
		$line = <IN_particle>;
	}
	
}

sub importDataFile {
	open (IN_data, "< data_${name}.dat");
	while (1) {
		$line = <IN_data>;
		($d1, $d2, $d3) = split(/\s+/, $line);
		last unless defined $line;
		if ($d1 != "#") {
			$time = $d1;
		}
	}
	close (IN_data);
}

sub InParticles {
	$radius_max = 0;
	$line = <IN_particle>;
	if (defined $line) {
		#($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress, $time) = split(/\s+/, $line);
		($buf, $shear_strain, $shear_disp, $shear_rate, $shear_stress, $time, $timeper, $epsilondot) = split(/\s+/, $line);
		if ($time > 0) {
			if ($magnetic_simulation == 1){
				$num_mag = 0;
				$num_nonmag = 0;
				if (1) {
					for ($i = 0; $i < $np; $i ++){
						$line = <IN_particle> ;
						#($ip, $a, $x, $y, $z, $mx, $my, $mz, $ms) = split(/\s+/, $line);
						($ip, $a, $x, $y, $z, $vx, $vy, $vz, $mx, $my, $mz, $ms) = split(/\s+/, $line);
						$posx[$i] = $x + $lxhalf;
						$posz[$i] = $z + $lxhalf;
						$radius[$i] = $a;
						$magsus[$i] = $ms;
						if ($ms > 0) {
							$num_mag ++;
						} else {
							$num_nonmag ++;
						}
					}
				} else {
					for ($i = 0; $i < $np; $i ++){
						$line = <IN_particle> ;
					}
				}
			} else {
				for ($i = 0; $i < $np; $i ++){
					$line = <IN_particle> ;
					#($ip, $a, $x, $y, $z, $mx, $my, $mz, $ms) = split(/\s+/, $line);
					($ip, $a, $x, $y, $z, $vx, $vy, $vz, $avx, $avy, $avz, $ang) = split(/\s+/, $line);
					$posx[$i] = $x + $lxhalf;
					$posz[$i] = $z + $lxhalf;
					$radius[$i] = $a;
				}
			}
		} else {
			for ($i = 0; $i < $np; $i ++){
				$line = <IN_particle> ;
			}
		}
	}
}

sub directOutput {
	my $timeinmin_int = int (100000*((${time}*${timeunit})/60));
	my $timeinmin = $timeinmin_int/100000;
	#	printf "$time $timeunit $timeinmin\n";
	#printf "$cnt $time $timeinmin $time\n";
	
	$cnt ++;
	
	$output_cnt++;
	$output = 0;
	if ($output_cnt == $output_skip) {
		$output_cnt = 0;
		$output = 1;
	}
	if ($output) {
		
		$frame ++;
		
		$epsilondot = 0.5;
		$ax = exp( $epsilondot*($timeper))*$cos_ma*$cos_ma+exp(-$epsilondot*($timeper))*$sin_ma*$sin_ma;
		$az = exp(-$epsilondot*($timeper))*$cos_ma*$sin_ma-exp( $epsilondot*($timeper))*$sin_ma*$cos_ma;
		$bx = exp(-$epsilondot*($timeper))*$cos_ma*$sin_ma-exp( $epsilondot*($timeper))*$sin_ma*$cos_ma;
		$bz = exp(-$epsilondot*($timeper))*$cos_ma*$cos_ma+exp( $epsilondot*($timeper))*$sin_ma*$sin_ma;
		
		
		if ($magnetic_simulation == 1) {
			printf "$cnt $timeinmin\n";
			printf OUT "frame $frame\n";
			printf OUT "time $timeinmin\n";
			printf OUT "number_of_magnetic_particles $num_mag\n";
			for ($i=0; $i < $num_mag; $i++) {
				my $xx_int = int (100000*$posx[$i]);
				my $xx = $xx_int/100000;
				my $zz_int = int (100000*$posz[$i]);
				my $zz = $zz_int/100000;
				$a = $radius[$i];
				$ms = $magsus[$i];
				printf OUT "$xx $zz $a $ms\n";
			}
			printf OUT "number_of_nonmagnetic_particles $num_nonmag\n";
			for ($i=$num_mag; $i < $np; $i++) {
				$xx = $posx[$i];
				$zz = $posz[$i];
				$a = $radius[$i];
				$ms = $magsus[$i];
				printf OUT "$xx $zz $a $ms\n";
			}
		} else {
			if (0) {
				printf "$cnt $time\n";
				printf OUT "frame $frame\n";
				printf OUT "strain $shear_strain \n";
				printf OUT "disp $shear_disp\n";
				printf OUT "np $np\n";
				for ($i=0; $i < $np; $i++) {
					my $xx_int = int (100000*$posx[$i]);
					my $xx = $xx_int/100000;
					my $zz_int = int (100000*$posz[$i]);
					my $zz = $zz_int/100000;
					$a = $radius[$i];
					$ms = $magsus[$i];
					printf OUT "$xx $zz $a \n";
				}
			} else {
				&trimmingKR;
				printf "$cnt $time\n";
				printf OUT "frame $frame\n";
				printf OUT "strain $shear_strain \n";
				printf OUT "disp $shear_disp\n";
				printf OUT "np $np_kr\n";
				for ($i=0; $i < $np_kr; $i++) {
					#	my $xx_int = int (100000*$posx2[$i]);
					#my $xx = $xx_int/100000;
					#my $zz_int = int (100000*$posz2[$i]);
					#my $zz = $zz_int/100000;
					my $xx = $posx2[$i];
					my $zz = $posz2[$i];
					$a = $radius2[$i];
					$ms = $magsus2[$i];
					printf OUT "$xx $zz $a \n";
				}
			}
		}
	}
}

sub trimmingKR {
	my $j = 0;
	my $axisrotation = 1;
	my $scale = 0.5;
	for ($ii = -3; $ii <= 3; $ii++) {
		for ($jj = -3; $jj <= 3; $jj++) {
			for ($i = 0; $i < $np; $i++) {
				$rr = $yap_radius*$radius[$i];
				
				$pd_xshift = $Lx*$ax*($ii)+$Lz*$bx*($jj);
				$pd_zshift = $Lx*$az*($ii)+$Lz*$bz*($jj);
				$xx = $posx[$i] + $pd_xshift;
				$yy = $posy[$i] - $Ly/2;
				$zz = $posz[$i] + $pd_zshift;
				if ($axisrotation) {
					$xx2 = cos($magicangle)*$xx - sin($magicangle)*$zz;
					$zz2 = sin($magicangle)*$xx + cos($magicangle)*$zz;
				} else {
					$xx2 = $xx;
					$zz2 = $zz;
				}
				#if ($xx2> 0 && $xx2 < $Lx &&  $zz2 > 0 && $zz2 < $Lz){
				if (abs($xx2)<$scale*$Lx && abs($zz2) < $scale*$Lz) {
					$posx2[$j] = $xx2 +$Lx/2;
					$posy2[$j] = $posy[$i];
					$posz2[$j] = $zz2 +$Lz/2;
					$radius2[$j] = $radius[$i];
					$ang2[$j] = $ang[$i];
					$j++;
					#printf OUT "c $xx2 $yy $zz2 \n";
				}
				# last LOOP2 if ($xx2<$Lx && $xx2 > 0 && $zz2 > 0 && $zz2 < $Lz);
			}
		}
	}
	
	$np_kr = $j;
	
	
}



sub averageOutput {
	if ($avecount == 0){
		for ($i=0; $i < $np; $i++) {
			$posx0[$i] = $posx[$i];
			$posz0[$i] = $posz[$i];
		}
	}
	
	for ($i=0; $i < $np; $i++) {
		$xpd[$i] = $posx[$i];
		$zpd[$i] = $posz[$i];
		if (abs($xpd[$i] - $posx0[$i]) > 0.5*$Lx) {
			if ($xpd[$i] > $posx0[$i]){
				$xpd[$i] -= $Lx;
			} else {
				$xpd[$i] += $Lx;
			}
		}
		if (abs($zpd[$i] - $posz0[$i]) > 0.5*$Lz) {
			if ($zpd[$i] > $posz0[$i]){
				$zpd[$i] -= $Lz;
			} else {
				$zpd[$i] += $Lz;
			}
		}
		$averagex[$i] += $xpd[$i];
		$averagez[$i] += $zpd[$i];
	}
	$avecount ++;
	if ($avecount == $num_average) {
		
		my $timeinmin_int = int (100000*((${time}*${timeunit})/60));
		my $timeinmin = $timeinmin_int/100000;
		
		
		printf "$cnt $timeinmin\n";
		printf OUT "frame $cnt\n";
		printf OUT "time $timeinmin\n";
		printf OUT "number_of_magnetic_particles $num_mag\n";
		for ($i=0; $i < $num_mag; $i++) {
			$xx = $averagex[$i]/$num_average;
			$zz = $averagez[$i]/$num_average;
			if ($xx < 0) {
				$xx += $Lx;
			} elsif ($xx > $Lx) {
				$xx -= $Lx;
			}
			if ($zz < 0) {
				$zz += $Lz;
			} elsif ($xx > $Lz) {
				$zz -= $Lz;
			}
			my $xx_int = int (100000*$xx);
			my $xx = $xx_int/100000;
			my $zz_int = int (100000*$zz);
			my $zz = $zz_int/100000;
			$a = $radius[$i];
			$ms = $magsus[$i];
			$averagex[$i] = 0;
			$averagez[$i] = 0;
			
			printf OUT "$xx $zz $a $ms\n";
		}
		printf OUT "number_of_nonmagnetic_particles $num_nonmag\n";
		for ($i=$num_mag; $i < $np; $i++) {
			$xx = $posx[$i];
			$zz = $posz[$i];
			
			$xx = $averagex[$i]/$num_average;
			$zz = $averagez[$i]/$num_average;
			if ($xx < 0) {
				$xx += $Lx;
			} elsif ($xx > $Lx) {
				$xx -= $Lx;
			}
			if ($zz < 0) {
				$zz += $Lz;
			} elsif ($xx > $Lz) {
				$zz -= $Lz;
			}
			my $xx_int = int (100000*$xx);
			my $xx = $xx_int/100000;
			my $zz_int = int (100000*$zz);
			my $zz = $zz_int/100000;
			$a = $radius[$i];
			$ms = $magsus[$i];
			$averagex[$i] = 0;
			$averagez[$i] = 0;
			
			
			$a = $radius[$i];
			$ms = $magsus[$i];
			printf OUT "$xx $zz $a $ms\n";
		}
		
		$cnt ++;
		$avecount = 0;
	}
}
