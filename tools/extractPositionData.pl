#!/usr/bin/perl
# Usage:
# $ extractPositionData.pl par_[...].dat [strain]
use Math::Trig;
use IO::Handle;
use Getopt::Long;

my $particle_data = $ARGV[0];
my $config_dir = $ARGV[1];
my $yap_radius = 1;
my $output_interval = 1;
my $xz_shift = 0;
my $axis = 0;
my $reversibility_test = 0;
my $monodisperse = 0;

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$i2 = index($particle_data, '__', $i)+1;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$configname = substr($particle_data, $i, $i2-$i);
open (IN_configdata, "< $config_dir/${configname}.dat");
$configline1 = <IN_configdata>;
$configline2 = <IN_configdata>;
close (IN_configdata);
printf "$configline1";
printf "$configline2";
open (IN_particle, "< ${particle_data}");
&readHeader;
$outputnum = 1;
$end_of_file = 0;
while (1) {
	&InParticles;
	if ($end_of_file eq 1) {
		last;
	}
}
$outputfilename = "${configname}rlx.dat";
open (OUT, "> ${outputfilename}");
($buf, $np1, $np2, $vf, $lx, $ly, $lz, $vf1, $vf2, $dispx0, $dispy0) = split(/\s+/, $configline2);
printf OUT "# np1 np2 vf lx ly lz vf1 vf2 dispx dispy\n";
printf OUT "# $np1 $np2 $vf $lx $ly $lz $vf1 $vf2 $dispx 0\n";
for ($i = 0; $i < $np; $i++) {
	$xx = $posx[$i] + $Lx/2;
	$yy = $posy[$i] + $Ly/2;
	$zz = $posz[$i] + $Lz/2;
	$rr = $radius[$i];
	printf OUT "$xx $yy $zz $rr\n";
	#	printf "$xx $yy $zz $rr\n";
}
close (OUT);
$outputnum ++;
close (IN_particle);
close (IN_interaction);
##################################################################
sub keepInitialConfig {
	for ($i = 0; $i < $np; $i ++){
		$posx_init[$i] = $posx[$i];
		$posy_init[$i] = $posy[$i];
		$posz_init[$i] = $posz[$i];
		$ang_init[$i] = $ang[$i];
		$radius_init[$i] = $radius[$i];
	}
}

sub readHeader {
	$line = <IN_particle>;
	$line = <IN_particle>; ($buf, $buf, $np) = split(/\s+/, $line);
	printf "np = $np\n";
	$line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
	printf "VF = $VF\n";

	$line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
	printf "Lx Ly Lz = $Lx $Ly $Lz\n";
	$line = <IN_particle>; ($buf, $buf, $flwtyp) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $buf, $dataunit) = split(/\s+/, $line);
	
	if ($Ly==0) {
		$number_of_header = 8;
	} else {
		$number_of_header = 7;
	}
	for ($i = 0; $i<$number_of_header; $i++) {
		$line = <IN_particle>;
		printf "$line";
	}
	#	$xo = $Lx/2;
	#	$yo = $Ly/2;
	#	$zo = $Lz/2;
	$xo = 0;
	$yo = 0;
	$zo = 0;
}

sub InParticles {
	$radius_max = 0;
	##  Snapshot Header
	$j = 0;
	while (1) {
		$line = <IN_particle>;
		if (!defined($line)) {
			$end_of_file = 1;
			last;
		}
		#		last unless defined $line;
		($buf, $val) = split(" : ", $line);
		($buf1) = split(/\s+/, $buf);
		if ($buf1 ne '#') {
			last;
		} else {
			$val =~ s/(\n|\r)//g;
			$ssHeader[$j++] = $val;
		}
	}
	if ($end_of_file eq 0) {
		$shear_strain = $ssHeader[0];
		$dispx = $ssHeader[1];
		$shear_rate = $ssHeader[2];
		$target_stress = $ssHeader[3];
		for ($i = 0; $i < $np; $i ++){
			if ($i > 0) {
				$line = <IN_particle>;
			}
			#($ip, $a, $x, $y, $z) = split(/\s+/, $line);
			#			($ip, $a, $x, $y, $z) = split(/\s+/, $line);
			($ip, $a, $x, $z) = split(/\s+/, $line);
			$y = 0;
			#($ip, $a, $x, $z, $vx, $vz, $vy, $ox, $oz, $oy, $angle) = split(/\s+/, $line);
			$ang[$i] = $angle;
			$radius[$i] = $a;
			$posx[$i] = $x-$xo;
			#	$posy[$i] = $y-$yo;
			$posz[$i] = $z-$zo;

		}
	}

}
