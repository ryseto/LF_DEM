#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;

$force_factor = 0.01;
$y_section = 0;

$particle_data = $ARGV[0];
if ($#ARGV >= 1){
	$force_factor = $ARGV[1];
}
if ($#ARGV == 2){
	$y_section = $ARGV[2];
	printf "section $y_section\n";
}

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);

$interaction_data = "int_${name}.dat";
printf "$interaction_data\n";
$output = "fc_$name.yap";
open (OUT, "> ${output}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
$first=1;
$c_traj=0;
$num = 0;
while (1){
	&InParticles;
	&InInteractions;
	last unless defined $line;

	&OutYaplotData;
	$num ++;
	printf "$shear_rate\n";
}
close (OUT);

sub readHeader{
	$line = <IN_particle>; ($buf, $np) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $VF) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $Lx) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $Ly) = split(/\s+/, $line);
	$line = <IN_particle>; ($buf, $Lz) = split(/\s+/, $line);
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
		#		$ang[$i] = $angle;
		#		if ($radius_max < $a){
		#	$radius_max = $a;
		#}
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
	($buf, $shear_rate, $num_interaction) = split(/\s+/, $line);
	if ($buf != "#"){
		exit;
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
	
	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction> ;
		($i, $j, $contact, $nx, $ny, $nz,
		$gap, $f_lub, $fc_n, $fc_tan, $fcol,
		$sxz_cont_xF, $n1_cont_xF, $n2_cont_xF) = split(/\s+/, $line);
		
		
		if ($num==$num_mathm){
			printf OUTMI "$line";
		}
		# $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
		$int0[$k] = $i;
		$int1[$k] = $j;
		$force[$k] = $f_lub + $fc_n + $fc_tan + $fcol;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		$Gap[$k] = $gap;
		printf OUTG "$gap ";
	}
	printf OUTG "\n";
}

sub analzeForceChain {
	$maxf = 0;
	$k_max = 0;
	for ($k=0; $k<$num_interaction; $k++) {
		if ($force[$k] > $maxf) {
			$maxf = $force[$k];
			$k_max = $k;
		}
	}
	$i0 = $int0[$k_max];
	$i1 = $int1[$k_max];
	
	printf OUT ("%3.5f %3.5f %3.5f ", posx[$i0], posy[$i0], posz[$i0]);

	
	
	
	
}



















