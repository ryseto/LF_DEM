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
$output = "y_$name.yap";
$output2 = "nvec_$name.dat";

open (OUT, "> ${output}");
open (OUT2, "> ${output2}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
$first=1;
$c_traj=0;
while (1){
	&InParticles;
	last unless defined $line;
	&OutYaplotData;
	#	printf "$shear_rate\n";
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
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_particle> ;
        ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
		$h_xzstress, $c_xzstressXF, $c_xzstressGU, $b_xzstress,  $angle ) = split(/\s+/, $line);
		$radius[$i] = $a;
        $posx[$i][$c_traj] = $x;
        $posy[$i][$c_traj] = $y;
        $posz[$i][$c_traj] = $z;
		if ($radius_max < $a){
			$radius_max = $a;
		}
    }
	$c_traj++;
}

sub calcsqdist {
    ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
	$dist = ($x1-$x2)*($x1-$x2);
	$dist += ($y1-$y2)*($y1-$y2);
	$dist += ($z1-$z2)*($z1-$z2);
	return $dist;
}

sub OutYaplotData{
	if ($first == 0){
		printf OUT "\n";
	} else {
		$first = 0;
	}
	printf OUT "y 7\n";
	
	printf OUT "r 0.1";

	for ($j = 0; $j < $np; $j++) {
		$col = $j % 6 + 2;
		printf OUT "@ $col\n";
		for ($i = 0; $i < $c_traj-3; $i += 3){
			$xs = $posx[$j][$i];
			$ys = $posy[$j][$i];
			$zs = $posz[$j][$i];
			$xe = $posx[$j][$i+3];
			$ye = $posy[$j][$i+3];
			$ze = $posz[$j][$i+3];
			if (abs($zs-$ze) < 2
				&& abs($xs-$xe) < 2
				&& abs($ys-$ye) < 2
				) {
					printf OUT "l $xs $ys $zs $xe $ye $ze\n";
				}
			#		$xs = $trajx[$i];
			#		$ys = $trajy[$i];
			#		$zs = $trajz[$i];
			#		printf OUT "c $xs $ys $zs\n";
		}
	}
}

sub OutBoundaryBox{
	$x0 = -$Lx/2;
	$x1 = -$Lx/2 + $shear_disp / 2;
	$z1 = $Lz/2;
	$x2 = $Lx/2;
	$z2 = 0;
	$x3 = $Lx/2 - $shear_disp / 2;
	$z3 = -$Lz/2;
	$lx2 = $Lx/2;
	
	printf OUT "y 7\n";
	printf OUT "@ 6\n";
	printf OUT "l -$lx2 0 0 $lx2 0 0\n";
	printf OUT "l $x0 0.01 0 $x1 0.01 $z1\n";
	printf OUT "l $x2 0.01 $z2 $x3 0.01 $z3\n";
	
}


sub OutNvec {
    ($k) = @_;
	$nx = $nrvec_x[$k];
	$ny = $nrvec_y[$k];
	$nz = $nrvec_z[$k];
	printf OUT2 "$nx $ny $nz\n";
	
}

sub OutString_width {
    ($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i] - 0.01;
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j] - 0.01;
    $zj = $posz[$j];
	$sq_dist = ($xi-$xj)**2 +	($yi-$yj)**2 + ($zi-$zj)**2;
	if (sqrt($sq_dist) < $radius[$i] + $radius[$j]+1){
		printf OUT "r ${string_width}\n";
		printf OUT "s $xi $yi $zi $xj $yj $zj\n";
	}
	
}


sub OutString {
    ($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i];
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j];
    $zj = $posz[$j];
	
	if (abs($xi-$xj) < $radius_max*5
		&&  abs($yi-$yj) < $radius_max*5
		&&  abs($zi-$zj) < $radius_max*5){
			if ( $y_section == 0
				|| abs($yi) < $y_section
				|| abs($yj) < $y_section){
					printf OUT "s $xi $yi $zi $xj $yj $zj\n";
				}
		}
}

sub OutCircle_middle {
    ($i, $j) = @_;
	$xi = $posx[$i];
    $yi = $posy[$i];
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j];
    $zj = $posz[$j];
	$ai = $radius[$i];
	$aj = $radius[$j];
    $xc = $posx[$i]*$aj/($ai+$aj) + $posx[$j]*$ai/($ai+$aj);
    $yc = $posy[$i]*$aj/($ai+$aj) + $posy[$j]*$ai/($ai+$aj) -0.01;
    $zc = $posz[$i]*$aj/($ai+$aj) + $posz[$j]*$ai/($ai+$aj);
	
	if (abs($xi-$xj) < $radius_max*5
		&&  abs($yi-$yj) < $radius_max*5
		&&  abs($zi-$zj) < $radius_max*5){
			printf OUT "c $xc $yc $zc\n";
		}
}

sub OutCross {
    ($i) = @_;
	$a = $radius[$i];
    $xi = $posx[$i];
    $yi = $posy[$i] - 0.01;
    $zi = $posz[$i];
    $angle = $ang[$i];
	$ux =  $a*cos($angle);
	$uz = -$a*sin($angle);
	$xa = $xi - $ux;
	$ya = $yi - 0.01;
	$za = $zi - $uz;
	$xb = $xi + $ux;
	$yb = $yi - 0.01;
	$zb = $zi + $uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
	$ux =  $a*sin($angle);
	$uz =  $a*cos($angle);
	$xa = $xi - $ux;
	$ya = $yi - 0.01;
	$za = $zi - $uz;
	$xb = $xi + $ux;
	$yb = $yi - 0.01;
	$zb = $zi + $uz;
	printf OUT "l $xa $ya $za $xb $yb $zb\n";
	
}














