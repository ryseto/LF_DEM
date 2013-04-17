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
$out_gaps = "gaps_$name.dat";

open (OUT, "> ${output}");
open (OUT2, "> ${output2}");
open (OUTG, "> ${out_gaps}");
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
&readHeader;
$first=1;
$c_traj=0;
while (1){
	
	&InParticles;
	&InInteractions;
	last unless defined $line;
	&OutYaplotData;
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
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_particle> ;
        ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz,
		$h_xzstress, $c_xzstressXF, $c_xzstressGU, $b_xzstress,  $angle ) = split(/\s+/, $line);
		$radius[$i] = $a;
        $posx[$i] = $x;
        $posy[$i] = $y;
        $posz[$i] = $z;
		$ang[$i] = $angle;
		if ($radius_max < $a){
			$radius_max = $a;
		}
    }
	
	if ($c_traj == 0) {
		$min_dist_origin = 100;
		$min_dist_origin1 = 100;
		$min_dist_origin2 = 100;
		for ($i = 0; $i < $np; $i ++) {
			#$sqdistance = ($posx[$i]*$posx[$i]+$posy[$i]*$posy[$i]+$posz[$i]*$posz[$i]);
			$sqdistance = &calcsqdist($posx[$i],$posy[$i],$posz[$i], 0,0,0);
			if ($sqdistance < $min_dist_origin) {
				$min_dist_origin = $sqdistance;
				$center = $i;
			}
			$sqdistance = &calcsqdist($posx[$i],$posy[$i],$posz[$i], 0,0,5);
			if ($sqdistance < $min_dist_origin1) {
				$min_dist_origin1 = $sqdistance;
				$uppder = $i;
			}
			$sqdistance = &calcsqdist($posx[$i],$posy[$i],$posz[$i], 0,0,-5);
			if ($sqdistance < $min_dist_origin2) {
				$min_dist_origin2 = $sqdistance;
				$lower = $i;
			}
			
		}
	}
	
	$trajx[$c_traj] = $posx[$center];
	$trajy[$c_traj] = $posy[$center];
	$trajz[$c_traj] = $posz[$center];
	$trajx2[$c_traj] = $posx[$uppder];
	$trajy2[$c_traj] = $posy[$uppder];
	$trajz2[$c_traj] = $posz[$uppder];
	$trajx3[$c_traj] = $posx[$lower];
	$trajy3[$c_traj] = $posy[$lower];
	$trajz3[$c_traj] = $posz[$lower];
	$c_traj++;
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
	
	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction> ;
		($i, $j, $f_lub, $fcol, $fc_n, $fc_tx, $fc_ty, $fc_tz,
		$nx, $ny, $nz, $gap, $sxz_lub, $contact) = split(/\s+/, $line);
		# $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
		$int0[$k] = $i;
		$int1[$k] = $j;
		$F_lub[$k] = $f_lub;
		$Sxz_lub[$k] = -($f_lub+$fc_n)*($radius[$i]+$radius[$j])*$nx*$nz;
		$Fc_n[$k] = $fc_n;
		$Ft_t[$k] = sqrt(($fc_tx)**2 + ($fc_ty)**2 + ($fc_tz)**2) ;
		$Fcol[$k] = $fcol;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		$Gap[$k] = $gap;
		printf OUTG "$gap ";
	}
	printf OUTG "\n";
}

sub OutYaplotData{
	if ($first == 0){
		printf OUT "\n";
	} else {
		$first = 0;
	}
	printf OUT "y 7\n";
	printf OUT "r 0.1\n";
    printf OUT "@ 5\n";
	for ($i = 0; $i < $c_traj; $i++){
		$xs = $trajx[$i];
		$ys = $trajy[$i];
		$zs = $trajz[$i];
		$xe = $trajx[$i+1];
		$ye = $trajy[$i+1];
		$ze = $trajz[$i+1];
		if (abs($zs-$ze) < 1
			&& abs($xs-$xe) < 1
			) {
			printf OUT "l $xs $ys $zs $xe $ye $ze\n";
		}
		#		$xs = $trajx[$i];
		#$ys = $trajy[$i];
		#$zs = $trajz[$i];
		#		printf OUT "c $xs $ys $zs\n";
	}
    printf OUT "@ 3\n";
	for ($i = 0; $i < $c_traj; $i++){
		$xs = $trajx2[$i];
		$ys = $trajy2[$i];
		$zs = $trajz2[$i];
		$xe = $trajx2[$i+1];
		$ye = $trajy2[$i+1];
		$ze = $trajz2[$i+1];
		if (abs($zs-$ze) < 1
			&& abs($xs-$xe) < 1
			) {
				printf OUT "l $xs $ys $zs $xe $ye $ze\n";
			}
		#$xs = $trajx2[$i];
		#$ys = $trajy2[$i];
		#$zs = $trajz2[$i];
		#		printf OUT "c $xs $ys $zs\n";
	}
    printf OUT "@ 4\n";
	for ($i = 0; $i < $c_traj; $i++){
		$xs = $trajx3[$i];
		$ys = $trajy3[$i];
		$zs = $trajz3[$i];
		$xe = $trajx3[$i+1];
		$ye = $trajy3[$i+1];
		$ze = $trajz3[$i+1];
		if (abs($zs-$ze) < 1
			&& abs($xs-$xe) < 1
			) {
				printf OUT "l $xs $ys $zs $xe $ye $ze\n";
			}
		#$xs = $trajx3[$i];
		#$ys = $trajy3[$i];
		#$zs = $trajz3[$i];
		#printf OUT "c $xs $ys $zs\n";
	}
	
	printf OUT "y 1\n";
    printf OUT "@ 2\n";
	$r = $radius[0];
	printf OUT "r $r\n";
    for ($i = 0; $i < $np; $i ++){
		if ($i >= 1 && $radius[$i] != $radius[$i-1]){
			$r = $radius[$i];
			printf OUT "r $r\n";
		}
#		if ($i % 100 == 0){
#			$col = $i/100 + 2;
#			printf OUT "@ $col\n";
#		} 
		if ($y_section == 0 ||
			abs($posy[$i]) < $y_section ){
				printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
			}
		
    }
	
	printf OUT "y 2\n";
	printf OUT "@ 0\n";
	for ($k = 0; $k < $num_interaction; $k ++){
		$force = $Fc_n[$k];
        if ($F_lub[$k] < 0) {
			$force += - $F_lub[$k];
		}
		if ($Gap[$k] < 0) {
			$string_width = 0.1;
			&OutString_width($int0[$k],  $int1[$k]);
		}
    }
    printf OUT "y 3\n";
    printf OUT "@ 3\n";
    for ($k=0; $k<$num_interaction; $k++){
		#$force = $F_lub[$k] + $Fc_n[$k] + $Fcol[$k];
		#$force = $Fcol[$k];
		#$force = $Fc_n[$k];
		$force = $F_lub[$k];
        if ($force < -5){
			$string_width = (-${force_factor})*${force};
			&OutString_width($int0[$k], $int1[$k]);
		}
    }
	printf OUT "y 4\n";
   printf OUT "@ 4\n";
   for ($k = 0; $k < $num_interaction; $k ++){
	   if ($Fc_n[$k] > 1){
		   $force = $Fc_n[$k] + $Fcol[$k];
	   } else {
		   $force = $F_lub[$k] + $Fcol[$k];
	   }
	   if ($force > 5){
		   $string_width = ${force_factor}*${force};
		   &OutString_width($int0[$k], $int1[$k]);
	   }
    }
	printf OUT "y 5\n";
	printf OUT "@ 5\n";
	for ($k = 0; $k < $num_interaction; $k ++){
		$radius = $ContVelo[$k]/10;
		printf OUT "r $radius\n";
		&OutCircle_middle($int0[$k],  $int1[$k]);
    }
	
	if ($Ly == 0){
		printf OUT "y 6\n";
		printf OUT "@ 0\n";
		for ($i = 0; $i < $np; $i ++){
			&OutCross($i);
		}
	}
	
	#	$maxS=0;
	#for ($k = 0; $k < $num_interaction; $k ++){
	#		if ($maxS < $Sxz_lub[$k]){
	#			$maxS = $Sxz_lub[$k];
	#		}
	#}
	
	#	printf OUT "y 4\n";
	#    printf OUT "@ 5\n";
	#    for ($k = 0; $k < $num_interaction; $k ++){
	#		$string_with = $Sxz_lub[$k]/$maxS;
	#		printf OUT "r ${string_with}\n";
	#		&OutCircle_middle($int0[$k],  $int1[$k]);
	#	}
	#	$zpos = $Lz / 2 + 1;
	#	printf OUT sprintf("t 0 0 %3.2f %2.4f\n", $zpos, $maxS);
	#
	#	printf OUT "y 5\n";
	#	printf OUT "@ 2\n";
	#printf OUT "r 0.3\n";
    #for ($k = 0; $k < $num_interaction; $k ++){
	#	if ($Gap[$k] < 0.01 ){
	#			&OutString($int0[$k],  $int1[$k]);
	#			&OutNvec($k);
	#		}
	#    }
	printf OUT2 "\n";
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
	$sq_dist = ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2;
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














