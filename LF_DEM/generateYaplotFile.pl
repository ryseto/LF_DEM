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
    ($buf, $shear_rate) = split(/\s+/, $line);
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_particle> ;
        ($ip, $a, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz, $angle ) = split(/\s+/, $line);
		$radius[$i] = $a;
        $posx[$i] = $x;
        $posy[$i] = $y;
        $posz[$i] = $z;
		$ang[$i] = $angle;
		if ($radius_max < $a){
			$radius_max = $a
		}
    }
}

sub InInteractions {
	$line = <IN_interaction>;
	($buf, $shear_rate, $num_interaction) = split(/\s+/, $line);
	if ($buf != "#"){
		exit;
	}
	for ($k = 0; $k < $num_interaction; $k ++){
		$line = <IN_interaction> ;
		($i, $j, $f_lub, $fc_n, $fc_tx, $fc_ty, $fc_tz,
		$nx, $ny, $nz, $gap, $sxz_lub, $neartime, $fric_st) = split(/\s+/, $line);
		$int0[$k] = $i;
		$int1[$k] = $j;
		$F_lub[$k] = $f_lub;
		$Sxz_lub[$k] = -($f_lub+$fc_n)*($radius[$i]+$radius[$j])*$nx*$nz;
		$Fc_n[$k] = $fc_n;
		$Ft_t[$k] = ($fc_tx)**2 + ($fc_ty)**2 + ($fc_tz)**2 ;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		
		$NearTime[$k] = $neartime;
		$Gap[$k] = $gap;
	}
}

sub OutYaplotData{
	if ($first == 0){
		printf OUT "\n";
	} else {
		$first = 0;
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
		if ($y_section == 0 ||
			abs($posy[$i]) < $y_section ){
				printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
			}
    }
	
    printf OUT "y 2\n";
    printf OUT "@ 2\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        if ($Fc_n[$k] > 0){
            $string_width = $force_factor*($Fc_n[$k]);
			&OutString_width($int0[$k],  $int1[$k]);
        }
    }
    printf OUT "y 3\n";
    printf OUT "@ 3\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        $force = $F_lub[$k] + $Fc_n[$k];
        if ($force < 0){
			$string_width = -${force_factor}*${force};
			&OutString_width($int0[$k], $int1[$k]);
		}
    }
    printf OUT "@ 4\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        $force = $F_lub[$k];
		if ($force > 0){
            $string_width = ${force_factor}*${force};
			&OutString_width($int0[$k], $int1[$k]);
		}
    }

	printf OUT "y 6\n";
	printf OUT "@ 0\n";
	for ($i = 0; $i < $np; $i ++){
		&OutCross($i);
		
	}
	
	$maxS=0;
	for ($k = 0; $k < $num_interaction; $k ++){
		if ($maxS < $Sxz_lub[$k]){
			$maxS = $Sxz_lub[$k];
		}
	}
	
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
	
	if (abs($xi-$xj) < $radius_max*5
		&&  abs($yi-$yj) < $radius_max*5
		&&  abs($zi-$zj) < $radius_max*5){
			if ( $y_section == 0
				|| abs($yi) < $y_section
				|| abs($yj) < $y_section){
					printf OUT "r ${string_width}\n";
					printf OUT "s $xi $yi $zi $xj $yj $zj\n";
				}
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
    $xc = ($posx[$i] + $posx[$j])/2;
    $yc = ($posy[$i] + $posy[$j])/2;
    $zc = ($posz[$i] + $posz[$j])/2;
	
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














