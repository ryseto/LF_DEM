#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
#use strict;
#use warnings;

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
$counter_fc = 0;
$istart0 = -1;
$istart1 = -1;
$k_max = -1;
$newchain = 1;
while (1){
	&InParticles;
	&InInteractions;
	last unless defined $line;
	for ($i=0; $i<$np; $i++) {
		@{int_in_p[$i]} = ();
	}
	for ($k=0; $k<$num_interaction; $k++) {
		$ii0 = $int0[$k];
		$ii1 = $int1[$k];
		# Particles memorize connecting bonds.
		# One bond is recorded in the connecting two particles.
		# Therefore, one way is rejected due to the comparison of thier forces,
		# another way should be accepted.
		if ($force[$k] > 1){
			push(@{int_in_p->[$ii0]}, $k);
			push(@{int_in_p->[$ii1]}, $k);
		}
	}
	$m = 0;
	FCLOOP: {
		do {
			&analyzeForceChain;
			$m ++;
			if ($k_max == -1){
				last FCLOOP;
			}
		} while ( $num_k >= 10);
	}

	printf OUT "@ 0\n";
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2,-$Ly/2,-$Lz/2, $Lx/2,-$Ly/2,-$Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2, $Ly/2,-$Lz/2, $Lx/2, $Ly/2,-$Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2,-$Ly/2, $Lz/2, $Lx/2,-$Ly/2, $Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2, $Ly/2, $Lz/2, $Lx/2, $Ly/2, $Lz/2);
	
	
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n",  $Lx/2,-$Ly/2, $Lz/2, $Lx/2, $Ly/2, $Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2,-$Ly/2, $Lz/2,-$Lx/2, $Ly/2, $Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n",  $Lx/2,-$Ly/2,-$Lz/2, $Lx/2, $Ly/2,-$Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2,-$Ly/2,-$Lz/2,-$Lx/2, $Ly/2,-$Lz/2);
	
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n",  $Lx/2, $Ly/2, $Lz/2, $Lx/2, $Ly/2,-$Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2, $Ly/2, $Lz/2,-$Lx/2, $Ly/2,-$Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n", -$Lx/2,-$Ly/2, $Lz/2,-$Lx/2,-$Ly/2,-$Lz/2);
	printf OUT ("l %3.3f %3.3f %3.3f  %3.3f %3.3f %3.3f\n",  $Lx/2,-$Ly/2, $Lz/2, $Lx/2,-$Ly/2,-$Lz/2);
	printf OUT ("\n");
	
	#&OutYaplotData;
	$num ++;
	printf "$shear_rate\n";
	printf ("--------------\n");
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
		$f_normal = $fc_n + $fcol + $f_lub;
		#	$force[$k] = sqrt($f_normal)
		if ($f_normal > 3){
			if ($gap < 0){
				#$force[$k] = sqrt($f_normal*$f_normal + $fc_tan*$fc_tan);
				$force[$k] = $f_normal;
			} else {
				$force[$k] = $f_normal;
			}
		} else {
			$force[$k] = 0;
		}
		
		#$force[$k] = $fc_n + $fc_tan;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		$Gap[$k] = $gap;
		$checked[$k] = 1;

		printf OUTG "$gap ";
	}
	printf OUTG "\n";
}

sub connection {
	($i) = @_;
	$accept = 0;
	if ($i == $istart1) {
		$accept = 1;
		$step = ($radius[$istart0]+$radius[$istart1]);
		$rootpos_x = $posx[$istart0]+$step*$nrvec_x[$k_start];
		$rootpos_y = $posy[$istart0]+$step*$nrvec_y[$k_start];
		$rootpos_z = $posz[$istart0]+$step*$nrvec_z[$k_start];
		$rootforce = $force[$k_start];
	} else {
		for ($l = $#node; $l>=0; $l--){
			if (@{node[$l]->[3]} == $i){
				$rootpos_x = @{node[$l]->[0]};
				$rootpos_y = @{node[$l]->[1]};
				$rootpos_z = @{node[$l]->[2]};
				$rootforce = @{fcsegment[$l]->[6]};
				break;
				
			}
		}
	}
	while(@{int_in_p->[$i]}){
		$k = pop(@{int_in_p->[$i]}); # interaction bond of particle i
		$i0 = $int0[$k];
		$i1 = $int1[$k];
		if ($i0 != $i){
			$i_next = $i0;
			$step = -($radius[$i0]+$radius[$i1]);
		} else {
			$i_next = $i1;
			$step = ($radius[$i0]+$radius[$i1]);
		}
		if ($accept ==1 || $force[$k] < $rootforce ){
			#			printf (" %3.3f --> %3.3f\n", $force[$k], $rootforce);
			#		if ($fc_max < $force[$k]){
						#	$fc_max = $force[$k];
						#		$k_max = $k;
						#}
			$dx = $step*$nrvec_x[$k];
			$dy = $step*$nrvec_y[$k];
			$dz = $step*$nrvec_z[$k];
			$xx = $rootpos_x + $dx;
			$yy = $rootpos_y + $dy;
			$zz = $rootpos_z + $dz;
			
			push(@ilist, $i_next);
			push(@node, ([$xx, $yy, $zz, $i_next]));
			push(@fcsegment, ([$rootpos_x, $rootpos_y,$rootpos_z, $xx, $yy, $zz, $force[$k]]));
			$checked[$k] = 0;

		}
	}
}

sub analyzeForceChain {
	$maxf = 0;
	$k_max = -1;
	$num_k = 0;
	for ($k=0; $k<$num_interaction; $k++) {
		if ($checked[$k] == 1 && $force[$k] >0){
			$num_k ++;
			if ($force[$k] > $maxf) {
				$maxf = $force[$k];
				$k_max = $k;
			}
		}
	}
	#	$maxf_previous = $maxf;

	$k_start = $k_max;
	$checked[$k_start] = 0;
	printf ("$k_start $maxf $num_k\n");
	$istart0 = $int0[$k_start];
	$istart1 = $int1[$k_start];
	@ilist = ($istart1, $istart0);
	@node = ([$posx[$istart0], $posy[$istart0], $posz[$istart0], $istart0]);

	$step =$radius[$istart0]+$radius[$istart1];
	$dx = $step*$nrvec_x[$k_start];
	$dy = $step*$nrvec_y[$k_start];
	$dz = $step*$nrvec_z[$k_start];
	$xx = $posx[$istart0]+$dx;
	$yy = $posy[$istart0]+$dy;
	$zz = $posz[$istart0]+$dz;
	@fcsegment = ([$posx[$istart0], $posy[$istart0], $posz[$istart0], $xx, $yy, $zz, $force[$k_start]]);
	$it=0;
	$fc_max = 0;
	$percolation = 0;
	while ($it < @ilist ) {
		&connection(@ilist[$it]);
		$it++;
	}
	printf ("%d\n", $it);
	
	if ($#fcsegment>=5){
		printf OUT "y 4\n";
		
		$col = $#fcsegment+8;
		if ($col >= 58){
			printf "size = $col\n";
			$col = 58;
		}
		
		printf OUT ("@ %d\n", $col);
		$maxf_cluster = 0;
		
		#		printf OUT ("r 0.5\n");
		#		printf OUT ("c  %3.6f %3.6f %3.6f\n",
		#		0.5*(@{fcsegment[0]->[0]} + @{fcsegment[0]->[3]}),
		#0.5*(@{fcsegment[0]->[1]} + @{fcsegment[0]->[4]}),
		#		0.5*(@{fcsegment[0]->[2]} + @{fcsegment[0]->[5]}));
		for ($l=0; $l< $#fcsegment; $l++){
			$forcevalue = 0.0002*10*@{fcsegment[$l]->[6]};
#			if ($forcevalue > 0.3){
#				$forcevalue = 0.3;
#			}
			printf OUT ("r  %3.6f\n", $forcevalue);
			printf OUT ("s  %3.6f %3.6f %3.6f  %3.6f %3.6f %3.6f\n",
			@{fcsegment[$l]->[0]}, @{fcsegment[$l]->[1]}, @{fcsegment[$l]->[2]},
			@{fcsegment[$l]->[3]}, @{fcsegment[$l]->[4]}, @{fcsegment[$l]->[5]});
		
			if (@{fcsegment[$l]->[6]} > $maxf_cluster){
				$maxf_cluster = @{fcsegment[$l]->[6]};
			}
			
		}

	
	}
}
