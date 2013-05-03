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
$counter_fc = 0;
$istart0 = -1;
$istart1 = -1;
$k_max = -1;
while (1){
	&InParticles;
	&InInteractions;
	last unless defined $line;
	&analyzeForceChain;
	
	#&OutYaplotData;
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
		#$force[$k] = $fc_n + $fc_tan;
		$nrvec_x[$k] = $nx;
		$nrvec_y[$k] = $ny;
		$nrvec_z[$k] = $nz;
		$Gap[$k] = $gap;
		printf OUTG "$gap ";
	}
	printf OUTG "\n";
}



sub connection {
	($i) = @_; # 分岐の粒子
	for ($l=$#chain; $l>=0; $l--){
		if (@{chain[$l]->[3]} == $i){
			$rootpos_x = @{chain[$l]->[0]};
			$rootpos_y = @{chain[$l]->[1]};
			$rootpos_z = @{chain[$l]->[2]};
			break;
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
		$exist = 0;
		if ($#ilist > 5){
			if ($i_next == $istart0){
				$percolation = 1;
			}
		}
		
		for ($l=0; $l<$#ilist; $l++){
			if ($i_next == @ilist[$l]){
				$exist = 1;
				break;
			}
		}
		if ($exist == 0){
			if ($fc_max < $force[$k]){
				$fc_max = $force[$k];
				$k_max = $k;
			}
			push(@ilist, $i_next);
			$dx = $step*$nrvec_x[$k];
			$dy = $step*$nrvec_y[$k];
			$dz = $step*$nrvec_z[$k];
			$xx = $rootpos_x + $dx;
			$yy = $rootpos_y + $dy;
			$zz = $rootpos_z + $dz;
			push(@chain, ([$xx,$yy,$zz,$i_next]));
			push(@fcsegment, ([$rootpos_x, $rootpos_y,$rootpos_z, $xx, $yy, $zz, $force[$k]]));
			
		}
	}
	
}


sub analyzeForceChain {
	for ($i=0; $i<$np; $i++) {
		@{int_in_p[$i]} = ();
	}


	$maxf = 0;
	$k_max = -1;
	for ($k=0; $k<$num_interaction; $k++) {
		if ($force[$k] > $maxf) {
			$maxf = $force[$k];
			$k_max = $k;
		}
	}
	$k_start = $k_max;
	$istart0 = $int0[$k_start];
	$istart1 = $int1[$k_start];
	
	for ($k=0; $k<$num_interaction; $k++) {
		$ii0 = $int0[$k];
		$ii1 = $int1[$k];
		if ($force[$k] > 0.333*$maxf){
			push(@{int_in_p->[$ii0]}, $k);
			push(@{int_in_p->[$ii1]}, $k);
		}
	}
	printf "$maxf\n";
	
	
	if ($istart0 != -1){
		@ilist = ($istart0);
		@chain = ([0 ,0, 0, $istart0]);
		
		$dx = $radius[$istart0]*$nrvec_x[$k_start];
		$dy = $radius[$istart0]*$nrvec_y[$k_start];
		$dz = $radius[$istart0]*$nrvec_z[$k_start];
		@fcsegment = ([0,0,0, $dx,$dy,$dz,$force[$k_start]] );
#		push(@ilist, $istart1);

#		push(@chain, ([$dx, $dy, $dz , $istart1]));
		$it=0;
		$fc_max = 0;
		$percolation = 0;
		while ($it < @ilist) {
			&connection(@ilist[$it]);
			$it++;
		}

		#	while(@ilist){
		#		print pop(@ilist);
		#		print ",";
		#}
		#print "\n";
		
		#	push(@{Pos->[0]}, 2);
		#push(@{Pos->[0]}, 3);
		
		#print @[Pos->[0]];
		#	for ($i=0; $i<$np; $i++) {
		#		if (@{Pos->[$i]}) {
		#			while(@{Pos->[$i]}){
		##				print pop(@{Pos->[$i]});
		#			print ",";
		#			}
		#			print "\n";
		#		}
		#}
#		if ($#chain>2){
#			for ($l=0; $l< $#chain; $l++){
#				printf OUT ("r %3.5f\n", $radius[@{chain[$l]->[3]}]);
#				printf OUT ("c %3.5f %3.5f %3.5f\n", @{chain[$l]->[0]}, @{chain[$l]->[1]}, @{chain[$l]->[2]});
#				
#			}
#			printf OUT ("t 0 0 28 %d\n", $#chain);
#			if ($percolation == 1){
#				printf OUT ("t 0 0 30 percolation\n");
#			}
#			if ($counter_fc == 0){
#				printf OUT ("t 0 0 32 renew origin\n");
#			}
#			printf OUT ("\n");
#		}
		if ($#chain>2){
			printf OUT "@ 4\n";
			for ($l=0; $l< $#fcsegment; $l++){
				printf OUT ("r  %3.5f\n", 0.002*@{fcsegment[$l]->[6]});
				printf OUT ("s  %3.5f %3.5f %3.5f  %3.5f %3.5f %3.5f\n",
				@{fcsegment[$l]->[0]}, @{fcsegment[$l]->[1]}, @{fcsegment[$l]->[2]},
				@{fcsegment[$l]->[3]}, @{fcsegment[$l]->[4]}, @{fcsegment[$l]->[5]});
			}
			printf OUT ("\n");
			#		if ($counter_fc++ > 1000){
			#			$istart0 = $int0[$k_max];
#			$istart1 = $int1[$k_max];
#			$k_start = $k_max;
#			$counter_fc = 0;
#			printf "k = $k_max\n";
#		}
		}
		#
		#		printf OUT ("@ 2\n");
		#		for ($l=0; $l<@ilist; $l++){
		#			$i = @ilist[$l];
		#			if ($i != $istart0 && $i != $istart1) {
		#				printf OUT ("r %3.5f\n", $radius[$i]);
		#				printf OUT ("c %3.5f %3.5f %3.5f\n", $posx[$i]-$xo, $posy[$i]-$yo, $posz[$i]-$zo);
		#			}
		#		}
		#
		#		printf OUT ("\n");
		#		printf OUT ("@ 5\n");
		#		printf OUT ("r %3.5f\n", $radius[$istart0]);
		#		printf OUT ("c %3.5f %3.5f %3.5f\n", $posx[$istart0]-$xo, $posy[$istart0]-$yo, $posz[$istart0]-$zo);
		#		printf OUT ("r %3.5f\n", $radius[$istart1]);
		#		printf OUT ("c %3.5f %3.5f %3.5f\n", $posx[$istart1]-$xo, $posy[$istart1]-$yo, $posz[$istart1]-$zo);
	}
}
