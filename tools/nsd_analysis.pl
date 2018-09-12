#!/usr/bin/perl
# Usage:
# $ nsd_analysis.pl par_[...].dat [output_interval] [xz_shift]
#
#

use Math::Trig;
use IO::Handle;
use Getopt::Long;

my $particle_data = $ARGV[0];

my $output_interval = 1;
my $xz_shift = 0;
my $axis = 0;
my $reversibility_test = 0;
my $monodisperse = 0;
my $twodim = 0;
my $pi = atan2(1, 1) * 4;
my $shear_strain_min = 1;
my $newdataformat = 1;
GetOptions(
'forcefactor=f' => \$force_factor,
'interval=i' => \$output_interval,
'shift=f' => \$xz_shift,
'axis' => \$axis,
'reversibility' => \$reversibility_test,
'monodisperse' => \$monodisperse);

printf "force_factor = $force_factor¥n";
printf "output_interval = $output_interval¥n";
printf "xz_shift = $xz_shift¥n";
printf "axis = $axis¥n";
printf "reversibility = $reversibility_test¥n";
printf "monodisperse = $monodisperse¥n";

# Create output file name
printf "$particle_data¥n";
$i = index($particle_data, "par_", 0)+4;
$iD3 = index($particle_data, "D3N");
$j = index($particle_data, ".dat", $i-1);
$name = substr($particle_data, $i, $j-$i);
my $dim = 3;
if ($iD3 eq -1){
    $dim = 2;
}
printf "$dim dimension¥n";;

$interaction_data = "int_${name}.dat";
$datafile = "data_${name}.dat";
$outputDataReconstruct = "r_$name.dat";
$outputDataContactforce = "cf_$name.dat";
$outputDataAngularDist = "ang_$name.dat";

open (OUTDR, "> ${outputDataReconstruct}");
open (OUTCF, "> ${outputDataContactforce}");
open (OUTANG, "> ${outputDataAngularDist}");

open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");
open (IN_data, "< ${datafile}");

&readHeader;

if ($Ly == 0) {
    $systemvolume = $Lx*$Lz;
} else {
    $systemvolume = $Lx*$Ly*$Lz;
}

$cnt_interval = 0;
$first = 1;
$first_int = 1;
$checkpoint = 1;
$shear_strain_previous = 0;
$shearrate_positive = 1;

$jbin = 100;
$delta_angle = $pi/$jbin;

for ($j = 0; $j < $jbin; $j ++){
    $n1bin[$j] = 0;
}

for ($j = 0; $j < $jbin; $j ++){
    $etabin[$j] = 0;
}

$n1_count = 0;
while (1) {
    $output = 1;
    &InParticles;
    last unless defined $line;
    while (1) {
        $line = <IN_data>;
        # ($time_, $css_, $sr_, $vis_all, $vis_con, $vis_dash, $vis_hydro, $vis_rep, $mf0, $mf3, $n1, $n2, $pressure) = split(/¥s+/, $line);
        ($time_, $css_) = split(/¥s+/, $line);
        #printf "$shear_strain == $css_¥n" ;
        if ($shear_strain == $css_) {
            if ($newdataformat) {
                ($time_, $css_, $sr_, $vis_all, $vis_con, $vis_dash, $vis_hydro, $vis_rep, $pressure,  $pressure_con,
                $mf0, $mf0c, $mf0d, $mf0h, $mf0r,
                $mf3, $mf3c, $mf3d, $mf3h, $mf3r,
                $ene, $mingap, $maxtandisp, $maxrollangle, $con_num, $fric_con_num, $num_int, $max_velo, $ max_ang_velo,
                $dt, $kn, $kt) = split(/¥s+/, $line);
            } else {
                ($time_, $css_, $sr_, $vis_all, $vis_con, $vis_dash, $vis_hydro, $vis_rep, $mf0, $mf3,
                $n1, $n2, $$pressure, $ene, $mingap, $maxtandisp, $maxrollangle, $con_num, $fric_con_num, $num_int, $max_velo, $ max_ang_velo,
                $dt, $kn, $kt) = split(/¥s+/, $line);
            }
            last;
        }
    }
    $viscosity = $vis_all;
    $normalstressdiff1 = -2*$mf3;
    printf "eta = $viscosity n1 = $normalstressdiff1 $mf3¥n";
    &InInteractions;
    if ($reversibility_test) {
        if ($first || $checkpoint == 1) {
            &keepInitialConfig;
        }
    }
    if ($output == 1) {
        &calcReducedN1;
    }
    $cnt_interval ++;
}

if (1) {
    $nangle = $pi + atan2($ny, $nx);
    $iangle = int($nangle/$delta_angle);
    $jbin = 100;
    for ($j = 0; $j < $jbin; $j ++){
        $angle = $delta_angle*$j + $delta_angle/2;
        $angle_diff = $angle - 3*$pi/4;
        if ($angle_diff < -$pi/2) {
            $angle_diff += $pi;
        }
        $n1 = $jbin*$n1bin[$j]/($n1_count);
        $eta = 6*$pi*$jbin*$etabin[$j]/($n1_count);
        $force = $forcebin[$j]/($n1_count);
        printf OUTANG "$angle_diff $n1 $eta $force¥n";
    }
}

close (IN_particle);
close (IN_interaction);
close (IN_data);
close (OUT);
close (OUTDR);
close (OUTCF);
close (OUTANG);

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
    $line = <IN_particle>; ($buf, $buf, $VF) = split(/\s+/, $line);
    $line = <IN_particle>; ($buf, $buf, $Lx) = split(/\s+/, $line);
    $line = <IN_particle>; ($buf, $buf, $Ly) = split(/\s+/, $line);
    $line = <IN_particle>; ($buf, $buf, $Lz) = split(/\s+/, $line);
    $line = <IN_particle>; ($buf, $buf, $flwtyp) = split(/\s+/, $line);
    $line = <IN_particle>; ($buf, $buf, $dataunit) = split(/\s+/, $line);
    
    if ($Ly == 0) {
        $number_of_header = 8;
    } else {
        $number_of_header = 7;
    }
    
    if ($newdataformat) {
        $number_of_header_data = 43;
    } else {
        $number_of_header_data = 37;
    }

    for ($i = 0; $i<$number_of_header; $i++) {
        $line = <IN_particle>;
        printf "$line";
    }
    if ($Ly == 0) {
        $number_of_header_int = 20;
    } else {
        $number_of_header_int = 20;
    }
    for ($i = 0; $i<$number_of_header_int; $i++) {
        $line = <IN_interaction>;
    }
    
    for ($i = 0; $i<$number_of_header_data; $i++) {
        $line = <IN_interaction>;
    }
    for ($i = 0; $i<$number_of_header_data; $i++) {
        $line = <IN_data>;
    }
    
    #    $xo = $Lx/2;
    #    $yo = $Ly/2;
    #    $zo = $Lz/2;
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
        ($buf, $val) = split(" : ", $line);
        
        $val =~ s/(\n|\r)//g;
        ($buf1) = split(/\s+/, $buf);
        if ($buf1 ne '#') {
            last;
        } else {
            $ssHeader[$j++] = $val;
        }
        last unless defined $line;
    }
    $shear_strain = $ssHeader[0];
    $shear_strain =~ s/(\n|\r)//g;
    printf "shear_strain = $shear_strainn";
    $shear_disp = $ssHeader[1];
    $shear_rate = $ssHeader[2];
    #    $target_stress = $ssHeader[3];
    #    $viscosity = $ssHeader[3];
    #    $normalstressdiff1 = $ssHeader[4];
    
    for ($i = 0; $i < $np; $i ++){
        if ($i > 0) {
            $line = <IN_particle>;
        }
        if ($output == 1) {
            if ($dim eq 3) {
                ($ip, $a, $x, $y, $z, $vx, $vz, $vy, $ox, $oz, $oy) = split(/\s+/, $line);
            } else {
                ($ip, $a, $x, $z, $vx, $vz, $vy, $ox, $oz, $oy, $angle) = split(/\s+/, $line);
            }
            #
            #
            $ang[$i] = $angle;
            $radius[$i] = $a;
            if ($xz_shift) {
                $x += $Lx/2;
                $z += $Lz/2;
                if ($x > $Lx/2) {
                    $x -= $Lx;
                }
                if ($z > $Lz/2) {
                    $z -= $Lz;
                }
            }
            $posx[$i] = $x-$xo;
            $posy[$i] = $y-$yo;
            $posz[$i] = $z-$zo;
            $velx[$i] = $vx;
            $vely[$i] = $vy;
            $velz[$i] = $vz;
            $omegax[$i] = $ox;
            $omegay[$i] = $oy;
            $omegaz[$i] = $oz;
            if ($radius_max < $a) {
                $radius_max = $a;
            }
        }
    }
}

sub InInteractions{
    #$line = <IN_interaction>;
    #($buf, $shear_strain_i, $num_interaction) = split(/s+/, $line);
    #printf "int $buf $shear_strain_i $num_interactionn";
    while (1) {
        $line = <IN_interaction>;
        ($buf, $val) = split(" : ", $line);
        ($buf1) = split(/\s+/, $buf);
        if ($buf1 ne '#') {
            last;
        } else {
            $val =~ s/(\n|\r)//g;
            $ssHeader[$j++] = $val;
        }
        last unless defined $line;
    }
    $k = 0;
    while (true) {
        if ($k > 0) {
            $line = <IN_interaction>;
        }
        ($i, $j, $contact, $nx, $nz, $ny, #1---6
        $gap, $f_lub_norm, # 7, 8
        $f_lub_tan_x, $f_lub_tan_z, $f_lub_tan_y, # 9, 10, 11
        $fc_norm, # 12
        $fc_tan_x, $fc_tan_z, $fc_tan_y, # 13, 14, 15
        $fr_norm, $s_xF) = split(/\s+/, $line); # 16, 17
        if ($i eq '#' || $i eq NU) {
            ($buf, $val) = split(" : ", $line);
            #            $val =~ s/(\n|\r)//g;
            #            $shear_strain = $val;
            last;
        }
        
        if (! defined $i) {
            last;
        }
        
        if ($output == 1) {
            $int0[$k] = $i;
            $int1[$k] = $j;
            $contactstate[$k] = $contact;
            $F_lub[$k] = $f_lub_norm;
            $Fc_n[$k] = $fc_norm;
            $Fc_t[$k] = sqrt($fc_tan_x**2+$fc_tan_y**2+$fc_tan_z**2);
            $Fc_tan_x[$k] = $fc_tan_x;
            $Fc_tan_y[$k] = $fc_tan_y;
            $Fc_tan_z[$k] = $fc_tan_z;
            $F_lub_t[$k] = sqrt($f_lub_tan_x**2+ $f_lub_tan_z**2 + $f_lub_tan_y**2);
            $S_bf[$k] =  $s_xF;
            #$force[$k] += $Fc_t[$k] + $F_lub_t[$k];
            $nrvec_x[$k] = $nx;
            $nrvec_y[$k] = $ny;
            $nrvec_z[$k] = $nz;
            $Gap[$k] = $gap;
            $distance[$k] = $radius[$i] + $radius[$j] + $gap;
            
            $k++;
        }
    }
    $num_interaction = $k;
}

sub calcReducedN1 {
    if ($first == 0) {
        printf OUT "\n";
    } else {
        $first = 0;
    }
    if ($shear_strain > $shear_strain_min) {
        $n1_count ++;
    }
    ## visualize force chain network
    printf OUT "y 6\n";
    printf OUT "@ 7\n";
    $etaTotal = 0;
    $eta_denominator = $systemvolume;
    $cf_normal_sum = 0;
    $cf_tangent_sum = 0;
    $cf_normal_max = 0;
    $cf_tangent_max = 0;
    
    $cf_cnt = 0;
    for ($k = 0; $k < $num_interaction; $k ++) {
        $force = $F_lub[$k] + $Fc_n[$k];
        $nx = $nrvec_x[$k];
        $ny = $nrvec_y[$k];
        $nz = $nrvec_z[$k];
        ### theta = atan2($ny, $nx) <--- tan theta = y/x
        #phi = atan2($ny, $nz);
        $eta_lub = -$F_lub[$k]*$distance[$k]*($nx*$ny);
        $eta_con = - $Fc_n[$k]*$distance[$k]*($nx*$ny);
        $eta = $eta_lub + $eta_con;
        $etaTotal += $eta;
        if ($Fc_n[$k] > 0) {
            #            printf "$Fc_n[$k] $Fc_t[$k] \n";
            $cf_ratio_sum += $Fc_t[$k]/$Fc_n[$k];
            $cf_normal_sum += $Fc_n[$k];
            $cf_tangent_sum += $Fc_t[$k];
            $cf_cnt ++;
            if ($cf_normal_max  < $Fc_n[$k]) {
                $cf_normal_max  = $Fc_n[$k];
                $cf_tangent_max = $Fc_t[$k];
            }
        }
        
        if ($shear_strain > $shear_strain_min) {
            $nangle = $pi + atan2($ny, $nx);
            $iangle = int($nangle/$delta_angle);
            $etabin[$iangle] += 6*$pi*$eta/$eta_denominator;
            $projxy = sqrt($nx*$nx +$ny*$ny);
            $forcebin[$iangle] += $force*$projxy;
        }
    }
    if ($cf_cnt > 0 ) {
        $cf_normal = $cf_normal_sum/$cf_cnt;
        $cf_tangent = $cf_tangent_sum/$cf_cnt;
    } else {
        printf "CF =  nocontact\n";
    }
    if ($shear_strain > 1) {
        printf OUTCF "$shear_strain $cf_normal $cf_tangent $cf_normal_max $cf_tangent_max\n"
    }
    $visapprox = $etaTotal/$systemvolume;
    $n1approx_pos_con = 0;
    $n1approx_neg_con = 0;
    $n1approx_pos_lub = 0;
    $n1approx_neg_lub = 0;
    $n1con_fric_sum = 0;
    # printf "viscosity = $viscosityn";
    #  $force_factor = 0.015/$viscosity;
    printf OUT "y 5\n";
    #    printf OUT "@ 7\n";
    $n1_denominator = $systemvolume*$viscosity;
    for ($k = 0; $k < $num_interaction; $k ++) {
        $nx = $nrvec_x[$k];
        $ny = $nrvec_y[$k];
        $n1lub = -$F_lub[$k]*$distance[$k]*($nx*$nx - $ny*$ny);
        $n1con = -$Fc_n[$k]*$distance[$k]*($nx*$nx - $ny*$ny);
        $n1con_fric_sum += $distance[$k]*($nx*$Fc_tan_x[$k] - $ny*$Fc_tan_y[$k]);
        
        if ($n1lub > 0) {
            $n1approx_pos_lub += $n1lub;
        } else {
            $n1approx_neg_lub += $n1lub;
        }
        if ($n1con > 0) {
            $n1approx_pos_con += $n1con;
        } else {
            $n1approx_neg_con += $n1con;
        }
        $n1 = $n1lub + $n1con;
        if ($shear_strain > $shear_strain_min) {
            $nangle = $pi + atan2($ny, $nx);
            $iangle = int($nangle/$delta_angle);
            $n1bin[$iangle] += 6*$pi*$n1/$n1_denominator;
        }
    }
    $n1approxCon =  ($n1approx_pos_con + $n1approx_neg_con)/$systemvolume;
    $n1approxLub = ($n1approx_pos_lub + $n1approx_neg_lub)/$systemvolume;
    if ($shear_strain >= $shear_strain_min) {
        $visapproxOut = 6*$pi*$visapprox;
        $normalstressdiff1ratio = $normalstressdiff1/$viscosity;
        $n1approxRatio = 6*$pi*($n1approxCon + $n1approxLub)/$viscosity;
        printf OUTDR "$shear_strain $viscosity $visapproxOut $normalstressdiff1ratio $n1approxRatio\n";
    }
}


