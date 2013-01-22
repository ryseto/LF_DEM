#! /usr/bin/perl
use Math::Trig;
$particle_data = $ARGV[0];
$interaction_data = $ARGV[1];
$force_factor = $ARGV[2];
open (IN_particle, "< ${particle_data}");
open (IN_interaction, "< ${interaction_data}");

# Read the header
#
$line = <IN_particle>; ($buf, $VF) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $Lx) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $Ly) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $Lz) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $np_a) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $np_b) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $radius_a) = split(/\s+/, $line);
$line = <IN_particle>; ($buf, $radius_b) = split(/\s+/, $line);
printf "$VF, $Lx, $Ly, $Lz, $np_a, $np_b, $radius_a, $radius_b\n";

# Create output file name
$i = index($particle_data, 'par_', 0)+4;
$j = index($particle_data, '.dat', $i-1);
$name = substr($particle_data, $i, $j-$i);
$output = "$name.yap";
open (OUT, "> ${output}");

$np = $np_a + $np_b;
while ( 1 ){
    $line = <IN_particle>; 
    ($buf, $shear_rate) = split(/\s+/, $line);
    printf "$shear_rate\n";
    for ($i = 0; $i < $np; $i ++){
        $line = <IN_particle> ;
        ($i, $x, $y, $z, $vx, $vy, $vz, $ox, $oy, $oz ) = split(/\s+/, $line);
        $posx[$i] = $x; 
        $posy[$i] = $y; 
        $posz[$i] = $z;
    }
    if ($buf != "#"){
        exit;
    }

    $line = <IN_interaction>; 
    ($buf, $shear_rate, $num_interaction) = split(/\s+/, $line);
    if ($buf != "#"){
        exit;
    }
    for ($k = 0; $k < $num_interaction; $k ++){
        $line = <IN_interaction> ;
        ($i, $j, $r, $f_lub, $fc_n, $fc_t, $fric_st) = split(/\s+/, $line);
        $int0[$k] = $i;
        $int1[$k] = $j;
        $F_lub[$k] = $f_lub;
        $Fc_n[$k] = $fc_n;
        $Ft_t[$k] = $fc_t;
    }
    printf OUT "y 1\n";
    printf OUT "@ 2\n";
    for ($i = 0; $i < $np; $i ++){
        printf OUT "c $posx[$i] $posy[$i] $posz[$i] \n";
    }
    printf OUT "y 2\n";
    printf OUT "@ 6\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        if ($Fc_n[$k] > 0){
            $string_with = $force_factor*($Fc_n[$k]);
            printf OUT "r ${string_with}\n";
            OutString($int0[$k],  $int1[$k]);
        }
    }
    printf OUT "y 3\n";
    printf OUT "@ 3\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        $force = $F_lub[$k] + $Fc_n[$k];
        if ($force > 0){
            $string_with = ${force_factor}*$force;
            printf OUT "r ${string_with}\n";
            OutString($int0[$k],  $int1[$k]);
        } 
    }    
    printf OUT "@ 4\n";
    for ($k = 0; $k < $num_interaction; $k ++){
        $force = $F_lub[$k] + $Fc_n[$k];
        if ($force < 0){
            $string_with = -${force_factor}*$force;
            printf OUT "r ${string_with}\n";
            OutString($int0[$k],  $int1[$k]);
        } 
    }
    printf OUT "\n";
}
close (OUT);
exit;    

sub OutString {
    ($i, $j) = @_;
    $xi = $posx[$i];
    $yi = $posy[$i]; 
    $zi = $posz[$i];
    $xj = $posx[$j];
    $yj = $posy[$j]; 
    $zj = $posz[$j];
    if (abs($xi-$xj) < 3
        &&  abs($yi-$yj) < 3
        &&  abs($zi-$zj) < 3){
        printf OUT "s $xi $yi $zi $xj $yj $zj\n";
    }
}

