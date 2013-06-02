#!/usr/bin/perl
use Math::Trig;
$initcongif_data = $ARGV[0];
# Create output file name
$j = index($initcongif_data, '.dat');
printf "$j\n";
$name = substr($initcongif_data , 0, $j);
printf "${name}\n";
$output_name = "${name}_s$scale_factor.dat";
printf "$output_name\n";
open (OUT, "> ${output_name}");
open (IN, "< $initcongif_data");
$line = <IN>;
printf OUT "$line";
$line = <IN>;

($buf1, $n1, $n2, $vol, $lx, $ly, $lz) = split(/\s+/, $line);
$lx2 = $lx*$scale_factor;
$ly2 = $ly*$scale_factor;
$lz2 = $lz*$scale_factor;
printf OUT "$buf1 0 $n1 $vol $lx2 $ly2 $lz2\n";
$n = $n1+$n2;
for ($i = 0; $i < $n; $i ++) {
    $line = <IN>;
    ($x, $y, $z, $a) = split(/\s+/, $line);
    $sx = $x*$scale_factor;
    $sy = $y*$scale_factor;
    $sz = $z*$scale_factor;
    $sa = $a*$scale_factor;
    printf OUT "$sx $sy $sz $sa\n";
}
exit;
    
