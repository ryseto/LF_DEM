#!/usr/bin/perl
use Math::Trig;
use IO::Handle;

$input = $ARGV[0];


$j = index($input, '.dat', $i-1);
$name = substr($input, $i, $j-$i);
$output = "${name}_inv.dat";
printf "$output\n";
open (IN_data, "< $input");
open (OUT_data, "> $output");
$line = <IN_data>;
printf OUT_data "$line";
$line = <IN_data>;
($buf, $np1, $np2, $vf, $lx, $ly, $lz, $vf1, $vf2) = split(/\s+/, $line);
printf OUT_data "$line";
while (1){
    $line = <IN_data>;
    last unless defined $line;
    ($x, $y, $z, $a) = split(/\s+/, $line);
    $xnew = $lx-$x;
    printf OUT_data "$xnew $y $z $a\n";
}









   








