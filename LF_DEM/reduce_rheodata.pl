#!/usr/bin/perl

# Usage:
# $ generateYaplotFile.pl par_[...].dat [force_factor] [y_section]
#
# force_factor: To change the widths of strings to exhibit forces.
# y_section: Visualize the trimed range of the y coordinates.

use Math::Trig;
use IO::Handle;
#use strict;
#use warnings;

$reducedtimes = 10;
if ($#ARGV >= 1){
	$reducedtimes = $ARGV[1];
}

$rheo_data = $ARGV[0];
$newname = "${rheo_data}_orig";
printf "$newname\n";
unless (rename $rheo_data, $newname) {
	printf "failed (%s => %s)\n", $rheo_data, $newname;
	printf "error(%d:%s)\n", $!, $!;
	exit 1;
}
open (IN_rheo, "< $newname");
open (OUT_rheo, "> $rheo_data");
my $i = 0;
while (1) {
	$line = <IN_rheo>;
	($buf, $val) = split(/\s+/, $line);
	last unless defined $line;
	if ($buf == "#") {
		printf OUT_rheo "$line";
	} else {
		if ($i % $reducedtimes == 0) {
			printf OUT_rheo "$line";
		}
		$i += 1;
	}
}



