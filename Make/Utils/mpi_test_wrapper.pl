#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my $input_file;

GetOptions(
    'i=s' => \$input_file,
 ) or die "Usage: $0 -i input_filename\n";

my $NT = `awk '/NP_T/ {printf "%d",\$3}' $input_file`;
if ( $NT<1 ) { $NT = 1 }

my $NX = `awk '/NP_X/ {printf "%d",\$3}' $input_file`;
if ( $NX<1 ) { $NX = 1 }
#
my $NY = `awk '/NP_Y/ {printf "%d",\$3}' $input_file`;
if ( $NY<1 ) { $NY = 1 }
#
my $NZ = `awk '/NP_Z/ {printf "%d",\$3}' $input_file`;
if ( $NZ<1 ) { $NZ = 1 }
#
my $NP=$NT*$NX*$NY*$NZ;
#print "NP=$NP ($NT,$NX,$NY,$NZ)\n";
exec("mpirun", "-np", $NP, @ARGV,"-i",$input_file);

