#!/usr/bin/perl -w
use strict;

print "#ifndef _APPROX_DATA\n#define _APPROX_DATA\n\n\n";
print "static const double approx_data[]={";

my $n=3;
my $tab="\n   ";
my $c=0;
while ($_=<STDIN>) {
	my @app;
	chomp $_;
	@app=split /,/;
	for (my $i=0;$i<@app;$i++){
		if($c%$n == 0) {
			print $tab;
		}
		$c++;
		print "$app[$i],";
	}
}
if($c%$n == 0) {
	print $tab;
}
print "0.\n";

print "};\n\n";

print "#endif\n";

