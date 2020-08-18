#!/usr/bin/perl -w
use strict;

my $line;
my ($tc,$DeltaS,$MVM,$rl,$rh,$plaq,$acc);
my $oldtc=0;
my @output;

undef($tc);
$acc=-1;
while($line=<>) {
        if(eof) {delerr(\@output); @output=(); $oldtc=0;}
        if ($line=~/.*?Trajectory #(\d+)\.\.\./) {$tc=$1;}
        if ($line=~/.*DeltaS = (.*?)\]/) {$DeltaS=$1; }
        if ($line=~/.*Configuration accepted/) {if ($acc!=-1) {print "Error acc\n";}  $acc=1; }
        if ($line=~/.*Configuration rejected/) {if ($acc!=-1) {print "Error acc\n";}  $acc=0; }
        if ($line=~/.*Range = \[(.*?),(.*?)\]/) {$rl=($1+1.e-4)/0.9; $rh=$2/1.1;} #normalize range
        if ($line=~/.*?Trajectory #(\d+): .*? MVM = (\d+)/) {
                $MVM=$2; 
                if($tc!=$1) {
                        die ("Error: trajectory numbers do not match\n");
                }
                if($tc<=$oldtc) { push @output, "Error\n"; }
                $oldtc=$tc;
        }
        if ($line=~/.*Plaquette: (.*)/) {
                my $l;
                $plaq=$1; 
                $l=sprintf("%d %1.8e %1.8e %1.8e %1.8e %1.8e %d\n",$tc,$plaq,$MVM,$rl,$rh,$DeltaS,$acc);
                push @output,$l; 
                undef($tc); 
                $acc=-1;
        }
}

sub delerr {
my $input=shift;
my @out;

while($_=shift(@$input)) {
        if ($_=~/Error/) {
                my ($next,$id,$prev,$id2);
                $next=shift(@$input);
                defined($next) or die("No line after Error!\n");
                $next=~/(\d+)\s+.*/ or die("Malformed line!\n[$next]");
                $id=$1;
                do{
                        $prev=pop(@out);
                        $prev=~/(\d+)\s+.*/ or die("Malformed line [2]!\n[$prev]");
                        $id2=$1;
                } while ($id2!=$id);
                push @out,$next;
        } else {
                push @out,$_;
        }
}

print @out;
}

