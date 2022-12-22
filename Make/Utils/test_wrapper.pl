#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long 'HelpMessage';
use File::Basename qw(dirname basename);

GetOptions(
    'test|t=s' => \(my $test_file),
    'mpi!'   => \(my $mpi = 0),
    'help'     =>   sub { HelpMessage(2) },
 ) or HelpMessage(1);

system("rm -f debug err_0 out_0");
my $input_file = "$test_file.in";

print ("[\e[1m $test_file \e[0m ] ... ");
my $mpicmd="";

if ($mpi) {
    open my $FH, $input_file or die "Could not open input file [$input_file]: $!\n";
    my ($NT, $NX, $NY, $NZ) = (1,1,1,1);
    my $count = 0;
    while( <$FH> )  {   
        if (/NP_T\s*=\s*(\d+)/) { $NT = $1; $count++; }    
        if (/NP_X\s*=\s*(\d+)/) { $NX = $1; $count++; }    
        if (/NP_Y\s*=\s*(\d+)/) { $NY = $1; $count++; }    
        if (/NP_Z\s*=\s*(\d+)/) { $NZ = $1; $count++; }    
        last if $count == 4;
    }
    close $FH;

    my $NP=$NT*$NX*$NY*$NZ;
    # print "NP=$NP ($NT,$NX,$NY,$NZ)\n";
    $ENV{"OMPI_MCA_btl_base_warn_component_unused"}=0;
    $mpicmd="mpirun -np $NP";
}

my $testdir = dirname($test_file);
my $testbn = basename($test_file);
my $ret = system("rm -f $testdir/.test_failed_$testbn && $mpicmd $test_file -i $input_file 2>&1 >/dev/null");
if ($ret) {
    print ("\e[1;31mFAILED\e[0m\n");
    system("touch $testdir/.test_failed_$testbn && tail -n 80 out_0");
    exit(0);
}
print ("PASSED\n");

exit(0);

=head1 NAME

test_wrapper - run test for HiRep

=head1 SYNOPSIS

  Option              Default     Description
  ----------------------------------------------------------------------------
  --help,-h                       Print this help

  --test,-t                       Name of test to run
  --[no]-mpi          [false]     Run test with mpirun. Number of processes is 
                                  read from input file

=head1 VERSION

1.0

=cut
