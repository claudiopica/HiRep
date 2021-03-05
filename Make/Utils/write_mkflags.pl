#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long 'HelpMessage';
use File::Copy qw(move);

GetOptions(
  'file|f=s' => \( my $file = 'MkFlags'),
  'ng|n=i' => \( my $NG = 2),
  'repr|r=s'   => \(my $repr = 'FUND'),
  'gauge=s'   => \(my $gauge = 'SUN'),
  'mpi!'   => \(my $mpi = 1),
  't=s'   => \(my $TBC = 'P'),
  'x=s'   => \(my $XBC = 'P'),
  'y=s'   => \(my $YBC = 'P'),
  'z=s'   => \(my $ZBC = 'P'),
  'twisted!'   => \(my $xyz_twist = 0),
  'sf!'   => \(my $sfbc = 0),
  'sfhalf!'   => \(my $sfhalfbc = 0),
  'sfrotated!'   => \(my $sfrotatedbc = 0),
  'smearing!'   => \(my $smearing = 0),
  'clover|c!'   => \(my $clover = 0),
  'expclover|e!'   => \(my $expclover = 0),
  'eo!'   => \(my $eoprec = 1),
  'quat|q!'   => \(my $quat = 0),
  'ndebug!'   => \(my $ndebug = 1),
  'dfloat!'   => \(my $dfloat = 0),
  'checkspinor!'   => \(my $scheck = 1),
  'mpitiming!'   => \(my $mpit = 0),
  'ioflush!'   => \(my $iof = 1),
  'unrollrepr!'   => \(my $unrollr = 0),
  'timing!'   => \(my $timing = 0),
  'bartiming!'   => \(my $btiming = 0),
  'memory!'   => \(my $mem = 0),
  'force!'   => \(my $force = 0),
  'cc=s'   => \(my $cc = "gcc"),
  'mpicc=s'   => \(my $mpicc = "mpicc"),
  'cflags=s'   => \(my $cflags = "-Wall -std=c99 -O3"),
  'ldflags=s'   => \(my $ldflags = ""),
  'include=s'   => \(my $include = ""),
  'ccache!'   => \(my $ccache = 0),
  'help'     =>   sub { HelpMessage(2) },
) or HelpMessage(1);

# validate parameters
validate_ng();
validate_repr();
validate_gauge();
validate_t();
validate_x();
validate_y();
validate_z();

sub validate_ng {
    if($NG<2) {
        print "Error: Number of colors (--ng|n) must be 2 or bigger\n";
        HelpMessage(1);
    }
}

sub validate_repr {
    if($repr eq "FUND") {
        $repr = "REPR_FUNDAMENTAL"
    } elsif ($repr eq "2S") { 
        $repr = "REPR_SYMMETRIC"        
    } elsif ($repr eq "2A") { 
        $repr = "REPR_ANTISYMMETRIC"
    } elsif ($repr eq "ADJ") { 
        $repr = "REPR_ADJOINT"
    } else {
        print "Error: The fermion representation must be one of the following: FUND, 2S, 2A, ADJ\n";
        HelpMessage(1);
    }
}

sub validate_gauge {
    if($gauge eq "SUN") {
        $gauge = "GAUGE_SUN"
    } elsif ($gauge eq "SON") { 
        $gauge = "GAUGE_SON"         
    } else {
        print "Error: The gauge group must be one of the following: SUN, SON\n";
        HelpMessage(1);
    }
}

sub validate_t {
    if($TBC eq "P") {
        $TBC = "BC_T_PERIODIC"
    } elsif ($TBC eq "A") { 
        $TBC = "BC_T_ANTIPERIODIC"
    } elsif ($TBC eq "T") { 
        $TBC = "BC_T_THETA"
    } elsif ($TBC eq "O") { 
        $TBC = "BC_T_OPEN"
    } else {
        print "Error: The T boundary condition representation must be one of the following: P, A, T, O\n";
        HelpMessage(1);
    }
}

sub validate_x {
    if($XBC eq "P") {
        $XBC = "BC_X_PERIODIC"
    } elsif ($XBC eq "A") { 
        $XBC = "BC_X_ANTIPERIODIC"
    } elsif ($XBC eq "T") { 
        $XBC = "BC_X_THETA"
    } else {
        print "Error: The X boundary condition representation must be one of the following: P, A, T\n";
        HelpMessage(1);
    }
}

sub validate_y {
    if($YBC eq "P") {
        $YBC = "BC_Y_PERIODIC"
    } elsif ($YBC eq "A") { 
        $YBC = "BC_Y_ANTIPERIODIC"
    } elsif ($YBC eq "T") { 
        $YBC = "BC_Y_THETA"
    } else {
        print "Error: The Y boundary condition representation must be one of the following: P, A, T\n";
        HelpMessage(1);
    }
}

sub validate_z {
    if($ZBC eq "P") {
        $ZBC = "BC_Z_PERIODIC"
    } elsif ($ZBC eq "A") { 
        $ZBC = "BC_Z_ANTIPERIODIC"
    } elsif ($ZBC eq "T") { 
        $ZBC = "BC_Z_THETA"
    } else {
        print "Error: The Z boundary condition representation must be one of the following: P, A, T\n";
        HelpMessage(1);
    }
}

sub backup_oldflags {
    my $suff = '.bak';
    my $bkpfile="$file$suff";
    move $file, $bkpfile;
}

backup_oldflags();
open(my $fh, '>', $file) or die "Could not open file '$file' !";
# write ng
print $fh "NG = $NG\n";
# write repr
print $fh "REPR = $repr\n"; 
# write gauge
print $fh "GAUGE_GROUP = $gauge\n";
# write T boundary condition
print $fh "MACRO += -D$TBC\n";
# write X boundary condition
print $fh "MACRO += -D$XBC\n";
# write Y boundary condition
print $fh "MACRO += -D$YBC\n";
# write Z boundary condition
print $fh "MACRO += -D$ZBC\n";
# write twisted boundary condition
$xyz_twist && print $fh "MACRO += -DBC_XYZ_TWISTED\n";
# write sf boundary condition
$sfbc && print $fh "MACRO += -DBASIC_SF\n";
# write sf half field boundary condition
$sfhalfbc && print $fh "MACRO += -DHALFBG_SF\n";
# write sf rotated boundary condition
$sfrotatedbc && print $fh "MACRO += -DROTATED_SF\n";
# write smearing
$smearing && print $fh "MACRO += -DWITH_SMEARING\n";
# write clover
$clover && print $fh "MACRO += -DWITH_CLOVER\n";
# write expclover
$expclover && print $fh "MACRO += -DWITH_EXPCLOVER\n";
# write eo preconditioning
$eoprec && print $fh "MACRO += -DWITH_EO\n";
# write quaternions 
$quat && print $fh "MACRO += -DWITH_QUATERNIONS\n";
# write ndebug
$ndebug && print $fh "MACRO += -DNDEBUG\n";
# write dphi float
$dfloat && print $fh "MACRO += -DDPHI_FLOAT\n";
# write check spinor
$scheck && print $fh "MACRO += -DCHECK_SPINOR_MATCHING\n";
# write mpitiming
$mpit && print $fh "MACRO += -DMPI_TIMING\n";
# write io flush
$iof && print $fh "MACRO += -DIO_FLUSH\n";
# write unroll representation
$unrollr && print $fh "MACRO += -DUNROLL_GROUP_REPRESENT\n";
# write timing
$timing && print $fh "MACRO += -DTIMING\n";
# write timing
$btiming && print $fh "MACRO += -DTIMING_WITH_BARRIERS\n";
# write memory
$mem && print $fh "MACRO += -DAMALLOC_MEASURE\n";
# write force
$force && print $fh "MACRO += -DMEASURE_FORCE\n";
# write mpi
$mpi && print $fh "MACRO += -DWITH_MPI\n";
# write compiler options
if ($ccache!=0) { $cc="ccache ".$cc; $mpicc="ccache ".$mpicc; }
print $fh "CC = $cc\n";
print $fh "MPICC = $mpicc\n";
print $fh "CFLAGS = $cflags\n";
print $fh "INCLUDE = $include\n";
print $fh "LDFLAGS = $ldflags\n";

close $fh;
#print "All done\n";

=head1 NAME

write_mkflags - write flags file for compilation of HiRep

=head1 SYNOPSIS

  Option              Default     Description
  ----------------------------------------------------------------------------
  --help,-h                       Print this help

  --file,-f           [MkFlags]   File name
  --ng,-n             [2]         Number of colors
  --repr,-r           [FUND]      Fermion representation (FUND, 2S, 2A, ADJ)
  --gauge,-g          [SUN]       Gauge group (SUN, SON)
  -t                  [P]         T boundary conditions (P, A, T, O)
  -x                  [P]         X boundary conditions (P, A, T)
  -y                  [P]         Y boundary conditions (P, A, T)
  -z                  [P]         Z boundary conditions (P, A, T)

  --[no-]mpi          [true]      Use MPI
  --cc                [gcc]       Compiler
  --mpicc             [mpicc]     MPI Compiler
  --cflags            [-O3]       Compilation options
  --include           []          Extra include headers
  --ldflags           []          Linking options
  --[no-]ndebug       [true]      Set ndebug flag
  --[no-]ccache       [false]     Use ccache

  --[no-]eo           [true]      Even-Odd preconditioning
  
  --[no-]twist        [false]     XYZ twisted boundary conditions
  --[no-]sf           [false]     Schrodinger functional b.c.
  --[no-]sfhalf       [false]     Schrodinger functional b.c., half field
  --[no-]sfrotate     [false]     Rotated Schrodinger functional b.c.
  --[no-]smearing     [false]     Smearing action
  --[no-]clover,-c    [false]     Clover improved action
  --[no-]expclover,-e [false]     ExpClover improved action

  --[no-]quat,-q      [false]     Use quaternion representation (only for SU2)
  --[no-]dfloat       [false]     Use single precision acceleration
  --[no-]unrollrepr   [false]     Unroll group representation functions

  --[no-]checkspinor  [true]      Check spinor field type
  --[no-]mpitiming    [false]     Enable timing of MPI calls
  --[no-]ioflush      [true]      Flush IO after each operations on logs
  --[no-]timing       [false]     Enable timing
  --[no-]bartiming    [false]     Enable MPI barriers in timing
  --[no-]memory       [false]     Print memory usage
  --[no-]force        [false]     Print statics for forces in molecular dynamics

=head1 VERSION

1.0

=cut