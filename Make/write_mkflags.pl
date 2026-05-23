#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long 'HelpMessage';
use File::Copy qw(move);

GetOptions(
  'file|f=s'     => \(my $file = 'MkFlags.ini'),
  'ng|n=i'       => \(my $NG = 2),
  'repr|r=s'     => \(my $repr = 'FUND'),
  'gauge=s'      => \(my $gauge = 'SUN'),
  'mpi!'         => \(my $mpi = 1),
  'gpu!'         => \(my $gpu = 0),
  'omp!'         => \(my $omp = 0),
  't=s'          => \(my $TBC = 'P'),
  'x=s'          => \(my $XBC = 'P'),
  'y=s'          => \(my $YBC = 'P'),
  'z=s'          => \(my $ZBC = 'P'),
  'twisted!'     => \(my $xyz_twist = 0),
  'sfhalf!'      => \(my $sfhalfbc = 0),
  'smearing!'    => \(my $smearing = 0),
  'clover|c!'    => \(my $clover = 0),
  'expclover|e!' => \(my $expclover = 0),
  'eo!'          => \(my $eoprec = 1),
  'newgeo!'      => \(my $newgeo = 0),
  'fixedstr!'    => \(my $fixedstr = 0),
  'quat|q!'      => \(my $quat = 0),
  'ndebug!'      => \(my $ndebug = 1),
  'dfloat!'      => \(my $dfloat = 0),
  'checkspinor!' => \(my $scheck = 1),
  'mpitiming!'   => \(my $mpit = 0),
  'ioflush!'     => \(my $iof = 1),
  'logallpids!'  => \(my $logallpids = 0),
  'unrollrepr!'  => \(my $unrollr = 0),
  'timing!'      => \(my $timing = 0),
  'bartiming!'   => \(my $btiming = 0),
  'memory!'      => \(my $mem = 0),
  'force!'       => \(my $force = 0),
  'avx2!'        => \(my $avx2 = 0),
  'vect!'        => \(my $vect = 0),
  'fuse!'        => \(my $fuse = 0),
  'probempi!'    => \(my $probempi = 0),
  'hwloc!'       => \(my $hwloc = 0),
  'env=s'        => \(my $env = ""),
  'cc=s'         => \(my $cc = "gcc"),
  'mpicc=s'      => \(my $mpicc = "mpicc"),
  'cflags=s'     => \(my $cflags = "-Wall -O3"),
  'nvcc=s'       => \(my $nvcc = "nvcc"),
  'gpuflags=s'   => \(my $gpuflags = ""),
  'ldflags=s'    => \(my $ldflags = ""),
  'include=s'    => \(my $include = ""),
  'ccache!'      => \(my $ccache = 0),
  'color!'       => \(my $color = 1),
  'help'         => sub { HelpMessage(2) },
) or HelpMessage(1);

# validate parameters
validate_ng();
validate_quat();
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

sub validate_quat {
    if($quat && ($NG != 2)) {
        print "Error: Quaternions (--quat|q) can only be used when colors (--ng|n) is 2\n";
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
    } elsif ($TBC eq "S") { 
        $TBC = "BC_T_SF"
    } elsif ($TBC eq "SR") { 
        $TBC = "BC_T_SF_ROTATED"
    } elsif ($TBC eq "M") { 
        $TBC = "BC_T_MIXED"
    } else {
        print "Error: The T boundary condition representation must be one of the following: P, A, T, O, S, SR, M\n";
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
open(my $fh, '>', $file) or die "Can not open '$file' $!";
# write ng
print $fh "NG = $NG\n";
# write repr
print $fh "REPR = $repr\n"; 
# write gauge
print $fh "GAUGE_GROUP = $gauge\n";
# write T boundary condition
print $fh "MACRO += $TBC\n";
# write X boundary condition
print $fh "MACRO += $XBC\n";
# write Y boundary condition
print $fh "MACRO += $YBC\n";
# write Z boundary condition
print $fh "MACRO += $ZBC\n";
# write twisted boundary condition
$xyz_twist && print $fh "MACRO += GAUGE_SPATIAL_TWIST\n";
# write sf half field boundary condition
$sfhalfbc && print $fh "MACRO += HALFBG_SF\n";
# write smearing
$smearing && print $fh "MACRO += WITH_SMEARING\n";
# write clover
$clover && print $fh "MACRO += WITH_CLOVER\n";
# write expclover
$expclover && print $fh "MACRO += WITH_EXPCLOVER\n";
# write eo preconditioning
$eoprec && print $fh "MACRO += UPDATE_EO\n";
# write new geometry
$newgeo && print $fh "MACRO += WITH_NEW_GEOMETRY\n";
# write quaternions 
$quat && print $fh "MACRO += WITH_QUATERNIONS\n";
# write ndebug
$ndebug && print $fh "MACRO += NDEBUG\n";
# write dphi float
$dfloat && print $fh "MACRO += DPHI_FLT\n";
# write check spinor
$scheck && print $fh "MACRO += CHECK_SPINOR_MATCHING\n";
# write mpitiming
$mpit && print $fh "MACRO += MPI_TIMING\n";
# write io flush
$iof && print $fh "MACRO += IO_FLUSH\n";
# log all pids
$logallpids && print $fh "MACRO += LOG_ALLPIDS\n";
# write unroll representation
$unrollr && print $fh "MACRO += UNROLL_GROUP_REPRESENT\n";
# write timing
$timing && print $fh "MACRO += TIMING\n";
# write timing
$btiming && print $fh "MACRO += TIMING_WITH_BARRIERS\n";
# write memory
$mem && print $fh "MACRO += AMALLOC_MEASURE\n";
# write avx2
$avx2 && print $fh "MACRO += AVX2_HIREP\n";
# write vect
$vect && print $fh "MACRO += SIMD_VECTOR_HIREP\n";
# write fuse
$fuse && print $fh "MACRO += WITH_FUSE_MASTER_FOR\n";
# write probe mpi
$probempi && print $fh "MACRO += WITH_PROBE_MPI\n";
# write hwloc
$hwloc && print $fh "MACRO += HWLOC\n";
if ($hwloc) {
    $ldflags .=" -lhwloc";
}
# write mpi
$mpi && print $fh "MACRO += WITH_MPI\n";
# write GPU
$gpu && print $fh "MACRO += WITH_GPU\n";
# openMP
if($omp) {
    $cflags .=" -fopenmp";
    $ldflags .=" -fopenmp -latomic";
}

# write compiler options
if ($ccache!=0) { $cc="ccache ".$cc; $mpicc="ccache ".$mpicc; }
if ($color!=1) { print $fh "NOCOLOR = 1\n\n"; }
print $fh "ENV = $env\n";
print $fh "CC = $cc\n";
print $fh "MPICC = $mpicc\n";
print $fh "CFLAGS = $cflags\n";

print $fh "NVCC = $nvcc\n";
print $fh "GPUFLAGS = $gpuflags\n";

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
                                  P=Periodic
                                  A=Antiperiodic
                                  T=Theta
                                  O=Open
                                  S=Schrodinger functional
                                  SR=Rotated schrodinger functional
                                  M=Open + Schrodinger functional
  -t                  [P]         T boundary conditions (P, A, T, O ,S ,SR , M)
  -x                  [P]         X boundary conditions (P, A, T)
  -y                  [P]         Y boundary conditions (P, A, T)
  -z                  [P]         Z boundary conditions (P, A, T)

  --[no-]mpi          [true]      Use MPI
  --[no-]gpu          [false]     Use GPU acceleration
  --[no-]omp          [false]     Use OpenMP acceleration
  --env               []          Environment variables used for compilation
  --cc                [gcc]       Compiler
  --mpicc             [mpicc]     MPI Compiler
  --cflags            [-Wall -O3] Compilation options
  --nvcc              [nvcc]      CUDA compiler
  --gpuflags          []          CUDA compilation options
  --include           []          Extra include headers
  --ldflags           []          Linking options
  --[no-]ndebug       [true]      Set ndebug flag
  --[no-]ccache       [false]     Use ccache
  --[no-]color        [true]      Color compilation output

  --[no-]eo           [true]      Even-Odd preconditioning
  --[no-]newgeo       [false]     Use new geometry
  --[no-]fixedstr     [false]     Access memory on GPU using a fixed stride of length 32 instead of piece length
  
  --[no-]gauge-twist  [false]     Spatial gauge twist
  --[no-]sfhalf       [false]     Schrodinger functional half field
  --[no-]smearing     [false]     Smearing action
  --[no-]clover,-c    [false]     Clover improved action
  --[no-]expclover,-e [false]     ExpClover improved action

  --[no-]quat,-q      [false]     Use quaternion representation (only valid for SU2)
  --[no-]dfloat       [false]     Use single precision acceleration
  --[no-]unrollrepr   [false]     Unroll group representation functions

  --[no-]avx2         [false]     Enable avx2 kernels of (some) linear algebra functions
  --[no-]vect         [false]     Enable SIMD vectorized kernels of (some) linear algebra functions
  --[no-]fuse         [false]     Enable use of the fused index in the masterfor
  --[no-]probempi     [false]     Enable probing of mpi communications in the dirac operator
  --[no-]checkspinor  [true]      Check spinor field type
  --[no-]mpitiming    [false]     Enable timing of MPI calls
  --[no-]ioflush      [true]      Flush IO after each operations on logs
  --[no-]logallpids   [false]     Write log output for all MPI processes
  --[no-]timing       [false]     Enable timing
  --[no-]bartiming    [false]     Enable MPI barriers in timing
  --[no-]memory       [false]     Print memory usage
  --[no-]force        [false]     Print statistics for forces in molecular dynamics
  --[no-]hwloc        [false]     Improved process (and possibly gpu) affinity using hwloc

=head1 VERSION

1.0

=cut
