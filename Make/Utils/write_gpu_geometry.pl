#!/usr/bin/perl -w
use strict;

## TODO: Prolog, Epilog, Macro docs, pl script docs
## Problems: Spinor field layers, scalar fields need separate def
## remove stride, because it is always half number of lattice points->just compile with this
## TODO: Possibly it is faster, if we do the cast only once (double*)
## Readout with assigning the complex from the doubles might not be ideal, because it 
## is a math operation internally while in reality it is just moving memory around.

# Read arguments from MkRules
my ($Ng,$rep,$su2quat,$gauge_group)=@ARGV;

open STDOUT, ">gpu_geometry.h";

my $Nf = 0;

if ($rep eq "REPR_FUNDAMENTAL") {
    $Nf=$Ng;
} elsif ($rep eq "REPR_SYMMETRIC") {
    $Nf=$Ng*($Ng+1)/2;
} elsif ($rep eq "REPR_ANTISYMMETRIC") {
    $Nf=$Ng*($Ng-1)/2;
} elsif ($rep eq "REPR_ADJOINT") {
    $Nf=$Ng*$Ng - 1;
} else {
    print "Fermion representation unspecified. Check Make/MkFlags.";
    exit(1);
}

## Single precision/double precision
my @precision_suffixes = ("", "_flt");
my @precision_types = ("double", "float");
my @precision_c_types = ("hr_complex", "hr_complex_flt");

# representations
my $basename = "suN";
my @rep_suffixes = ("g", "f");

## Write functions
my @dim_vector = ($Ng, $Nf);
write_gpu_geometry_functions("_spinor", @dim_vector); # this is the same as vector, but has more components (4)
write_gpu_geometry_functions("_vector", @dim_vector);

my @dim_suN = ($Ng*$Ng, $Nf*$Nf);
write_gpu_geometry_functions("", @dim_suN); # suNg, suNf, ...
write_gpu_geometry_functions("c", @dim_suN); # suNgc, suNfc, ...

my @dim_algebra_vector = (($Ng*$Ng)-1, ($Nf*$Nf)-1);
write_gpu_geometry_functions("_algebra_vector", @dim_algebra_vector);

write_gpu_geometry_function_scalar();

sub write_gpu_geometry_functions {
    my ($site_element, @size) = @_;
    my @a = (0,1);
    my $i;
    for my $prec (@a) { # generate for different precisions
        for my $repr (@a) { # generate for different representations
            my $dataname = $basename.$rep_suffixes[$repr];
            my $precision_suffix = $precision_suffixes[$prec];
            my $type = $precision_types[$prec];
            my $N = 2*$size[$repr]; # factor of two due to complex

            my $typename = $dataname.$site_element.$precision_suffix;

            print "#define read_gpu_${typename}(stride, v, out, ix, comp) \\\n";
            print "\tdo { \\\n";
            print "\t\t${type} real_pt, imag_pt; \\\n";
            print "\t\tint __iz = (ix) + ((comp)*$N)*(stride); \\\n";
            for ($i=0; $i<($N-2)/2; $i++) {
                print "\t\treal_pt = ((${type}*)(in))\[__iz\]; __iz+=(stride); \\\n";
                print "\t\timag_pt = ((${type}*)(in))\[__iz\]; __iz+=(stride); \\\n";
                print "\t\t(v).c\[$i\]=hr_complex(real_pt, imag_pt); \\\n\t\t\\\n";
            }
                print "\t\treal_pt = ((${type}*)(in))\[__iz\]; __iz+=(stride); \\\n";
                print "\t\timag_pt = ((${type}*)(in))\[__iz\]; \\\n";
                print "\t\t(v).c\[$i\]=hr_complex(real_pt, imag_pt); \\\n";
            print "\t} while (0) \n\n";

            print "#define write_gpu_${typename}(stride, v, out, ix, comp) \\\n";
            print "\tdo { \\\n";
            print "\t\tint __iz = (ix) + ((comp)*$N)*(stride); \\\n";
            for ($i=0; $i<($N-2)/2; $i++) {
                print "\t\t((${type}*)(out))\[__iz\]=_complex_re((v).c\[$i\]); __iz+=(stride);\\\n";
                print "\t\t((${type}*)(out))\[__iz\]=_complex_im((v).c\[$i\]); __iz+=(stride);\\\n";
            }
            print "\t\t((${type}*)(out))\[__iz\]=_complex_re((v).c\[$i\]); __iz+=(stride);\\\n";
            print "\t\t((${type}*)(out))\[__iz\]=_complex_im((v).c\[$i\]);\\\n";
            print "\t} while (0) \n\n";
        }
    }
}

sub write_gpu_geometry_function_scalar {
    my @a = (0,1);
    my $i;
    for my $prec (@a) {
        my $typename = $precision_types[$prec];

        print "#define read_gpu_${typename}(stride, v, out, ix, comp) \\\n";
        print "\tdo { \\\n";
        print "\t\t (v)=(in)[ix]; \\\n";
        print "\t} while (0) \n\n";

        print "#define write_gpu_${typename}(stride, v, out, ix, comp) \\\n";
        print "\tdo { \\\n";
        print "\t\t (in)[ix]=(v); \\\n";
        print "\t} while (0) \n\n";
    }
}


