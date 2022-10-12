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
my @precision_desc = ("Double Precision", "Single Precision");

# representations
my $basename = "suN";
my @rep_suffixes = ("g", "f");

### WRITING FILE CONTENTS
write_prolog();
my @a = (0,1);
for my $prec (@a) {
    for my $repr (@a) {
        write_gpu_spinor($prec, $repr);
    }
}

write_epilog();
###

## Write functions
#write_gpu_spinor();
#write_gpu_geometry_functions("_vector", @dim_vector);

#my @dim_suN = ($Ng*$Ng, $Nf*$Nf);
#write_gpu_geometry_functions("", @dim_suN); # suNg, suNf, ...
#write_gpu_geometry_functions("c", @dim_suN); # suNgc, suNfc, ...

#my @dim_algebra_vector = (($Ng*$Ng)-1, ($Nf*$Nf)-1);
#write_gpu_geometry_functions("_algebra_vector", @dim_algebra_vector);

#write_gpu_geometry_function_scalar();

sub write_prolog {
    print "/*******************************************************************************\n";
    print "*\n";
    print "* File gpu_geometry.h\n";
    print "*\n";
    print "* Reading and writing according to the GPU geometry \n";
    print "*\n";
    print "*******************************************************************************/\n";
    print "\n";
    print "/**\n";
    print " * \@brief Memory access patterns are crucial to achieve best performance.     \n";
    print " *        Due to this, we store lattice field data differently in device memory\n";
    print " *        than in host memory. Documentation on this can be found in the HiRep \n";
    print " *        Development Guide, section GPU Geometry.\n";
    print "*/\n\n";
    print "#ifndef GPU_GEOMETRY_H\n";
    print "#define GPU_GEOMETRY_H\n\n";
}

sub write_epilog {
    print "\n\n\n#endif";
}

sub write_gpu_spinor {
    my ($prec, $repr) = @_;
    my @dim_vector = ($Ng, $Nf);
    my $i;

    # Generate basename for given representation
    my $dataname = $basename.$rep_suffixes[$repr]."_spinor";
    
    # Generate precision suffix
    my $precision_suffix = $precision_suffixes[$prec];

    # Complete typename with suffixes
    my $typename = $dataname.$precision_suffix;
    my $component_type = $basename.$rep_suffixes[$repr]."_vector".$precision_suffix;

    # Complex vector components to be separated by stride
    # here we need elementary complex types.
    my $type = $precision_c_types[$prec];

    # representation dimension
    my $N = $dim_vector[$repr];

    # Generate read macro
    print "/**\n";
    print " * \@brief Read ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${component_type} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tSpinor component of the ${typename} to read.\n";
    print " */\n";
    print "#define read_gpu_${typename}(_stride, _v, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\t${type} real_pt, imag_pt; \\\n";
    print "\t\tint __iz = (_ix) + ((_comp)*$N)*(_stride); \\\n";
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t(_v).c\[$i\]=((${type}*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
    }
    print "\t\t(_v).c\[$i\]=((${type}*)(_in))\[__iz\]; \\\n";
    print "\t} while (0) \n\n";

    # Generate write macro
    print "/**\n";
    print " * \@brief Write ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${component_type} target to write to the field _out\n";
    print " * \@param _out   \t\tInput field to write to\n";
    print " * \@param _ix    \t\tIndex at which to write \n";
    print " * \@param _comp  \t\tComponent of the ${typename} to write. \n";
    print " */\n";
    print "#define write_gpu_${typename}(_stride, _v, _out, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\tint __iz = (_ix) + ((_comp)*$N)*(_stride); \\\n";
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t((${type}*)(_out))\[__iz\]=(_v).c\[$i\]; __iz+=(_stride);\\\n";
    }
    print "\t\t((${type}*)(_out))\[__iz\]=(_v).c\[$i\]; \\\n";
    print "\t} while (0) \n\n";
}