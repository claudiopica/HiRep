#!/usr/bin/perl -w
use strict;

## TODO: pl script docs clean up
## TODO: Possibly it is faster, if we do the cast only once (double*)
## Readout with assigning the complex from the doubles might not be ideal, because it 
## is a math operation internally while in reality it is just moving memory around.

# Read arguments from MkRules
# Ng ... Number of Colors
# rep ... Fermion Representation String Descriptor
# su2quat ... SU(2) with Quaternions true/false
# Gauge group ... SU(N) or SO(N), possible strings GAUGE_SUN, GAUGE_SON
my ($Ng,$rep,$su2quat,$gauge_group)=@ARGV;

# open STDOUT, ">gpu_geometry.h";
open STDOUT, ">strided_reads.h";


my $Nf = 0;
my $complex = "C"; # Need this for suN-matrices

if ($rep eq "REPR_FUNDAMENTAL") {
    $Nf=$Ng;
} elsif ($rep eq "REPR_SYMMETRIC") {
    $Nf=$Ng*($Ng+1)/2;
} elsif ($rep eq "REPR_ANTISYMMETRIC") {
    $Nf=$Ng*($Ng-1)/2;
} elsif ($rep eq "REPR_ADJOINT") {
    $Nf=$Ng*$Ng - 1;
    $complex="R";
} else {
    print "Fermion representation unspecified. Check Make/MkFlags.";
    exit(1);
}

## SO(N) are always real
if ($gauge_group eq "GAUGE_SON") {
    $complex="R";
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
        write_gpu_vector($prec, $repr);
        write_gpu_scalar($prec, $repr);
        write_gpu_suN($prec, $repr);
        write_gpu_suN_av($prec, $repr);
        write_gpu_clover_term($prec, $repr);
    }

    write_su2quat_redefinitions($prec);
}

write_gpu_ldl_field();

write_epilog();

sub write_prolog {
    print <<END
/*******************************************************************************
*
* File gpu_geometry.h
*
* Reading and writing according to the GPU geometry 
*
*******************************************************************************/

/**
 * \@file
 * \@brief Memory access patterns are crucial to achieve best performance.     
 *        Due to this, we store lattice field data differently in device memory
 *        than in host memory. Documentation on this can be found in the HiRep 
 *        Development Guide, section GPU Geometry.
*/
#ifndef GPU_GEOMETRY_H
#define GPU_GEOMETRY_H

#ifdef WITH_GPU
END
}

sub write_epilog {
    print "\n\n#endif\n#endif";
}

sub write_gpu_spinor {
    my ($prec, $repr) = @_;
    my @dim_vector = ($Ng, $Nf);
    my $i;
    my $comp;

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
    print " * \@param _s     \t\t${typename} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tComponent to read, choose 0 for spinor fields.\n";
    print " */\n";
    print "#define read_gpu_${typename}(_stride, _s, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\tint __iz = (_ix); \\\n";
    for ($comp=0; $comp<3; $comp++) {
        for ($i=0; $i<$N; $i++) {
            print "\t\t(_s).c\[$comp\].c\[$i\]=(($type*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
        }
    }
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t(_s).c\[$comp\].c\[$i\]=(($type*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
    }
    print "\t\t(_s).c\[$comp\].c\[$i\]=(($type*)(_in))\[__iz\]; \\\n";
    print "\t} while (0) \n\n";

    # Generate write macro
    print "/**\n";
    print " * \@brief Write ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _s     \t\t${component_type} target to write to the field _out\n";
    print " * \@param _out   \t\tInput field to write to\n";
    print " * \@param _ix    \t\tIndex at which to write \n";
    print " * \@param _comp  \t\tComponent to write, choose 0 for spinor fields. \n";
    print " */\n";
    print "#define write_gpu_${typename}(_stride, _s, _out, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\tint __iz = (_ix); \\\n";
    for ($comp=0; $comp<3; $comp++) {
        for ($i=0; $i<$N; $i++) {
            print "\t\t(($type*)(_out))\[__iz\]=(_s).c\[$comp\].c\[$i\]; __iz+=(_stride); \\\n";
        }
    }
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t(($type*)(_out))\[__iz\]=(_s).c\[$comp\].c\[$i\]; __iz+=(_stride); \\\n";
    }
    print "\t\t(($type*)(_out))\[__iz\]=(_s).c\[$comp\].c\[$i\]; \\\n";
    print "\t} while (0) \n\n";
}

sub write_gpu_vector {
    my ($prec, $repr) = @_;
    my @dim_vector = ($Ng, $Nf);
    my $i;

    # Generate basename for given representation
    my $dataname = $basename.$rep_suffixes[$repr]."_vector";
    
    # Generate precision suffix
    my $precision_suffix = $precision_suffixes[$prec];

    # Complete typename with suffixes
    my $typename = $dataname.$precision_suffix;

    # Complex vector components to be separated by stride
    # here we need elementary complex types.
    my $type = $precision_c_types[$prec];

    # representation dimension
    my $N = $dim_vector[$repr];

    # Generate read macro
    print "/**\n";
    print " * \@brief Read ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${typename} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tComponent of the ${typename} to read.\n";
    print " */\n";
    print "#define read_gpu_${typename}(_stride, _v, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
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
    print " * \@param _v     \t\t${typename} target to write to the field _out\n";
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

sub write_gpu_suN {
    my ($prec, $repr) = @_;
    my @dim_suN = ($Ng*$Ng, $Nf*$Nf);
    my $i;

    # Generate basename for given representation
    my $dataname = $basename.$rep_suffixes[$repr];
    
    # Generate precision suffix
    my $precision_suffix = $precision_suffixes[$prec];

    # Complete typename with suffixes
    my $typename = $dataname.$precision_suffix;

    # If the representation is real, then we want to separated the 
    # real components (double or float) by the stride otherwise, the components 
    # will be complex numbers, which occupy twice as much memory.
    my $type;
    if ($repr == 1) {
        if ($complex eq "C") {
            $type = $precision_c_types[$prec];
        } else {
            $type = $precision_types[$prec];
        }
    } else {
        $type = $precision_c_types[$prec];
    }

    #$type = $precision_c_types[$prec];
    my $N = $dim_suN[$repr];

    # Generate read macro
    print "/**\n";
    print " * \@brief Read ${typename} matrix according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${typename} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tLink direction to read.\n";
    print " */\n";
    print "#define read_gpu_${typename}(_stride, _v, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\tint __iz = (_ix) + ((_comp)*$N)*(_stride); \\\n";
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t(_v).c\[$i\]=((${type}*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
    }
    print "\t\t(_v).c\[$i\]=((${type}*)(_in))\[__iz\]; \\\n";
    print "\t} while (0) \n\n";

    # Generate write macro
    print "/**\n";
    print " * \@brief Write ${typename} matrix according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${typename} target to write to the field _out\n";
    print " * \@param _out   \t\tInput field to write to\n";
    print " * \@param _ix    \t\tIndex at which to write \n";
    print " * \@param _comp  \t\tLink direction to write. \n";
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

## TODO: Adjust this to work for real representations, where
## suN*c is not aliased to suN* but suN*_FMAT
sub write_gpu_clover_term {

    my ($prec, $repr) = @_;
    my @dim_suN = ($Ng*$Ng, $Nf*$Nf);
    my $N = $dim_suN[$repr];
    my $i;

    # Generate basename for given representation
    my $dataname = $basename.$rep_suffixes[$repr];
    my $alias_dataname = $dataname."c";
    
    # Generate precision suffix
    my $precision_suffix = $precision_suffixes[$prec];

    # Complete typename with suffixes
    my $typename = $dataname.$precision_suffix;
    my $typename_alias = $alias_dataname.$precision_suffix;
    my $type = $precision_c_types[$prec];

    if ($complex eq "C") { 
        print "/**\n";
        print " * \@brief Read ${typename_alias} according to device geometry structure \n";
        print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
        print " * \@param _v     \t\t${typename_alias} target to read to from the field _in\n";
        print " * \@param _in    \t\tInput field to read from \n";
        print " * \@param _ix    \t\tIndex at which to read \n";
        print " * \@param _comp  \t\tComponent to read for argument consistency between different GPU read/write functions.\\\n";
        print "                  \t\tUse this macro here always with _comp=0, because this is for a scalar field!\n";
        print "*/\n";
        print "#define read_gpu_${typename_alias}(_stride, _v, _in, _ix, _comp) \\\n";
        print "\t\t\tread_gpu_${typename}((_stride), (_v), (_in), (_ix), (_comp))\n";

        # Generate write macro
        print "/**\n";
        print " * \@brief Write ${typename_alias} according to device geometry structure \n";
        print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
        print " * \@param _v     \t\t${typename_alias} target to write to the field _out\n";
        print " * \@param _out   \t\tInput field to write to\n";
        print " * \@param _ix    \t\tIndex at which to write \n";
        print " * \@param _comp  \t\tComponent to write for argument consistency between different GPU read/write functions.\\\n";
        print "                  \t\tUse this macro here always with _comp=0, because this is for a scalar field!\n";
        print " */\n";
        print "#define write_gpu_${typename_alias}(_stride, _v, _out, _ix, _comp) \\\n";
        print "\t\t\twrite_gpu_${typename}((_stride), (_v), (_out), (_ix), (_comp))\n";
    } else {
        # Generate read macro
        print "/**\n";
        print " * \@brief Read ${typename_alias} matrix according to device geometry structure \n";
        print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
        print " * \@param _v     \t\t${typename_alias} target to read to from the field _in\n";
        print " * \@param _in    \t\tInput field to read from \n";
        print " * \@param _ix    \t\tIndex at which to read \n";
        print " * \@param _comp  \t\tLink direction to read.\n";
        print " */\n";
        print "#define read_gpu_${typename_alias}(_stride, _v, _in, _ix, _comp) \\\n";
        print "\tdo { \\\n";
        print "\t\tint __iz = (_ix) + ((_comp)*$N)*(_stride); \\\n";
        for ($i=0; $i<$N-1; $i++) {
            print "\t\t(_v).c\[$i\]=((${type}*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
        }
        print "\t\t(_v).c\[$i\]=((${type}*)(_in))\[__iz\]; \\\n";
        print "\t} while (0) \n\n";

        # Generate write macro
        print "/**\n";
        print " * \@brief Write ${typename_alias} matrix according to device geometry structure \n";
        print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
        print " * \@param _v     \t\t${typename_alias} target to write to the field _out\n";
        print " * \@param _out   \t\tInput field to write to\n";
        print " * \@param _ix    \t\tIndex at which to write \n";
        print " * \@param _comp  \t\tLink direction to write. \n";
        print " */\n";
        print "#define write_gpu_${typename_alias}(_stride, _v, _out, _ix, _comp) \\\n";
        print "\tdo { \\\n";
        print "\t\tint __iz = (_ix) + ((_comp)*$N)*(_stride); \\\n";
        for ($i=0; $i<$N-1; $i++) {
            print "\t\t((${type}*)(_out))\[__iz\]=(_v).c\[$i\]; __iz+=(_stride);\\\n";
        }
        print "\t\t((${type}*)(_out))\[__iz\]=(_v).c\[$i\]; \\\n";
        print "\t} while (0) \n\n";
    }
}

sub write_su2quat_redefinitions {
    my $prec = @_;
    if ($su2quat==1) {
        # Basenames for both representations
        my $dataname = $basename.$rep_suffixes[0];
        my $datanameR = $basename.$rep_suffixes[1];
        
        # Generate precision suffix
        my $precision_suffix = $precision_suffixes[$prec];

        # Complete typenames with suffixes
        my $typename = $dataname.$precision_suffix;
        my $typenameR = $datanameR.$precision_suffix;
    
        print "#define write_gpu_${typename}(_stride, _v, _out, _ix, _comp) ";
        print "write_gpu_${typenameR}(_stride, _v, _out, _ix, _comp)\n\n";
    }
}

sub write_gpu_suN_av {
    my ($prec, $repr) = @_;

    # Depending on gauge group type, different dimensions
    my @dim_vector = ($Ng*$Ng-1, $Nf*$Nf-1);
    if ($gauge_group eq "GAUGE_SON") {
        @dim_vector = ($Ng*($Ng-1)/2, $Nf*($Nf-1)/2);
    }
    my $i;
    my $N = $dim_vector[$repr];

    # Generate basename for given representation
    my $dataname = $basename.$rep_suffixes[$repr]."_algebra_vector";
    
    # Generate precision suffix
    my $precision_suffix = $precision_suffixes[$prec];

    # Complete typename with suffixes
    my $typename = $dataname.$precision_suffix;

    # Algebra vectors have real components -> non-c type
    my $type = $precision_types[$prec];

    # Generate read macro
    print "/**\n";
    print " * \@brief Read ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${typename} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tComponent of the ${typename} to read.\n";
    print " */\n";
    print "#define read_gpu_${typename}(_stride, _v, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
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
    print " * \@param _v     \t\t${typename} target to write to the field _out\n";
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

sub write_gpu_scalar {
    my ($prec, $repr) = @_;
    my @dim_vector = ($Ng, $Nf);
    my $i;

    # Complex vector components to be separated by stride
    # here we need elementary complex types.
    my $type = $precision_types[$prec];

    # Generate read macro
    print "/**\n";
    print " * \@brief Read ${type} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${type} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tComponent of the ${type} to read.\n";
    print " */\n";
    print "#define read_gpu_${type}(_stride, _v, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\t(_v)=*((_in)+(_ix));\\\n";
    print "\t} while (0) \n\n";

    # Generate write macro
    print "/**\n";
    print " * \@brief Write ${type} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${type} target to write to the field _out\n";
    print " * \@param _out   \t\tInput field to write to\n";
    print " * \@param _ix    \t\tIndex at which to write \n";
    print " * \@param _comp  \t\tComponent of the ${type} to write. \n";
    print " */\n";
    print "#define write_gpu_${type}(_stride, _v, _out, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\t(*((_out)+(_ix)))=(_v);\\\n";
    print "\t} while (0) \n\n";
}


sub write_gpu_ldl_field {
    my $i;
    my $N = $Nf * (2 * $Nf + 1);

    # Complete typename with suffixes
    my $typename = "ldl_t";
    my $type = "hr_complex";

    # Generate read macro
    print "/**\n";
    print " * \@brief Read ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${typename} target to read to from the field _in\n";
    print " * \@param _in    \t\tInput field to read from \n";
    print " * \@param _ix    \t\tIndex at which to read \n";
    print " * \@param _comp  \t\tComponent to read for consistency (put 0 for this type).\n";
    print " */\n";
    print "#define read_gpu_${typename}(_stride, _v, _in, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\tint __iz = (_ix) + ((_comp)*2*$N)*(_stride); \\\n";
    for ($i=0; $i<$N; $i++) {
        print "\t\t(_v).up\[$i\]=((${type}*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
    }
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t(_v).dn\[$i\]=((${type}*)(_in))\[__iz\]; __iz+=(_stride); \\\n";
    } 
    print "\t\t(_v).dn\[$i\]=((${type}*)(_in))\[__iz\]; \\\n";

    print "\t} while (0) \n\n";

    # Generate write macro
    print "/**\n";
    print " * \@brief Write ${typename} according to device geometry structure \n";
    print " * \@param _stride\t\tInteger valued stride with which the components are stored\n";
    print " * \@param _v     \t\t${typename} target to write to the field _out\n";
    print " * \@param _out   \t\tInput field to write to\n";
    print " * \@param _ix    \t\tIndex at which to write \n";
    print " * \@param _comp  \t\tComponent to write for consistency (put 0 for this type).\n";
    print " */\n";
    print "#define write_gpu_${typename}(_stride, _v, _out, _ix, _comp) \\\n";
    print "\tdo { \\\n";
    print "\t\tint __iz = (_ix) + ((_comp)*2*$N)*(_stride); \\\n";
    for ($i=0; $i<$N; $i++) {
        print "\t\t((${type}*)(_out))\[__iz\]=(_v).up\[$i\]; __iz+=(_stride);\\\n";
    }
    for ($i=0; $i<$N-1; $i++) {
        print "\t\t((${type}*)(_out))\[__iz\]=(_v).dn\[$i\]; __iz+=(_stride);\\\n";
    }
    print "\t\t((${type}*)(_out))\[__iz\]=(_v).dn\[$i\]; \\\n";
    print "\t} while (0) \n\n";
}
