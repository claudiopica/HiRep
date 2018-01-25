#!/usr/bin/perl -w
use strict;

(@ARGV==2 or @ARGV==3 or @ARGV==4) or die("Usage: $0 Ng rep [su2_quaternion gauge_group]\nsu2 quaternion: 0 = 2x2 complex matrix, 1 = 4 reals\n");

my ($Nmax,$unroll)=(5,4);
my ($vd,$vr); #for vectors
my ($avd,$avr); #for algebra vectors
my ($md,$mr); #for matrices
my ($md2,$mr2); #for matrices


my ($Ng,$rep,$su2quat,$gauge_group)=@ARGV;
if (not defined ($gauge_group)) { $gauge_group = "SUN"; }
if (not defined ($su2quat)) { $su2quat = "0"; }
if (!($su2quat eq "0") and !($su2quat eq "1")) {
   die("Invalid option for su2 quaternion [$su2quat] specification. Exiting...\n");
}
if (!($su2quat eq "0") and $Ng!=2) {
   die("su2 quaternion option can only be used with Ng=2. Exiting...\n");
}

my ($Nf,$c1,$c2);
$c1="C"; #default gauge field complex
if ($rep eq "REPR_FUNDAMENTAL") {
	$Nf=$Ng;
	$c2="C";
} elsif ($rep eq "REPR_SYMMETRIC") {
	$Nf=$Ng*($Ng+1)/2;
	$c2="C";
} elsif ($rep eq "REPR_ANTISYMMETRIC") {
	$Nf=$Ng*($Ng-1)/2;
	$c2="C";
} elsif ($rep eq "REPR_ADJOINT") {
	$Nf=$Ng*$Ng-1;
	$c2="R";
} else {
	print "Please specify:\nFUN => fundamental rep\nSYM => symmetric rep\nASY => antisymmetric rep\nADJ => adjoint rep\n";
	exit(1);
}

#debug
#print "Ng=$Ng complex=$c1\n";
#print "Nf=$Nf complex=$c2\n";

#check for SO(N) gauge group
if ($gauge_group eq "GAUGE_SON"){ #all represenatations real
    $c1="R";
    $c2="R";
}

my ($N,$suff,$complex,$to);

my $dataname;
my $rdataname;
my $cname="c";
my $structdef="typedef struct\n{\n";

my ($basename,$fundsuff,$repsuff)=("suN","g","f"); #basename for types and suffix for fundamental representation


open STDOUT, ">suN_types.h";

write_prolog_suN_types();

print "#define NG $Ng\n";
#system("./write_suN_def.pl $Ng g $c1 T");
write_suN_h($Ng,$fundsuff,$c1,"T");

print "#define NF $Nf\n";
#system("./write_suN_def.pl $Nf f $c2 T");
write_suN_h($Nf,$repsuff,$c2,"T");

write_epilog();

open STDOUT, ">suN.h";

write_prolog_suN();

#system("./write_suN_def.pl $Ng g $c1 O");
write_suN_h($Ng,$fundsuff,$c1,"O");

#system("./write_suN_def.pl $Nf f $c2 O");
write_suN_h($Nf,$repsuff,$c2,"O");

write_epilog();


#end main program

sub write_prolog_suN {
  print <<END
/*******************************************************************************
*
* File suN.h
*
* Type definitions and macros for SU(N) matrices and spinors  
*
*******************************************************************************/

#ifndef SUN_H
#define SUN_H

#include "suN_types.h"

END
}

sub write_prolog_suN_types {
  print <<END
/*******************************************************************************
*
* File suN_types.h
*
* Type definitions for SU(N) matrices and spinors  
*
*******************************************************************************/

#ifndef SUN_TYPES_H
#define SUN_TYPES_H

#include "complex.h"

END
}

sub write_epilog {
  print <<END

#endif
END
}

sub write_suN_h {

($N,$suff,$complex,$to)=@_;
($vd,$vr)=(int($N/$unroll)*$unroll,$N%$unroll);
($avd,$avr)=(int(($N*$N-1)/(2*$unroll))*2*$unroll,($N*$N-1)%(2*$unroll));
($md,$mr)=(int(($N-1)/$unroll)*$unroll,($N-1)%$unroll);
($md2,$mr2)=(int(($N*$N)/$unroll)*$unroll,($N*$N)%$unroll);

if ((not ($complex eq "C")) and (not ($complex eq "R"))) {
    die("Error: type must be C or R!\n");
}
if ((not ($to eq "T")) and (not ($to eq "O"))) {
    die("Error: Specify T for data types or O for operations!\n");
}

$dataname="$basename$suff";
$rdataname=$dataname;
if ($complex eq "R") {
    $dataname.="c";
}

if ($to eq "T") {

print <<END
/*******************************************************************************
*
* Definitions of Data Structures
*
*******************************************************************************/

END
;
write_suN_vector();

if ($su2quat==0) {
	write_suN();
    if ($complex eq "R") {
	    write_suNr();
		print "typedef $rdataname ${rdataname}_FMAT;\n\n";
		print "typedef ${rdataname}_flt ${rdataname}_FMAT_flt;\n\n";
    } else {
		print "typedef $dataname ${dataname}c;\n\n";
		print "typedef ${dataname}_flt ${dataname}c_flt;\n\n";
		print "typedef $dataname ${dataname}_FMAT;\n\n";
		print "typedef ${dataname}_flt ${dataname}_FMAT_flt;\n\n";
	}
	
} else {
    write_su2($su2quat);
    my ($ldn,$lrdn)=($dataname,$rdataname);
    $dataname="${rdataname}c";
	$rdataname="${rdataname}_FMAT";
    write_suN();
    if ($complex eq "R") {
        write_suNr();
    } else {
		print "typedef $dataname ${rdataname};\n\n";
		print "typedef ${dataname}_flt ${rdataname}_flt;\n\n";
	}
	$dataname=$ldn;
    $rdataname=$lrdn;
    
}



write_spinor();
if ($suff eq $fundsuff) { #algebra operations only for gauge
	write_suN_algebra_vector();
}

} else {

print <<END
/*******************************************************************************
*
* The following macros are the same for single and double precision types
*
* Depending on the macro, arguments are variables of type suN_vector and suN
* (or suN_vector_flt and suN_flt)
*
*******************************************************************************/

END
;
##write_vector_copy();
write_vector_zero();
write_vector_minus();
write_vector_i_plus();
write_vector_i_minus();
write_vector_mul();
write_vector_mulc();
write_vector_mulc_star();
write_vector_add();
write_vector_sub();
write_vector_i_add();
write_vector_i_sub();
write_vector_add_assign();
write_vector_sub_assign();
write_vector_i_add_assign();
write_vector_i_sub_assign();
write_vector_prod_re();
write_vector_prod_im();
write_vector_mulc_add_assign();
write_vector_mul_add_assign();
write_vector_lc();
write_vector_lc_add_assign();
write_vector_clc();
write_vector_clc_add_assign();
write_vector_prod_assign();
write_vector_prod_add_assign_re();
write_vector_prod_add_assign_im();
write_vector_prod_sub_assign_re();
write_vector_prod_sub_assign_im();

write_vector_project();

if ($su2quat==0) {
  write_suN_multiply();
  write_suN_inverse_multiply();
  write_suN_zero();
  if ($complex eq "R") {
    write_suNr_multiply();
    write_suNr_inverse_multiply();
 	write_suNr_zero();
  } else {
	  print "#define _suNfc_multiply(a,b,c) _suNf_multiply(a,b,c)\n\n";
	  print "#define _suNfc_inverse_multiply(a,b,c) _suNf_inverse_multiply(a,b,c)\n\n";
	  print "#define _suNfc_zero(a) _suNf_zero(a)\n\n";
  }
} else {
    #write_su2_decode($su2quat);
    write_su2_multiply();
    write_su2_inverse_multiply();
	write_su2_zero();
 
    my ($ldn,$lrdn)=($dataname,$rdataname);
  	if ($complex eq "C") {
    	$dataname="${dataname}c";
  	}
    write_suN_multiply();
    write_suN_inverse_multiply();
	write_suN_zero();
    $dataname=$ldn;
}

if ($suff eq "g") { #algebra operations only for gauge
	write_algebra_vector_add_assign();
	write_algebra_vector_sub_assign();
	write_algebra_vector_mul_add_assign();
	write_algebra_vector_mul();
	write_algebra_vector_zero();
	write_algebra_vector_sqnorm();
}

print <<END
/*******************************************************************************
*
* Macros for SU(N) matrices
*
* Arguments are variables of type suN
*
*******************************************************************************/

END
;

if ($su2quat==0) {
  write_suN_dagger();
  write_suN_times_suN();
  write_suN_times_suN_dagger();
  write_suN_dagger_times_suN();
 ## write_suN_zero();
  write_suN_unit();
  write_suN_minus();
# write_suN_copy();
  write_suN_mul();
  write_suN_mulc();
  write_suN_add_assign();
  write_suN_sub_assign();
  write_suN_sqnorm();
  write_suN_sqnorm_m1();
  write_suN_trace_re();
  write_suN_trace_im();
  #write_suN_2TA();
  #write_suN_TA();
  write_suN_FMAT();

  if ($complex eq "R") { # we only need these functions at the moment...
 ##   write_suNr_zero();
    write_suNr_FMAT();
    write_suNr_unit();
    write_suNr_dagger();
    write_suNr_times_suNr();
    write_suNr_times_suNr_dagger();
    write_suNr_dagger_times_suNr();
    write_suNr_add_assign();
    write_suNr_sub_assign();
    write_suNr_mul();
    write_suNr_trace_re();
    write_suNr_sqnorm();
    write_suNr_minus();
  }

} else {
    write_su2_dagger();
    write_su2_times_su2();
    write_su2_times_su2_dagger();
    write_su2_dagger_times_su2();
##	write_su2_zero();
    write_su2_unit();
    write_su2_minus();
    write_su2_mul();
    write_su2_add_assign();
    write_su2_sub_assign();
    write_su2_sqnorm();
    write_su2_sqnorm_m1();
    write_su2_trace_re();
    write_su2_trace_im();
    #write_su2_2TA();
    #write_su2_TA();
    write_suN_FMAT(); #this is the same as before
    if ($complex eq "R") {
        write_suNr_FMAT(); #this is the same as before
    }
    
    write_su2_exp();
}

my ($ldn,$lrdn)=($dataname,$rdataname);
$dataname="${dataname}_FMAT";
$rdataname="${rdataname}_FMAT";
write_suN_zero();
if ($complex eq "R") { # we only need these functions at the moment...
	write_suNr_zero();
}
$dataname=$ldn;
$rdataname=$lrdn;

print <<END
/*******************************************************************************
*
* Macros for spinors
*
* Arguments are variables of type spinors
*
*******************************************************************************/

END
;
#write_spinor_copy();
write_spinor_zero();
write_spinor_g5();
write_spinor_minus();
write_spinor_mul();
write_spinor_mulc();
write_spinor_mulc_add_assign();
write_spinor_mul_add_assign();
write_spinor_lc();
write_spinor_lc_add_assign();
write_spinor_clc();
write_spinor_clc_add_assign();
write_spinor_add();
write_spinor_sub();
write_spinor_add_assign();
write_spinor_sub_assign();
write_spinor_i_add_assign();
write_spinor_i_sub_assign();
write_spinor_prod_re();
write_spinor_prod_im();
write_spinor_prod();
write_spinor_prod_assign();
write_spinor_g5_prod_re();
write_spinor_g5_prod_im();
write_spinor_project();
write_spinor_pminus();
write_spinor_pplus();
    
    #GPU READ/WRITE Functions
    write_read_spinor_gpu();    
    write_write_spinor_gpu();
    write_suN_av_read_gpu();
    write_suN_av_write_gpu();
    write_suN_av_mul_add_assign_gpu();
    
    if ($su2quat==0) {
        if ($complex eq "R") {
            write_suNr_read_gpu();
            write_suNr_write_gpu();
        } 
        write_suN_read_gpu();
        write_suN_write_gpu();
    } else {
        write_su2_read_gpu();
        write_su2_write_gpu();
    }
    
# COMMENTATO
# print <<END
# /*******************************************************************************
# *
# * Debug utilities
# *
# *******************************************************************************/

# END
# ;
# write_vector_myrandom();
# write_vector_iszero();
# write_su3_myrandom();



#write_epilog();

}

} #end of write_suN

#
# DATA TYPES
#

sub write_suN_vector {
  print $structdef;
  print "   complex $cname\[$N\];\n";
  print "} ${rdataname}_vector;\n\n";
  print $structdef;
  print "   complex_flt $cname\[$N\];\n";
  print "} ${rdataname}_vector_flt;\n\n";
}

sub write_suN_algebra_vector {
  print $structdef;
  my $d=($N*$N)-1;
  if ($gauge_group eq "GAUGE_SON"){ #N*(N-1)/2 generators
      $d=$N*($N-1)/2;
  }
  print "   double $cname\[$d\];\n";
  print "} ${rdataname}_algebra_vector;\n\n";
  print $structdef;
  print "   float $cname\[$d\];\n";
  print "} ${rdataname}_algebra_vector_flt;\n\n";
}

sub write_suN {
  print $structdef;
	my $d=($N*$N);
  print "   complex $cname\[$d\];\n";
  print "} $dataname;\n\n";
  print $structdef;
  print "   complex_flt $cname\[$d\];\n";
  print "} ${dataname}_flt;\n\n";
}

sub write_su2 {
  my ($su2repr)=@_;
  my $d;
  if ($su2repr==1) {
    $d=4
  } elsif ($su2repr==2) {
    $d=3;
  } else {
    die("Unknown su2 quaternionic form. Exiting...\n");
  }

  print $structdef;
  print "   double $cname\[$d\];\n";
  print "} $rdataname;\n\n";
  print $structdef;
  print "   float $cname\[$d\];\n";
  print "} ${rdataname}_flt;\n\n";
}

sub write_suNr {
  print $structdef;
	my $d=($N*$N);
  print "   double $cname\[$d\];\n";
  print "} $rdataname;\n\n";
  print $structdef;
  print "   float $cname\[$d\];\n";
  print "} ${rdataname}_flt;\n\n";
}

sub write_spinor {
  print $structdef;
  my $slen=4;
  print "   ${rdataname}_vector $cname\[$slen\];\n";
  print "} ${rdataname}_spinor;\n\n";
  print $structdef;
  print "   ${rdataname}_vector_flt $cname\[$slen\];\n";
  print "} ${rdataname}_spinor_flt;\n\n";
}

#
# VECTOR OPERATIONS
#

sub write_vector_copy {
  print "/* r=s */\n";
  print "#define _vector_copy_${suff}(r,s) \\\n";
  print "   (r)=(s)\n\n";
}

sub write_vector_zero {
  print "/* r=0 */\n";
  print "#define _vector_zero_${suff}(r) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_0((r).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_0((r).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_0((r).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_minus {
  print "/* r=-s */\n";
  print "#define _vector_minus_${suff}(r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all
		for(my $i=0;$i<$N;$i++){
			print "   _complex_minus((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_minus((r).$cname\[_i\],(s).$cname\[_i\]); ++_i; \\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_minus((r).$cname\[_i\],(s).$cname\[_i\]); ++_i; \\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_i_plus {
  print "/* r= i*s */\n";
  print "#define _vector_i_plus_${suff}(r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all
		for(my $i=0;$i<$N;$i++){
			print "   _complex_i_plus((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_i_plus((r).$cname\[_i\],(s).$cname\[_i\]); ++_i; \\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_i_plus((r).$cname\[_i\],(s).$cname\[_i\]); ++_i; \\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_i_minus {
    print "/* r=-i*s */\n";
    print "#define _vector_i_minus_${suff}(r,s) \\\n";
    if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all
	for(my $i=0;$i<$N;$i++){
	    print "   _complex_i_minus((r).$cname\[$i\],(s).$cname\[$i\])";
	    if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
	}
    } else { #partial unroll
	print "   do { \\\n";
	print "      int _i;for (_i=0; _i<$vd; ){\\\n";
	for(my $i=0;$i<$unroll;$i++){
	    print "         _complex_i_minus((r).$cname\[_i\],(s).$cname\[_i\]); ++_i; \\\n";
	}
	print "      }\\\n";
	for(my $i=0;$i<$vr;$i++){
	    print "      _complex_i_minus((r).$cname\[_i\],(s).$cname\[_i\]); ++_i; \\\n";
	}
	print "   } while(0) \n\n";
    }
}


sub write_vector_mul {
  print "/* r=k*s (k real) */\n";
  print "#define _vector_mul_${suff}(r,k,s) \\\n";
	if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all
		for(my $i=0;$i<$N;$i++){
			print "   _complex_mulr((r).$cname\[$i\],(k),(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mulr((r).$cname\[_i\],(k),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_mulr((r).$cname\[_i\],(k),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_mulc {
  print "/* r=z*s (z complex) */\n";
  print "#define _vector_mulc_${suff}(r,z,s) \\\n";
	if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all
		for(my $i=0;$i<$N;$i++){
			print "   _complex_mul((r).$cname\[$i\],(z),(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mul((r).$cname\[_i\],(z),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_mul((r).$cname\[_i\],(z),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_mulc_star {
    print "/* r=(z^+)*s (z complex) */\n";
    print "#define _vector_mulc_star_${suff}(r,z,s) \\\n";
	if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all
		for(my $i=0;$i<$N;$i++){
			print "   _complex_mul_star((r).$cname\[$i\],(s).$cname\[$i\],(z))";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mul_star((r).$cname\[_i\],(s).$cname\[_i\],(z)); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_mul_star((r).$cname\[_i\],(s).$cname\[_i\],(z)); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_add {
  print "/* r=s1+s2 */\n";
  print "#define _vector_add_${suff}(r,s1,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_add((r).$cname\[$i\],(s1).$cname\[$i\],(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_add((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_add((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_sub {
  print "/* r=s1-s2 */\n";
  print "#define _vector_sub_${suff}(r,s1,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_sub((r).$cname\[$i\],(s1).$cname\[$i\],(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_sub((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_sub((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_i_add {
  print "/* r=s1+i*s2 */\n";
  print "#define _vector_i_add_${suff}(r,s1,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_i_add((r).$cname\[$i\],(s1).$cname\[$i\],(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_i_add((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_i_add((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_i_sub {
  print "/* r=s1-i*s2 */\n";
  print "#define _vector_i_sub_${suff}(r,s1,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_i_sub((r).$cname\[$i\],(s1).$cname\[$i\],(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_i_sub((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_i_sub((r).$cname\[_i\],(s1).$cname\[_i\],(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_add_assign {
  print "/* r+=s */\n";
  print "#define _vector_add_assign_${suff}(r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_add_assign((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_add_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_add_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_sub_assign {
  print "/* r-=s */\n";
  print "#define _vector_sub_assign_${suff}(r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_sub_assign((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_sub_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_sub_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_i_add_assign {
  print "/* r+=i*s */\n";
  print "#define _vector_i_add_assign_${suff}(r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_i_add_assign((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_i_add_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_i_add_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_i_sub_assign {
  print "/* r-=i*s */\n";
  print "#define _vector_i_sub_assign_${suff}(r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_i_sub_assign((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_i_sub_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_i_sub_assign((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_prod_re {
	print "/* k=Re(r^*s) */\n";
  print "#define _vector_prod_re_${suff}(k,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)=_complex_prod_re((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		for(my $i=1;$i<$N;$i++){
			print "   (k)+=_complex_prod_re((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)=_complex_prod_re((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		print "      for (_i=1; _i<$md; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)+=_complex_prod_re((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr;$i++){
			print "      (k)+=_complex_prod_re((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_prod_im {
  print "/* k=Im(r*s) */\n";
  print "#define _vector_prod_im_${suff}(k,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)=_complex_prod_im((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		for(my $i=1;$i<$N;$i++){
			print "   (k)+=_complex_prod_im((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)=_complex_prod_im((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		print "      for (_i=1; _i<$md; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)+=_complex_prod_im((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr;$i++){
			print "      (k)+=_complex_prod_im((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}



sub write_vector_prod_add_assign_re {
	print "/* k+=Re(r^*s) */\n";
  print "#define _vector_prod_add_assign_re_${suff}(k,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)+=_complex_prod_re((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		for(my $i=1;$i<$N;$i++){
			print "   (k)+=_complex_prod_re((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)+=_complex_prod_re((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		print "      for (_i=1; _i<$md; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)+=_complex_prod_re((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr;$i++){
			print "      (k)+=_complex_prod_re((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}


sub write_vector_prod_add_assign_im {
  print "/* k+=Im(r*s) */\n";
  print "#define _vector_prod_add_assign_im_${suff}(k,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)+=_complex_prod_im((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		for(my $i=1;$i<$N;$i++){
			print "   (k)+=_complex_prod_im((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)+=_complex_prod_im((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		print "      for (_i=1; _i<$md; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)+=_complex_prod_im((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr;$i++){
			print "      (k)+=_complex_prod_im((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_prod_sub_assign_re {
	print "/* k-=Re(r^*s) */\n";
  print "#define _vector_prod_sub_assign_re_${suff}(k,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)-=_complex_prod_re((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		for(my $i=1;$i<$N;$i++){
			print "   (k)-=_complex_prod_re((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)-=_complex_prod_re((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		print "      for (_i=1; _i<$md; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)-=_complex_prod_re((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr;$i++){
			print "      (k)-=_complex_prod_re((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}


sub write_vector_prod_sub_assign_im {
  print "/* k-=Im(r*s) */\n";
  print "#define _vector_prod_sub_assign_im_${suff}(k,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)-=_complex_prod_im((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		for(my $i=1;$i<$N;$i++){
			print "   (k)-=_complex_prod_im((r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)-=_complex_prod_im((r).$cname\[0\],(s).$cname\[0\]);\\\n";
		print "      for (_i=1; _i<$md; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)-=_complex_prod_im((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr;$i++){
			print "      (k)-=_complex_prod_im((r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_mulc_add_assign {
  print "/* r+=z*s (z complex) */\n";
  print "#define _vector_mulc_add_assign_${suff}(r,z,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_mul_assign((r).$cname\[$i\],(z),(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mul_assign((r).$cname\[_i\],(z),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_mul_assign((r).$cname\[_i\],(z),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_mul_add_assign {
  print "/* r+=k*s (k real) */\n";
  print "#define _vector_mul_add_assign_${suff}(r,k,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_mulr_assign((r).$cname\[$i\],(k),(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mulr_assign((r).$cname\[_i\],(k),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_mulr_assign((r).$cname\[_i\],(k),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}


sub write_algebra_vector_add_assign {
  print "/* r+=s */\n";
  print "#define _algebra_vector_add_assign_${suff}(r,s) \\\n";
  my $d=$N*$N-1;
  if ($gauge_group eq "GAUGE_SON"){ #N(N-1)/2 generators
      $d=$N*($N-1)/2;
  }
  if ($N<$Nmax or $d<(4*$unroll+1) ) { #unroll all 
                for(my $i=0;$i<$d;$i++){
                        print "      (r).$cname\[$i\]+=(s).$cname\[$i\]";
                        if($i==$d-1) { print "\n\n"; } else { print "; \\\n"; }
                }
        } else { #partial unroll
                print "   do { \\\n";
                print "      int _i;for (_i=0; _i<$avd; ){\\\n";
                for(my $i=0;$i<2*$unroll;$i++){
                        print "         (r).$cname\[_i\]+=(s).$cname\[_i\]; ++_i;\\\n";
                }
                print "      }\\\n";
                for(my $i=0;$i<$avr;$i++){
                        print "      (r).$cname\[_i\]+=(s).$cname\[_i\]; ++_i;\\\n";
                }
                print "   } while(0) \n\n";
        }
}

sub write_algebra_vector_sub_assign {
  print "/* r-=s */\n";
  print "#define _algebra_vector_sub_assign_${suff}(r,s) \\\n";
  my $d=$N*$N-1;
  if ($gauge_group eq "GAUGE_SON"){ #N(N-1)/2 generators
      $d=$N*($N-1)/2;
  }
  if ($N<$Nmax or $d<(4*$unroll+1) ) { #unroll all 
                for(my $i=0;$i<$d;$i++){
                        print "      (r).$cname\[$i\]-=(s).$cname\[$i\]";
                        if($i==$d-1) { print "\n\n"; } else { print "; \\\n"; }
                }
        } else { #partial unroll
                print "   do { \\\n";
                print "      int _i;for (_i=0; _i<$avd; ){\\\n";
                for(my $i=0;$i<2*$unroll;$i++){
                        print "         (r).$cname\[_i\]-=(s).$cname\[_i\]; ++_i;\\\n";
                }
                print "      }\\\n";
                for(my $i=0;$i<$avr;$i++){
                        print "      (r).$cname\[_i\]-=(s).$cname\[_i\]; ++_i;\\\n";
                }
                print "   } while(0) \n\n";
        }
}

sub write_algebra_vector_mul_add_assign {
  print "/* r+=k*s (k real) */\n";
  print "#define _algebra_vector_mul_add_assign_${suff}(r,k,s) \\\n";
  my $d=$N*$N-1;
  if ($gauge_group eq "GAUGE_SON"){ #N(N-1)/2 generators
      $d=$N*($N-1)/2;
  }
  if ($N<$Nmax or $d<(4*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$d;$i++){
			print "      (r).$cname\[$i\]+=(k)*(s).$cname\[$i\]";
			if($i==$d-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$avd; ){\\\n";
		for(my $i=0;$i<2*$unroll;$i++){
			print "         (r).$cname\[_i\]+=(k)*(s).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$avr;$i++){
			print "      (r).$cname\[_i\]+=(k)*(s).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_algebra_vector_mul {
  print "/* r=k*s (k real) */\n";
  print "#define _algebra_vector_mul_${suff}(r,k,s) \\\n";
  my $d=$N*$N-1;
  if ($gauge_group eq "GAUGE_SON"){ #N(N-1)/2 generators
      $d=$N*($N-1)/2;
  }
  if ($N<$Nmax or $d<(4*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$d;$i++){
			print "      (r).$cname\[$i\]=(k)*(s).$cname\[$i\]";
			if($i==$d-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$avd; ){\\\n";
		for(my $i=0;$i<2*$unroll;$i++){
			print "         (r).$cname\[_i\]=(k)*(s).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$avr;$i++){
			print "      (r).$cname\[_i\]=(k)*(s).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_algebra_vector_zero {
  print "/* r=0  */\n";
  print "#define _algebra_vector_zero_${suff}(r) \\\n";
  my $d=$N*$N-1;
  if ($gauge_group eq "GAUGE_SON"){ #N(N-1)/2 generators
      $d=$N*($N-1)/2;
  }
  if ($N<$Nmax or $d<(4*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$d;$i++){
			print "      (r).$cname\[$i\]=0.";
			if($i==$d-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$avd; ){\\\n";
		for(my $i=0;$i<2*$unroll;$i++){
			print "         (r).$cname\[_i\]=0.; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$avr;$i++){
			print "      (r).$cname\[_i\]=0.; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_algebra_vector_sqnorm {
  print "/* k=|v|^2  */\n";
  print "#define _algebra_vector_sqnorm_${suff}(k,r) \\\n";
  my $last=$N*$N-1;
  if ($gauge_group eq "GAUGE_SON"){ #N(N-1)/2 generators
      $last=$N*($N-1)/2;
  }
  if ($N<$Nmax or $last<(4*$unroll+1) ) { #unroll all 
		print "   (k)=";
		my $n=0;
		for(my $i=0;$i<$last;$i++){
			print "((r).$cname\[$i\]*(r).$cname\[$i\])";
			$n+=$N+1;
			if($i==$last-1) { print "\n\n"; } else { print "+ \\\n       "; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_n=0;\\\n";
		print "      (k)=0.;\\\n";
		print "      for (_i=0; _i<$avd; ){\\\n";
		print "         (k)+=";
		my $n=2*$unroll;
		for(my $i=0;$i<2*$unroll;$i++){
			if ($i==0) { print "((r).$cname\[_i\]*(r).$cname\[_i\])"; }
			else { print "((r).$cname\[_i+$i\]*(r).$cname\[_i+$i\])"; }
			if($i==2*$unroll-1) { print ";\\\n"; } else { print "+ \\\n              "; }
		}
		print "         _i+=$n;\\\n";
		print "      }\\\n";
		print "      (k)+=" unless ($avr==0);
		for(my $i=0;$i<$avr;$i++){
			if ($i==0) { print "((r).$cname\[_i\]*(r).$cname\[_i\])"; }
			else { print "((r).$cname\[_i+$i\]*(r).$cname\[_i+$i\])"; }
			if($i==$avr-1) { print ";\\\n"; } else { print "+ \\\n           "; }
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_lc {
  print "/* r=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */\n";
  print "#define _vector_lc_${suff}(r,k1,s1,k2,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_rlc((r).$cname\[$i\],(k1),(s1).$cname\[$i\],(k2),(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_rlc((r).$cname\[_i\],(k1),(s1).$cname\[_i\],(k2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_rlc((r).$cname\[_i\],(k1),(s1).$cname\[_i\],(k2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_lc_add_assign {
  print "/* r+=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */\n";
  print "#define _vector_lc_add_assign_${suff}(r,k1,s1,k2,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_rlc_assign((r).$cname\[$i\],(k1),(s1).$cname\[$i\],(k2),(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_rlc_assign((r).$cname\[_i\],(k1),(s1).$cname\[_i\],(k2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_rlc_assign((r).$cname\[_i\],(k1),(s1).$cname\[_i\],(k2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_clc {
  print "/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */\n";
  print "#define _vector_clc_${suff}(r,z1,s1,z2,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_clc((r).$cname\[$i\],(z1),(s1).$cname\[$i\],(z2),(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_clc((r).$cname\[_i\],(z1),(s1).$cname\[_i\],(z2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_clc((r).$cname\[_i\],(z1),(s1).$cname\[_i\],(z2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_clc_add_assign {
  print "/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */\n";
  print "#define _vector_clc_add_assign_${suff}(r,z1,s1,z2,s2) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_clc_assign((r).$cname\[$i\],(z1),(s1).$cname\[$i\],(z2),(s2).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_clc_assign((r).$cname\[_i\],(z1),(s1).$cname\[_i\],(z2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_clc_assign((r).$cname\[_i\],(z1),(s1).$cname\[_i\],(z2),(s2).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_prod_assign {
  print "/* z+=r^*s (c complex) */\n";
  print "#define _vector_prod_assign_${suff}(z,r,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_prod_assign((z),(r).$cname\[$i\],(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_prod_assign((z),(r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_prod_assign((z),(r).$cname\[_i\],(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_vector_project {
  print "/* r-=z*s (z complex) */\n";
  print "#define _vector_project_${suff}(r,z,s) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		for(my $i=0;$i<$N;$i++){
			print "   _complex_mul_sub_assign((r).$cname\[$i\],(z),(s).$cname\[$i\])";
			if($i==$N-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$vd; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mul_sub_assign((r).$cname\[_i\],(z),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$vr;$i++){
			print "      _complex_mul_sub_assign((r).$cname\[_i\],(z),(s).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

#
# MATRIX VECTOR OPERATIONS
#

sub write_suN_multiply {
  print "/* SU(N) matrix u times SU(N) vector s */\n";
  print "/* r=u*s */\n";
  print "#define _${dataname}_multiply(r,u,s) \\\n";
	if ($N<$Nmax) { #unroll all 
		my ($k)=(0);
		for(my $i=0;$i<$N;$i++){
			print "      _complex_mul((r).$cname\[$i\],(u).$cname\[$k\],(s).$cname\[0\]);\\\n";
			$k++;
			for(my $j=1;$j<$N;$j++){
				print "      _complex_mul_assign((r).$cname\[$i\],(u).$cname\[$k\],(s).$cname\[$j\])";
				$k++;
				if($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_k=0;for (_i=0; _i<$N; ++_i){\\\n";
		print "         _complex_mul((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[0\]); ++_k;\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "         _complex_mul_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[$j\]); ++_k;\\\n";
			}
		} else {
			print "         int _j=1; for (; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            _complex_mul_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[_j\]); ++_k; ++_j;\\\n";
			}
			print "         } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "         _complex_mul_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[_j\]); ++_k; ++_j;\\\n";
			}
		}
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suNr_multiply {
  print "/* SU(N) matrix u times SU(N) vector s */\n";
  print "/* r=u*s */\n";
  print "#define _${rdataname}_multiply(r,u,s) \\\n";
	if ($N<$Nmax) { #unroll all 
		my ($k)=(0);
		for(my $i=0;$i<$N;$i++){
			print "      _complex_mulr((r).$cname\[$i\],(u).$cname\[$k\],(s).$cname\[0\]);\\\n";
			$k++;
			for(my $j=1;$j<$N;$j++){
				print "      _complex_mulr_assign((r).$cname\[$i\],(u).$cname\[$k\],(s).$cname\[$j\])";
				$k++;
				if($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_k=0;for (_i=0; _i<$N; ++_i){\\\n";
		print "         _complex_mulr((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[0\]); ++_k;\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "         _complex_mulr_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[$j\]); ++_k;\\\n";
			}
		} else {
			print "         int _j=1; for (; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            _complex_mulr_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[_j\]); ++_k; ++_j;\\\n";
			}
			print "         } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "         _complex_mulr_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[_j\]); ++_k; ++_j;\\\n";
			}
		}
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_su2_decode {
  my ($su2repr)=@_; #this is unused at the moment, since we only have 1 q. repr

  print "/* SU(2) quaternion u => matrix decode */\n";
  print "/* this declares and uses the temporary array _tmp  */\n";

if ($N==2) {
  print "/* _tmp is a complex 2x2 matrix  */\n";
  print "#define _${dataname}_fun_decode(_tmp,u) \\\n";
  print "      _COMPLEX _tmp[4]; \\\n";
  print "      _tmp[0].re = (u).$cname\[0\]; \\\n";
  print "      _tmp[0].im = (u).$cname\[3\]; \\\n";
  print "      _tmp[1].re = -(u).$cname\[2\]; \\\n";
  print "      _tmp[1].im = (u).$cname\[1\]; \\\n";
  print "      _tmp[2].re = (u).$cname\[2\]; \\\n";
  print "      _tmp[2].im = (u).$cname\[1\]; \\\n";
  print "      _tmp[3].re = (u).$cname\[0\]; \\\n";
  print "      _tmp[3].im = -(u).$cname\[3\]";
  print "\n\n";
} elsif ($N==3) {
  print "/* _tmp is a real 3x3 matrix  */\n";
  print "#define _${rdataname}_adj_decode(_tmp,u) \\\n";
  print "      _REAL _tmp[10];  \\\n";
  print "      _tmp[4]=-2.*(u).$cname\[3\]*(u).$cname\[3\];  \\\n";
  print "      _tmp[8]=-2.*(u).$cname\[1\]*(u).$cname\[1\];  \\\n";
  print "      _tmp[0]=1.+_tmp[8]+_tmp[4];  \\\n";
  print "      _tmp[1]=2.*(u).$cname\[2\]*(u).$cname\[2\];  \\\n";
  print "      _tmp[4]+=1.-_tmp[1];  \\\n";
  print "      _tmp[8]+=1.-_tmp[1];  \\\n";
  print "      _tmp[3]=-2.*(u).$cname\[0\]*(u).$cname\[3\];  \\\n";
  print "      _tmp[7]=2.*(u).$cname\[1\]*(u).$cname\[2\];  \\\n";
  print "      _tmp[1]=_tmp[7]-_tmp[3];  \\\n";
  print "      _tmp[3]+=_tmp[7];  \\\n";
  print "      _tmp[6]=2.*(u).$cname\[0\]*(u).$cname\[1\];  \\\n";
  print "      _tmp[7]=2.*(u).$cname\[2\]*(u).$cname\[3\];  \\\n";
  print "      _tmp[2]=_tmp[7]-_tmp[6];  \\\n";
  print "      _tmp[6]+=_tmp[7];  \\\n";
  print "      _tmp[9]=-2.*(u).$cname\[0\]*(u).$cname\[2\];  \\\n";
  print "      _tmp[7]=2.*(u).$cname\[1\]*(u).$cname\[3\];  \\\n";
  print "      _tmp[5]=_tmp[7]-_tmp[9];  \\\n";
  print "      _tmp[7]+=_tmp[9]";
  print "\n\n";
}

}

sub old_write_su2_multiply {
  print "/* SU(2) matrix u times SU(2) vector s */\n";
  print "/* r=u*s */\n";
  print "/* using quaternionic representations for SU(2) matrices */\n";

  my $localdn;
  my $real="";
  my $localrn;

if ($N==2) { #fundamental representation
  $localrn="fun";
  $localdn=$dataname;
  $real="";
} elsif ($N==3) { #adjoint representation
  $localrn="adj";
  $localdn=$rdataname;
  $real="r";
} else {
  die("Undefined fermion representation in quaternionic code. Exiting...\n");
}

  print "#define _${localdn}_multiply(r,u,s) \\\n";
  print "   do { \\\n";
  print "      _${localdn}_${localrn}_decode(_tmp,(u));  \\\n";
  my ($k)=(0);
    for(my $i=0;$i<$N;$i++){
      print "      _complex_mul${real}((r).$cname\[$i\],_tmp\[$k\],(s).$cname\[0\]);\\\n";
      $k++;
      for(my $j=1;$j<$N;$j++){
        print "      _complex_mul${real}_assign((r).$cname\[$i\],_tmp\[$k\],(s).$cname\[$j\]);\\\n";
	$k++;
    }
  }
  print "   } while(0) \n\n";
}

sub old_write_su2_inverse_multiply {
  print "/* SU(2) matrix u^dagger times SU(2) vector s */\n";
  print "/* r=u^(dagger)*s */\n";
  print "/* using quaternionic representations for SU(2) matrices */\n";

  my $localdn;
  my $real="";
  my $localrn;	
  my $shift=$N*$N-$N-1;

if ($N==2) { #fundamental representation
  $localrn="fun";
  $localdn=$dataname;
  $real="_star";

  print "#define _${localdn}_inverse_multiply(r,u,s) \\\n";
  print "   do { \\\n";
  print "      _${localdn}_${localrn}_decode(_tmp,(u));  \\\n";
  my ($k)=(0);
  for(my $i=0;$i<$N;$i++){
    print "      _complex_mul${real}((r).$cname\[$i\],(s).$cname\[0\],_tmp\[$k\]);\\\n";
    for(my $j=1;$j<$N;$j++){
      $k+=$N;
      print "      _complex_mul${real}_assign((r).$cname\[$i\],(s).$cname\[$j\],_tmp\[$k\]);\\\n";
    }
    $k-=$shift;
  }

} elsif ($N==3) { #adjoint representation
  $localrn="adj";
  $localdn=$rdataname;
  $real="r";

  print "#define _${localdn}_inverse_multiply(r,u,s) \\\n";
  print "   do { \\\n";
  print "      _${localdn}_${localrn}_decode(_tmp,(u));  \\\n";
  my ($k)=(0);
  for(my $i=0;$i<$N;$i++){
    print "      _complex_mul${real}((r).$cname\[$i\],_tmp\[$k\],(s).$cname\[0\]);\\\n";
    for(my $j=1;$j<$N;$j++){
      $k+=$N;
      print "      _complex_mul${real}_assign((r).$cname\[$i\],_tmp\[$k\],(s).$cname\[$j\]);\\\n";
    }
    $k-=$shift;
  }

} else {
  die("Undefined fermion representation in quaternionic code. Exiting...\n");
}

  print "   } while(0) \n\n";

}



sub write_suN_inverse_multiply {
	print "/* SU(N) matrix u^dagger times SU(N) vector s */\n";
	print "/* r=(u^dagger)*s */\n";
	print "#define _${dataname}_inverse_multiply(r,u,s) \\\n";
	my $shift=$N*$N-$N-1;
	if ($N<$Nmax) { #unroll all 
		my ($k)=(0);
		for(my $i=0;$i<$N;$i++){
			print "      _complex_mul_star((r).$cname\[$i\],(s).$cname\[0\],(u).$cname\[$k\]);\\\n";
			for(my $j=1;$j<$N;$j++){
				$k+=$N;
				print "      _complex_mul_star_assign((r).$cname\[$i\],(s).$cname\[$j\],(u).$cname\[$k\])";
				if($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
			$k-=$shift;
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_k=0;for (_i=0; _i<$N; ++_i){\\\n";
		print "         _complex_mul_star((r).$cname\[_i\],(s).$cname\[0\],(u).$cname\[_k\]);\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "         _k+=$N; _complex_mul_star_assign((r).$cname\[_i\],(s).$cname\[$j\],(u).$cname\[_k\]);\\\n";
			}
		} else {
			print "         int _j=1; for (; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            _k+=$N; _complex_mul_star_assign((r).$cname\[_i\],(s).$cname\[_j\],(u).$cname\[_k\]); ++_j;\\\n";
			}
			print "         } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "         _k+=$N; _complex_mul_star_assign((r).$cname\[_i\],(s).$cname\[_j\],(u).$cname\[_k\]); ++_j;\\\n";
			}
		}
		print "         _k-=$shift;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suNr_inverse_multiply {
  print "/* SU(N) matrix u^dagger times SU(N) vector s */\n";
  print "/* r=(u^dagger)*s */\n";
  print "#define _${rdataname}_inverse_multiply(r,u,s) \\\n";
	my $shift=$N*$N-$N-1;
	if ($N<$Nmax) { #unroll all 
		my ($k)=(0);
		for(my $i=0;$i<$N;$i++){
			print "      _complex_mulr((r).$cname\[$i\],(u).$cname\[$k\],(s).$cname\[0\]);\\\n";
			for(my $j=1;$j<$N;$j++){
				$k+=$N;
				print "      _complex_mulr_assign((r).$cname\[$i\],(u).$cname\[$k\],(s).$cname\[$j\])";
				if($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
			$k-=$shift;
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_k=0;for (_i=0; _i<$N; ++_i){\\\n";
		print "         _complex_mulr((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[0\]);\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "         _k+=$N; _complex_mulr_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[$j\]);\\\n";
			}
		} else {
			print "         int _j=1; for (; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            _k+=$N; _complex_mulr_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[_j\]); ++_j;\\\n";
			}
			print "         } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "         _k+=$N; _complex_mulr_assign((r).$cname\[_i\],(u).$cname\[_k\],(s).$cname\[_j\]); ++_j;\\\n";
			}
		}
		print "         _k-=$shift;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

#
# MATRIX-MATRIX OPERATIONS
#

sub write_suN_dagger {
  print "/* u=v^dagger */\n";
  print "#define _${dataname}_dagger(u,v) \\\n";
	my $shift=$N*$N-$N-1;
	if ($N<$Nmax) { #unroll all 
		my ($n,$k)=(0,0);
		for(my $i=1;$i<=$N;$i++){
			for(my $j=1;$j<=$N;$j++){
				print "   _complex_star((u).$cname\[$n\],(v).$cname\[$k\])";
				if ($j!=$N) {$n++; $k+=$N;}
				if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
			}
			$n++; $k-=$shift;
		}
	} else { #partial unroll
		print "   do { \\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_n=0,_k=0;\\\n";
		} else {
		    print "      int _i,_j,_n=0,_k=0;\\\n";
		}
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         _complex_star((u).$cname\[_n\],(v).$cname\[_k\]);\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "         ++_n; _k+=$N; _complex_star((u).$cname\[_n\],(v).$cname\[_k\]);\\\n";
			}
		} else {
			print "         for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            ++_n; _k+=$N; _complex_star((u).$cname\[_n\],(v).$cname\[_k\]); ++_j;\\\n";
			}
			print "         } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "         ++_n; _k+=$N; _complex_star((u).$cname\[_n\],(v).$cname\[_k\]);\\\n";
			}
		}
		print "         ++_n; _k-=$shift;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}


sub write_suNr_dagger {
    print "/* u=v^dagger */\n";
    print "#define _${rdataname}_dagger(u,v) \\\n";
	my $shift=$N*$N-$N-1;
	if ($N<$Nmax) { #unroll all 
		my ($n,$k)=(0,0);
		for(my $i=1;$i<=$N;$i++){
			for(my $j=1;$j<=$N;$j++){
				print "   (u).$cname\[$n\]=(v).$cname\[$k\]";
				if ($j!=$N) {$n++; $k+=$N;}
				if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
			}
			$n++; $k-=$shift;
		}
	} else { #partial unroll
		print "   {\\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_n=0,_k=0;\\\n";
		} else {
		    print "      int _i,_j,_n=0,_k=0;\\\n";
		}
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         (u).$cname\[_n\]=(v).$cname\[_k\];\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "         ++_n; _k+=$N; (u).$cname\[_n\]=(v).$cname\[_k\];\\\n";
			}
		} else {
			print "         for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            ++_n; _k+=$N; (u).$cname\[_n\]=(v).$cname\[_k\]; ++_j;\\\n";
			}
			print "         } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "         ++_n; _k+=$N; (u).$cname\[_n\]=(v).$cname\[_k\];\\\n";
			}
		}
		print "         ++_n; _k-=$shift;\\\n";
		print "      }\\\n";
		print "   }((void)0) \n\n";
	}
}

sub write_suN_times_suN {
  print "/* u=v*w */\n";
  print "#define _${dataname}_times_${dataname}(u,v,w) \\\n";
	my $shift=$N*$N-$N-1;
	my $shift2=$N-1;
	if ($N<$Nmax) { #unroll all 
	my ($n,$k,$l)=(0,0,0);
	for(my $i=0;$i<$N;$i++){
		for(my $y=0;$y<$N;$y++){
			print "      _complex_mul((u).$cname\[$n\],(v).$cname\[$k\],(w).$cname\[$l\]);\\\n";
			for(my $j=1;$j<$N;$j++){
				$k++; $l+=$N;
				print "      _complex_mul_assign((u).$cname\[$n\],(v).$cname\[$k\],(w).$cname\[$l\])";
				if($i==$N-1 and $y==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
			$n++;
			$k-=$shift2;
			$l-=$shift;
		}
		$k+=$N;
		$l=0;
	}
	} else { #partial unroll
		print "   do { \\\n";

		if($N<(2*$unroll+1)) {	
		    print "      int _i,_y,_n=0,_k=0,_l=0;\\\n"; 
		} else {
		    print "      int _i,_y,_j,_n=0,_k=0,_l=0;\\\n"; 
		}

		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         for (_y=0; _y<$N; ++_y){\\\n";
		print "            _complex_mul((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]);\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "            ++_k; _l+=$N; _complex_mul_assign((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]);\\\n";
			}
		} else {
			print "            for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "               ++_k; _l+=$N; _complex_mul_assign((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]); ++_j;\\\n";
			}
			print "            } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "            ++_k; _l+=$N; _complex_mul_assign((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]);\\\n";
			}
		}
		print "            ++_n; _k-=$shift2; _l-=$shift;\\\n";
		print "         } _k+=$N; _l=0;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suNr_times_suNr {
  print "/* u=v*w */\n";
  print "#define _${rdataname}_times_${rdataname}(u,v,w) \\\n";
	my $shift=$N*$N-$N-1;
	my $shift2=$N-1;
	if ($N<$Nmax) { #unroll all 
	my ($n,$k,$l)=(0,0,0);
	for(my $i=0;$i<$N;$i++){
		for(my $y=0;$y<$N;$y++){
			print "      (u).$cname\[$n\]=(v).$cname\[$k\]*(w).$cname\[$l\];\\\n";
			for(my $j=1;$j<$N;$j++){
				$k++; $l+=$N;
				print "      (u).$cname\[$n\]+=(v).$cname\[$k\]*(w).$cname\[$l\]";
				if($i==$N-1 and $y==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
			$n++;
			$k-=$shift2;
			$l-=$shift;
		}
		$k+=$N;
		$l=0;
	}
	} else { #partial unroll
		print "   do { \\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_y,_n=0,_k=0,_l=0;\\\n";
		} else {
		    print "      int _i,_y,_j,_n=0,_k=0,_l=0;\\\n";
		}		    
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         for (_y=0; _y<$N; ++_y){\\\n";
		print "            (u).$cname\[_n\]=(v).$cname\[_k\]*(w).$cname\[_l\];\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "            ++_k; _l+=$N; (u).$cname\[_n\]+=(v).$cname\[_k\]*(w).$cname\[_l\];\\\n";
			}
		} else {
			print "            for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "               ++_k; _l+=$N; (u).$cname\[_n\]+=(v).$cname\[_k\]*(w).$cname\[_l\]; ++_j;\\\n";
			}
			print "            } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "            ++_k; _l+=$N; (u).$cname\[_n\]+=(v).$cname\[_k\]*(w).$cname\[_l\];\\\n";
			}
		}
		print "            ++_n; _k-=$shift2; _l-=$shift;\\\n";
		print "         } _k+=$N; _l=0;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suN_times_suN_dagger {
  print "/* u=v*w^+ */\n";
  print "#define _${dataname}_times_${dataname}_dagger(u,v,w) \\\n";
#	my $shift=$N*$N-$N-1;
	my $shift2=$N-1;
	if ($N<$Nmax) { #unroll all 
	my ($n,$k,$l)=(0,0,0);
	for(my $i=0;$i<$N;$i++){
		for(my $y=0;$y<$N;$y++){
			print "      _complex_mul_star((u).$cname\[$n\],(v).$cname\[$k\],(w).$cname\[$l\]);\\\n";
			for(my $j=1;$j<$N;$j++){
				$k++; $l++;
				print "      _complex_mul_star_assign((u).$cname\[$n\],(v).$cname\[$k\],(w).$cname\[$l\])";
				if($i==$N-1 and $y==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
			$n++;
			$k-=$shift2;
			$l++;
		}
		$k+=$N;
		$l=0;
	}
	} else { #partial unroll
		print "   do { \\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_y,_n=0,_k=0,_l=0;\\\n";
		} else {
		    print "      int _i,_y,_j,_n=0,_k=0,_l=0;\\\n";
		}		    
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         for (_y=0; _y<$N; ++_y){\\\n";
		print "            _complex_mul_star((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]);\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "            ++_k; ++_l; _complex_mul_star_assign((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]);\\\n";
			}
		} else {
			print "            for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "               ++_k; ++_l; _complex_mul_star_assign((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]); ++_j;\\\n";
			}
			print "            } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "            ++_k; ++_l; _complex_mul_star_assign((u).$cname\[_n\],(v).$cname\[_k\],(w).$cname\[_l\]);\\\n";
			}
		}
		print "            ++_n; _k-=$shift2; ++_l;\\\n";
		print "         } _k+=$N; _l=0;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suN_dagger_times_suN {
  print "/* u=v^+*w */\n";
  print "#define _${dataname}_dagger_times_${dataname}(u,v,w) \\\n";
	my $shift=$N*$N-$N-1;
	my $shift2=$N-1;
	if ($N<$Nmax) { #unroll all 
	my ($n,$k,$l)=(0,0,0);
	my ($v,$w,$u)=(0,0,0);
	for(my $i=0;$i<$N;$i++){
		for(my $y=0;$y<$N;$y++){
		    $u=$i*$N+$y;
		    $v=$i;
		    $w=$y;		    
			print "      _complex_mul_star((u).$cname\[$u\],(w).$cname\[$w\],(v).$cname\[$v\]);\\\n";
			for(my $j=1;$j<$N;$j++){
			    $v=$j*$N+$i;
			    $w=$j*$N+$y;
			    print "      _complex_mul_star_assign((u).$cname\[$u\],(w).$cname\[$w\],(v).$cname\[$v\])";
			    if($i==$N-1 and $y==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
		}
	}
	} else { #partial unroll
		print "   do { \\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_y,_n=0,_k=0,_l=0;\\\n";
		} else {
		    print "      int _i,_y,_j,_n=0,_k=0,_l=0;\\\n";
		}		    
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         for (_y=0; _y<$N; ++_y){\\\n";
 		print "            _k=_y; _l=_i;\\\n";
		print "            _complex_mul_star((u).$cname\[_n\],(w).$cname\[_k\],(v).$cname\[_l\]);\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "            _k+=$N; _l+=$N; _complex_mul_star_assign((u).$cname\[_n\],(w).$cname\[_k\],(v).$cname\[_l\]);\\\n";
			}
		} else {
			print "            for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "               _k+=$N; _l+=$N; _complex_mul_star_assign((u).$cname\[_n\],(w).$cname\[_k\],(v).$cname\[_l\]); ++_j;\\\n";
			}
			print "            } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "            _k+=$N; _l+=$N; _complex_mul_star_assign((u).$cname\[_n\],(w).$cname\[_k\],(v).$cname\[_l\]);\\\n";
			}
		}
		print "            ++_n;\\\n";
		print "         }\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suNr_times_suNr_dagger {
  print "/* u=v*w^t */\n";
  print "#define _${rdataname}_times_${rdataname}_dagger(u,v,w) \\\n";
	##my $shift=$N*$N-$N-1;
	my $shift2=$N-1;
	if ($N<$Nmax) { #unroll all 
	my ($n,$k,$l)=(0,0,0);
	for(my $i=0;$i<$N;$i++){
		for(my $y=0;$y<$N;$y++){
			print "      (u).$cname\[$n\]=(v).$cname\[$k\]*(w).$cname\[$l\];\\\n";
			for(my $j=1;$j<$N;$j++){
				$k++; $l++;
				print "      (u).$cname\[$n\]+=(v).$cname\[$k\]*(w).$cname\[$l\]";
				if($i==$N-1 and $y==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
			$n++;
			$k-=$shift2;
			$l++;
		}
		$k+=$N;
		$l=0;
	}
	} else { #partial unroll
		print "   do { \\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_y,_n=0,_k=0,_l=0;\\\n";
		} else { 
		    print "      int _i,_y,_j,_n=0,_k=0,_l=0;\\\n";
		}
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         for (_y=0; _y<$N; ++_y){\\\n";
		print "            (u).$cname\[_n\]=(v).$cname\[_k\]*(w).$cname\[_l\];\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "            ++_k; ++_l; (u).$cname\[_n\]+=(v).$cname\[_k\]*(w).$cname\[_l\];\\\n";
			}
		} else {
			print "            for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "               ++_k; ++_l; (u).$cname\[_n\]+=(v).$cname\[_k\]*(w).$cname\[_l\]; ++_j;\\\n";
			}
			print "            } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "            ++_k; ++_l; (u).$cname\[_n\]+=(v).$cname\[_k\]*(w).$cname\[_l\];\\\n";
			}
		}
		print "            ++_n; _k-=$shift2; ++_l;\\\n";
		print "         } _k+=$N; _l=0;\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suNr_dagger_times_suNr {
  print "/* u=v^+*w */\n";
  print "#define _${rdataname}_dagger_times_${rdataname}(u,v,w) \\\n";
	my $shift=$N*$N-$N-1;
	my $shift2=$N-1;
	if ($N<$Nmax) { #unroll all 
	my ($n,$k,$l)=(0,0,0);
	my ($v,$w,$u)=(0,0,0);
	for(my $i=0;$i<$N;$i++){
		for(my $y=0;$y<$N;$y++){
		    $u=$i*$N+$y;
		    $v=$i;
		    $w=$y;		    
			print "     (u).$cname\[$u\]=(w).$cname\[$w\]*(v).$cname\[$v\];\\\n";
			for(my $j=1;$j<$N;$j++){
			    $v=$j*$N+$i;
			    $w=$j*$N+$y;
			    print "      (u).$cname\[$u\]+=(w).$cname\[$w\]*(v).$cname\[$v\]";
			    if($i==$N-1 and $y==$N-1 and $j==$N-1) { print "\n\n"; } else { print "; \\\n"; }
			}
		}
	}
	} else { #partial unroll
		print "   do { \\\n";
		if($N<(2*$unroll+1)) {
		    print "      int _i,_y,_n=0,_k=0,_l=0;\\\n";
		} else {
		    print "      int _i,_y,_j,_n=0,_k=0,_l=0;\\\n";
		}		    
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		print "         for (_y=0; _y<$N; ++_y){\\\n";
 		print "            _k=_y; _l=_i;\\\n";
		print "            (u).$cname\[_n\]=(w).$cname\[_k\]*(v).$cname\[_l\];\\\n";
		if($N<(2*$unroll+1)) {
			for(my $j=1;$j<$N;$j++){
				print "            _k+=$N; _l+=$N; (u).$cname\[_n\]+=(w).$cname\[_k\]*(v).$cname\[_l\];\\\n";
			}
		} else {
			print "            for (_j=0; _j<$md; ){ \\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "               _k+=$N; _l+=$N; (u).$cname\[_n\]+=(w).$cname\[_k\]*(v).$cname\[_l\]; ++_j;\\\n";
			}
			print "            } \\\n";
			for(my $i=0;$i<$mr;$i++){
				print "            _k+=$N; _l+=$N; (u).$cname\[_n\]+=(w).$cname\[_k\]*(v).$cname\[_l\];\\\n";
			}
		}
		print "            ++_n;\\\n";
		print "         }\\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}


sub write_suN_zero {
  print "/* u=0 */\n";
  print "#define _${dataname}_zero(u) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "    _complex_0((u).$cname\[$i\])";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_0((u).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      _complex_0((u).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_zero {
  print "/* u=0 */\n";
  print "#define _${rdataname}_zero(u) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "    (u).$cname\[$i\]=0.";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (u).$cname\[_i\]=0.; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (u).$cname\[_i\]=0.; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_unit {
  print "/* u=1 */\n";
  print "#define _${dataname}_unit(u) \\\n";
	my $shift=$N+1;
	if ($N<$Nmax) { #unroll all 
		my $n=0;
		for (my $i=0; $i<$N; $i++) {
			for (my $j=0; $j<$N; $j++) {
				if ($i==$j) { 
					print "   _complex_1((u).$cname\[$n\])";
				} else {
					print "   _complex_0((u).$cname\[$n\])";
				}
				if ($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print ";\\\n"; }
				$n++;
			}
		}
	} else {
		print "   do { \\\n";
		print "      _${dataname}_zero((u));\\\n";
		if ($N<(2*$unroll+1)) {
			my $n=0;
			for (my $i=0; $i<$N; $i++) {
				print "      _complex_1((u).$cname\[$n\]);\\\n";
				$n+=$shift;
			}
		} else {
			print "      int _i,_n=0; for (_i=0; _i<$vd; ){\\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "         _complex_1((u).$cname\[_n\]); _n+=$shift; ++_i;\\\n";
			}
			print "      }\\\n";
			for(my $i=0;$i<$vr;$i++){
				print "      _complex_1((u).$cname\[_n\]); _n+=$shift;\\\n";
			}
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_unit {
  print "/* u=1 */\n";
  print "#define _${rdataname}_unit(u) \\\n";
	my $dim=$N*$N;
	my $shift=$N+1;
	if ($N<$Nmax) { #unroll all 
		my $n=0;
		for (my $i=0; $i<$N; $i++) {
			for (my $j=0; $j<$N; $j++) {
				if ($i==$j) { 
					print "   (u).$cname\[$n\]=1.";
				} else {
					print "   (u).$cname\[$n\]=0.";
				}
				if ($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print ";\\\n"; }
				$n++;
			}
		}
	} else {
		print "   do { \\\n";
		print "      _${rdataname}_zero((u));\\\n";
		if ($N<(2*$unroll+1)) {
			my $n=0;
			for (my $i=0; $i<$N; $i++) {
				print "      (u).$cname\[$n\]=1.;\\\n";
				$n+=$shift;
			}
		} else {
			print "      int _i,_n=0; for (_i=0; _i<$vd; ){\\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "         (u).$cname\[_n\]=1.; _n+=$shift; ++_i;\\\n";
			}
			print "      }\\\n";
			for(my $i=0;$i<$vr;$i++){
				print "      (u).$cname\[_n\]=1.; _n+=$shift;\\\n";
			}
		}
		print "   } while(0) \n\n";}
}

sub write_suN_copy {
  print "/* u=v */\n";
  print "#define _${dataname}_copy(u,v) \\\n";
  print "   (u)=(v)\n\n";
}

sub write_suN_minus {
  print "/* u=-v */\n";
  print "#define _${dataname}_minus(u,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   _complex_minus((u).$cname\[$i\],(v).$cname\[$i\])";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_minus((u).$cname\[_i\],(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      _complex_minus((u).$cname\[_i\],(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_minus {
  print "/* u=-v */\n";
  print "#define _${rdataname}_minus(u,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   (u).$cname\[$i\]=-(v).$cname\[$i\]";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (u).$cname\[_i\]=-(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (u).$cname\[_i\]=-(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_mul {
  print "/* u=r*v (u,v mat, r real) */\n";
  print "#define _${dataname}_mul(u,r,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   _complex_mulr((u).$cname\[$i\],(r),(v).$cname\[$i\])";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mulr((u).$cname\[_i\],(r),(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      _complex_mulr((u).$cname\[_i\],(r),(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_mulc {
  print "/* u=r*v (u,v mat, r complex) */\n";
  print "#define _${dataname}_mulc(u,r,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   _complex_mul((u).$cname\[$i\],(r),(v).$cname\[$i\])";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_mul((u).$cname\[_i\],(r),(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      _complex_mul((u).$cname\[_i\],(r),(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}


sub write_suNr_mul {
  print "/* u=r*v (u,v mat, r real) */\n";
  print "#define _${rdataname}_mul(u,r,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   (u).$cname\[$i\]=(r)*(v).$cname\[$i\]";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (u).$cname\[_i\]=(r)*(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (u).$cname\[_i\]=(r)*(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_add_assign {
  print "/* u+=v */\n";
  print "#define _${dataname}_add_assign(u,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   _complex_add_assign((u).$cname\[$i\],(v).$cname\[$i\])";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_add_assign((u).$cname\[_i\],(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      _complex_add_assign((u).$cname\[_i\],(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_add_assign {
  print "/* u+=v */\n";
  print "#define _${rdataname}_add_assign(u,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   (u).$cname\[$i\]+=(v).$cname\[$i\]";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (u).$cname\[_i\]+=(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (u).$cname\[_i\]+=(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_sub_assign {
  print "/* u-=v */\n";
  print "#define _${dataname}_sub_assign(u,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   _complex_sub_assign((u).$cname\[$i\],(v).$cname\[$i\])";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         _complex_sub_assign((u).$cname\[_i\],(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      _complex_sub_assign((u).$cname\[_i\],(v).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_sub_assign {
  print "/* u-=v */\n";
  print "#define _${rdataname}_sub_assign(u,v) \\\n";
	my $dim=$N*$N;
	if ($N<$Nmax or $dim<(2*$unroll+1)) { #unroll all 
		for(my $i=0; $i<$dim; $i++) {
			print "   (u).$cname\[$i\]-=(v).$cname\[$i\]";
			if($i==$dim-1) {print "\n\n";} else { print ";\\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (u).$cname\[_i\]-=(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (u).$cname\[_i\]-=(v).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_sqnorm {
  print "/* k=| u |2 ) */\n";
	print "#define _${dataname}_sqnorm(k,u) \\\n";
	my $dim=$N*$N;
  if ($N<$Nmax or $dim<(2*$unroll+1) ) { #unroll all 
		print "   (k)=0.;\\\n";
		for(my $i=0;$i<$dim;$i++){
			print "   (k)+=_complex_prod_re((u).$cname\[$i\],(u).$cname\[$i\])";
			if($i==$dim-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)=0.;\\\n";
		print "      for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)+=_complex_prod_re((u).$cname\[_i\],(u).$cname\[_i\]); ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (k)+=_complex_prod_re((u).$cname\[_i\],(u).$cname\[_i\]); ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_sqnorm {
  print "/* k= | u |2 ) */\n";
	print "#define _${rdataname}_sqnorm(k,u) \\\n";
	my $dim=$N*$N;
  if ($N<$Nmax or $dim<(2*$unroll+1) ) { #unroll all 
		print "   (k)=0.;\\\n";
		for(my $i=0;$i<$dim;$i++){
			print "   (k)+=(u).$cname\[$i\]*(u).$cname\[$i\]";
			if($i==$dim-1) { print "\n\n"; } else { print "; \\\n"; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i;\\\n";
		print "      (k)=0.;\\\n";
		print "      for (_i=0; _i<$md2; ){\\\n";
		for(my $i=0;$i<$unroll;$i++){
			print "         (k)+=(u).$cname\[_i\]*(u).$cname\[_i\]; ++_i;\\\n";
		}
		print "      }\\\n";
		for(my $i=0;$i<$mr2;$i++){
			print "      (k)+=(u).$cname\[_i\]*(u).$cname\[_i\]; ++_i;\\\n";
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_sqnorm_m1 {
  print "/* k=| 1 - u |2 ) */\n";
  print "#define _${dataname}_sqnorm_m1(k,u) \\\n";
	my $shift=$N*$N-1;
	if ($N<2*$Nmax) { #unroll all : here we use an higher Nmax because we cannot unroll this
	    print "   (k)=\\\n    ";
	    my $n=0;
		for(my $i=1;$i<=$N;$i++){
			for(my $j=1;$j<=$N;$j++){
				if ($i==$j) {
					print "+_complex_prod_m1_re((u).$cname\[$n\],(u).$cname\[$n\])";
				} else {
					print "+_complex_prod_re((u).$cname\[$n\],(u).$cname\[$n\])";
				}
				$n++;
				if($j==$N and $i==$N) {print "\n\n";} else {print "\\\n    ";}
			}
		}
	} else {
		print "   do { \\\n";
		print "      (k)=0.;\\\n";
		print "      int _i,_j,_n=0,_l=0,_s2=0;\\\n";
		print "      for(_i=0;_i<$N;){\\\n";
		print "         (k)+=_complex_prod_m1_re((u).$cname\[_n\],(u).$cname\[_n\]);\\\n";
		print "         ++_n; _l+=$N; \\\n";
		print "         for(_j=_i+1;_j<$N;++_j){\\\n";
		print "            (k)+=_complex_prod_re((u).$cname\[_n\],(u).$cname\[_n\]);\\\n";
		print "            (k)+=_complex_prod_re((u).$cname\[_l\],(u).$cname\[_l\]);\\\n";
		print "            ++_n; _l+=$N; \\\n";
		print "         }\\\n";
		print "         ++_i; _s2+=$N; _n+=_i; _l-=$shift-_s2; \\\n";
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suN_trace_re {
  print "/* k=Re Tr (u) */\n";
  print "#define _${dataname}_trace_re(k,u) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)=";
		my $n=0;
		for(my $i=0;$i<$N;$i++){
			print "_complex_re((u).$cname\[$n\])";
			$n+=$N+1;
			if($i==$N-1) { print "\n\n"; } else { print "+ \\\n       "; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_n=0;\\\n";
		print "      (k)=0.;\\\n";
		print "      for (_i=0; _i<$vd; _i+=$unroll){\\\n";
		print "         (k)+=";
		my $n=0;
		for(my $i=0;$i<$unroll;$i++){
			if ($i==0) { print "_complex_re((u).$cname\[_n\])"; }
			else { print "_complex_re((u).$cname\[_n+$n\])"; }
			$n+=$N+1;
			if($i==$unroll-1) { print ";\\\n"; } else { print "+ \\\n              "; }
		}
		print "         _n+=$n;\\\n";
		print "      }\\\n";
		$n=0;
		print "      (k)+=" unless ($vr==0);
		for(my $i=0;$i<$vr;$i++){
			if ($i==0) { print "_complex_re((u).$cname\[_n\])"; }
			else { print "_complex_re((u).$cname\[_n+$n\])"; }
			$n+=$N+1;
			if($i==$vr-1) { print ";\\\n"; } else { print "+ \\\n           "; }
		}
		print "   } while(0) \n\n";
	}
}

sub write_suNr_trace_re {
  print "/* k=Re Tr (u) */\n";
  print "#define _${rdataname}_trace_re(k,u) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)=";
		my $n=0;
		for(my $i=0;$i<$N;$i++){
			print "(u).$cname\[$n\]";
			$n+=$N+1;
			if($i==$N-1) { print "\n\n"; } else { print "+ \\\n       "; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_n=0;\\\n";
		print "      (k)=0.;\\\n";
		print "      for (_i=0; _i<$vd; _i+=$unroll){\\\n";
		print "         (k)+=";
		my $n=0;
		for(my $i=0;$i<$unroll;$i++){
			if ($i==0) { print "(u).$cname\[_n\]"; }
			else { print "(u).$cname\[_n+$n\]"; }
			$n+=$N+1;
			if($i==$unroll-1) { print ";\\\n"; } else { print "+ \\\n              "; }
		}
		print "         _n+=$n;\\\n";
		print "      }\\\n";
		$n=0;
		print "      (k)+=" unless ($vr==0);
		for(my $i=0;$i<$vr;$i++){
			if ($i==0) { print "(u).$cname\[_n\]"; }
			else { print "(u).$cname\[_n+$n\]"; }
			$n+=$N+1;
			if($i==$vr-1) { print ";\\\n"; } else { print "+ \\\n           "; }
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_trace_im {
  print "/* k=Im Tr (u) */\n";
  print "#define _${dataname}_trace_im(k,u) \\\n";
  if ($N<$Nmax or $N<(2*$unroll+1) ) { #unroll all 
		print "   (k)=";
		my $n=0;
		for(my $i=0;$i<$N;$i++){
			print "_complex_im((u).$cname\[$n\])";
			$n+=$N+1;
			if($i==$N-1) { print "\n\n"; } else { print "+ \\\n       "; }
		}
	} else { #partial unroll
		print "   do { \\\n";
		print "      int _i,_n=0;\\\n";
		print "      (k)=0.;\\\n";
		print "      for (_i=0; _i<$vd; _i+=$unroll){\\\n";
		print "         (k)+=";
		my $n=0;
		for(my $i=0;$i<$unroll;$i++){
			if ($i==0) { print "_complex_im((u).$cname\[_n\])"; }
			else { print "_complex_im((u).$cname\[_n+$n\])"; }
			$n+=$N+1;
			if($i==$unroll-1) { print ";\\\n"; } else { print "+ \\\n              "; }
		}
		print "         _n+=$n;\\\n";
		print "      }\\\n";
		$n=0;
		print "      (k)+=" unless ($vr==0);
		for(my $i=0;$i<$vr;$i++){
			if ($i==0) { print "_complex_im((u).$cname\[_n\])"; }
			else { print "_complex_im((u).$cname\[_n+$n\])"; }
			$n+=$N+1;
			if($i==$vr-1) { print ";\\\n"; } else { print "+ \\\n           "; }
		}
		print "   } while(0) \n\n";
	}
}

sub write_suN_2TA {
  print "/* u=v - v^+ -1/N Tr(v - v^+)*I */\n";
  print "#define _${dataname}_2TA(u,v) \\\n";
  print "   do { \\\n";
	my $shift=$N*$N-1;
	if ($N<2*$Nmax) { #unroll all : here we use an higher Nmax because we cannot unroll this
		print "      double _trim; _${dataname}_trace_im(_trim,(v)); _trim*=(2./$N.);\\\n";
		my ($n,$k)=(0,0);
		my $s2=0;
		for(my $i=0;$i<$N;){
			for(my $j=$i;$j<$N;$j++){
				if($i==$j) {
					print "      (u).$cname\[$n\].re= 0.; \\\n";
					print "      (u).$cname\[$n\].im= 2.*(v).$cname\[$k\].im-_trim; \\\n";
				} else {
					print "      (u).$cname\[$n\].re= (v).$cname\[$n\].re-(v).$cname\[$k\].re; \\\n";
					print "      (u).$cname\[$n\].im= (v).$cname\[$k\].im+(v).$cname\[$n\].im; \\\n";
					print "      (u).$cname\[$k\].re= -(u).$cname\[$n\].re; \\\n";
					print "      (u).$cname\[$k\].im= (u).$cname\[$n\].im; \\\n";
				}
				$n++; $k+=$N;
			}
			$i++;
			$s2+=$N;
			$n+=$i; $k-=$shift-$s2;
		}
	} else {
	  print "      int _i,_j,_n=0,_k=0,_s2=0;\\\n";
		print "      double _trim; _${dataname}_trace_im(_trim,(v)); _trim*=(2./$N.);\\\n";
    print "      for(_i=0;_i<$N;){\\\n";
		print "         (u).$cname\[_n\].re= 0.; \\\n";
		print "         (u).$cname\[_n\].im= 2.*(v).$cname\[_k\].im-_trim; \\\n";
		print "         ++_n; _k+=$N; \\\n";
    print "         for(_j=_i+1;_j<$N;++_j){\\\n";
		print "            (u).$cname\[_n\].re= (v).$cname\[_n\].re-(v).$cname\[_k\].re; \\\n";
		print "            (u).$cname\[_n\].im= (v).$cname\[_k\].im+(v).$cname\[_n\].im; \\\n";
		print "            (u).$cname\[_k\].re= -(u).$cname\[_n\].re; \\\n";
		print "            (u).$cname\[_k\].im= (u).$cname\[_n\].im; \\\n";
		print "            ++_n; _k+=$N; \\\n";
		print "         }\\\n";
		print "         ++_i; _s2+=$N; _n+=_i; _k-=$shift-_s2; \\\n";
		print "      }\\\n";
	}
	print "   } while(0) \n\n";
}

sub write_suN_TA {
  print "/* u=0.5(v - v^+) -1/(2N) Tr(v - v^+)*I */\n";
  print "#define _${dataname}_TA(u,v) \\\n";
  print "   do { \\\n";
	my $shift=$N*$N-1;
	if ($N<2*$Nmax) { #unroll all : here we use an higher Nmax because we cannot unroll this
		print "      double _trim; _${dataname}_trace_im(_trim,(v)); _trim*=(1./$N.);\\\n";
		my ($n,$k)=(0,0);
		my $s2=0;
		for(my $i=0;$i<$N;){
			for(my $j=$i;$j<$N;$j++){
				if($i==$j) {
					print "      (u).$cname\[$n\].re= 0.; \\\n";
					print "      (u).$cname\[$n\].im= (v).$cname\[$k\].im-_trim; \\\n";
				} else {
					print "      (u).$cname\[$n\].re= 0.5*((v).$cname\[$n\].re-(v).$cname\[$k\].re); \\\n";
					print "      (u).$cname\[$n\].im= 0.5*((v).$cname\[$k\].im+(v).$cname\[$n\].im); \\\n";
					print "      (u).$cname\[$k\].re= -(u).$cname\[$n\].re; \\\n";
					print "      (u).$cname\[$k\].im= (u).$cname\[$n\].im; \\\n";
				}
				$n++; $k+=$N;
			}
			$i++;
			$s2+=$N;
			$n+=$i; $k-=$shift-$s2;
		}
	} else {
	  print "      int _i,_j,_n=0,_k=0,_s2=0;\\\n";
		print "      double _trim; _${dataname}_trace_im(_trim,(v)); _trim*=(1./$N.);\\\n";
    print "      for(_i=0;_i<$N;){\\\n";
		print "         (u).$cname\[_n\].re= 0.; \\\n";
		print "         (u).$cname\[_n\].im= (v).$cname\[_k\].im-_trim; \\\n";
		print "         ++_n; _k+=$N; \\\n";
    print "         for(_j=_i+1;_j<$N;++_j){\\\n";
		print "            (u).$cname\[_n\].re= 0.5*((v).$cname\[_n\].re-(v).$cname\[_k\].re); \\\n";
		print "            (u).$cname\[_n\].im= 0.5*((v).$cname\[_k\].im+(v).$cname\[_n\].im); \\\n";
		print "            (u).$cname\[_k\].re= -(u).$cname\[_n\].re; \\\n";
		print "            (u).$cname\[_k\].im= (u).$cname\[_n\].im; \\\n";
		print "            ++_n; _k+=$N; \\\n";
		print "         }\\\n";
		print "         ++_i; _s2+=$N; _n+=_i; _k-=$shift-_s2; \\\n";
		print "      }\\\n";
	}
	print "   } while(0) \n\n";
}

sub write_suN_FMAT {
  print "/* This fuction compute the hmc force matrix */\n";
  print "/* this fuction accumulates on the original matrix u */\n";
  print "#define _${dataname}_FMAT(u,s) \\\n";
	if ($N<$Nmax) { #unroll all 
		my $n=0;
		for(my $i=0;$i<$N;$i++){
			for(my $j=0;$j<$N;$j++){
				print "   _complex_mul_star_assign((u).$cname\[$n\],(s).${cname}\[0\].$cname\[$i\],(s).${cname}\[2\].$cname\[$j\]); \\\n";
				print "   _complex_mul_star_assign((u).$cname\[$n\],(s).${cname}\[1\].$cname\[$i\],(s).${cname}\[3\].$cname\[$j\])";
				$n++;
				if ($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print ";\\\n"; }
			}
		}
	} else {
		print "   do { \\\n";
		print "      int _i,_j,_n=0;\\\n";
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		if ($N<(2*$unroll+1)) {
			my $n=0;
			for (my $j=0; $j<$N; $j++) {
				print "         _complex_mul_star_assign((u).$cname\[_n\],(s).${cname}\[0\].$cname\[_i\],(s).${cname}\[2\].$cname\[$j\]); \\\n";
				print "         _complex_mul_star_assign((u).$cname\[_n\],(s).${cname}\[1\].$cname\[_i\],(s).${cname}\[3\].$cname\[$j\]); \\\n";
				print "         ++_n; \\\n";
			}
		} else {
			print "         for (_j=0; _j<$vd; ){\\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            _complex_mul_star_assign((u).$cname\[_n\],(s).${cname}\[0\].$cname\[_i\],(s).${cname}\[2\].$cname\[_j\]); \\\n";
				print "            _complex_mul_star_assign((u).$cname\[_n\],(s).${cname}\[1\].$cname\[_i\],(s).${cname}\[3\].$cname\[_j\]); \\\n";
				print "            ++_n; ++_j; \\\n";
			}
			print "         }\\\n";
			for(my $i=0;$i<$vr;$i++){
				print "         _complex_mul_star_assign((u).$cname\[_n\],(s).${cname}\[0\].$cname\[_i\],(s).${cname}\[2\].$cname\[_j\]); \\\n";
				print "         _complex_mul_star_assign((u).$cname\[_n\],(s).${cname}\[1\].$cname\[_i\],(s).${cname}\[3\].$cname\[_j\]); \\\n";
				print "         ++_n; ++_j; \\\n";
			}
		}
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

sub write_suNr_FMAT {
  print "/* This fuction compute the hmc force matrix */\n";
  print "/* this fuction accumulates on the original matrix u */\n";
  print "#define _${rdataname}_FMAT(u,s) \\\n";
	if ($N<$Nmax) { #unroll all 
		my $n=0;
		for(my $i=0;$i<$N;$i++){
			for(my $j=0;$j<$N;$j++){
				print "   _complex_mul_star_assign_re((u).$cname\[$n\],(s).${cname}\[0\].$cname\[$i\],(s).${cname}\[2\].$cname\[$j\]); \\\n";
				print "   _complex_mul_star_assign_re((u).$cname\[$n\],(s).${cname}\[1\].$cname\[$i\],(s).${cname}\[3\].$cname\[$j\])";
				$n++;
				if ($i==$N-1 and $j==$N-1) { print "\n\n"; } else { print ";\\\n"; }
			}
		}
	} else {
		print "   do { \\\n";
		print "      int _i,_j,_n=0;\\\n";
		print "      for (_i=0; _i<$N; ++_i){\\\n";
		if ($N<(2*$unroll+1)) {
			my $n=0;
			for (my $j=0; $j<$N; $j++) {
				print "         _complex_mul_star_assign_re((u).$cname\[_n\],(s).${cname}\[0\].$cname\[_i\],(s).${cname}\[2\].$cname\[$j\]); \\\n";
				print "         _complex_mul_star_assign_re((u).$cname\[_n\],(s).${cname}\[1\].$cname\[_i\],(s).${cname}\[3\].$cname\[$j\]); \\\n";
				print "         ++_n; \\\n";
			}
		} else {
			print "         for (_j=0; _j<$vd; ){\\\n";
			for(my $i=0;$i<$unroll;$i++){
				print "            _complex_mul_star_assign_re((u).$cname\[_n\],(s).${cname}\[0\].$cname\[_i\],(s).${cname}\[2\].$cname\[_j\]); \\\n";
				print "            _complex_mul_star_assign_re((u).$cname\[_n\],(s).${cname}\[1\].$cname\[_i\],(s).${cname}\[3\].$cname\[_j\]); \\\n";
				print "            ++_n; ++_j; \\\n";
			}
			print "         }\\\n";
			for(my $i=0;$i<$vr;$i++){
				print "         _complex_mul_star_assign_re((u).$cname\[_n\],(s).${cname}\[0\].$cname\[_i\],(s).${cname}\[2\].$cname\[_j\]); \\\n";
				print "         _complex_mul_star_assign_re((u).$cname\[_n\],(s).${cname}\[1\].$cname\[_i\],(s).${cname}\[3\].$cname\[_j\]); \\\n";
				print "         ++_n; ++_j; \\\n";
			}
		}
		print "      }\\\n";
		print "   } while(0) \n\n";
	}
}

#
# SPINOR OPERATIONS
#

sub write_spinor_copy {
  print "/*  r=s (r,s spinors) */\n";
  print "#define _spinor_copy_${suff}(r,s) \\\n";
  print "  (r)=(s)\n\n";
}

sub write_spinor_zero {
  print "/*  r=0  (r spinor) */\n";
  print "#define _spinor_zero_${suff}(r) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_zero_${suff}((r).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_g5 {
  print "/*  s=g5*r (r,s spinors, g5 matrix) */\n";
  print "#define _spinor_g5_${suff}(s,r) \\\n";
  for (my $k=0; $k<2; $k++){ # g5 acts only on 3,4 components
    print "  (s).$cname\[$k\]=(r).$cname\[$k\]; \\\n";
  }
  for (my $k=2; $k<4; $k++){ # g5 acts only on 3,4 components
    print "  _vector_minus_${suff}((s).$cname\[$k\],(r).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
  print "/*  r=g5*r (r,s spinors, g5 matrix) */\n";
  print "#define _spinor_g5_assign_${suff}(r) \\\n";
  for (my $k=2; $k<4; $k++){ # g5 acts only on 3,4 components
    print "  _vector_minus_${suff}((r).$cname\[$k\],(r).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_minus {
  print "/*  s=-r (r,s spinors) */\n";
  print "#define _spinor_minus_${suff}(s,r) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_minus_${suff}((s).$cname\[$k\],(r).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mul {
  print "/*  r=k*s (k real; r,s spinors) */\n";
  print "#define _spinor_mul_${suff}(r,k,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_mul_${suff}((r).$cname\[$k\],k,(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mulc {
  print "/*  r=z*s (z complex; r,s spinors) */\n";
  print "#define _spinor_mulc_${suff}(r,z,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_mulc_${suff}((r).$cname\[$k\],z,(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mulc_add_assign {
  print "/*  r+=z*s (z complex; r,s spinors) */\n";
  print "#define _spinor_mulc_add_assign_${suff}(r,z,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_mulc_add_assign_${suff}((r).$cname\[$k\],(z),(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mul_add_assign {
  print "/*  r+=k*s (k real; r,s spinors) */\n";
  print "#define _spinor_mul_add_assign_${suff}(r,k,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_mul_add_assign_${suff}((r).$cname\[$k\],(k),(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_lc {
  print "/*  r=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */\n";
  print "#define _spinor_lc_${suff}(r,k1,s1,k2,s2) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_lc_${suff}((r).$cname\[$k\],(k1),(s1).$cname\[$k\],(k2),(s2).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_lc_add_assign {
  print "/*  r+=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */\n";
  print "#define _spinor_lc_add_assign_${suff}(r,k1,s1,k2,s2) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_lc_add_assign_${suff}((r).$cname\[$k\],(k1),(s1).$cname\[$k\],(k2),(s2).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_clc {
  print "/*  r=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */\n";
  print "#define _spinor_clc_${suff}(r,z1,s1,z2,s2) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_clc_${suff}((r).$cname\[$k\],(z1),(s1).$cname\[$k\],(z2),(s2).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_clc_add_assign {
  print "/*  r+=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */\n";
  print "#define _spinor_clc_add_assign_${suff}(r,z1,s1,z2,s2) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_clc_add_assign_${suff}((r).$cname\[$k\],(z1),(s1).$cname\[$k\],(z2),(s2).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_add {
  print "/*  r=s1+s2 (r,s1,s2 spinors) */\n";
  print "#define _spinor_add_${suff}(r,s1,s2) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_add_${suff}((r).$cname\[$k\],(s1).$cname\[$k\],(s2).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_sub {
  print "/*  r=s1-s2 (r,s1,s2 spinors) */\n";
  print "#define _spinor_sub_${suff}(r,s1,s2) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_sub_${suff}((r).$cname\[$k\],(s1).$cname\[$k\],(s2).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_add_assign {
  print "/*  r+=s (r,s spinors) */\n";
  print "#define _spinor_add_assign_${suff}(r,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_add_assign_${suff}((r).$cname\[$k\],(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_sub_assign {
  print "/*  r-=s (r,s spinors) */\n";
  print "#define _spinor_sub_assign_${suff}(r,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_sub_assign_${suff}((r).$cname\[$k\],(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_i_add_assign {
  print "/*  r+=i*s (r,s spinors) */\n";
  print "#define _spinor_i_add_assign_${suff}(r,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_i_add_assign_${suff}((r).$cname\[$k\],(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_i_sub_assign {
  print "/*  r-=i*s (r,s spinors) */\n";
  print "#define _spinor_i_sub_assign_${suff}(r,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_i_sub_assign_${suff}((r).$cname\[$k\],(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_prod_re {
  print "/* k=Real part of the scalar product r*s (r,s spinors) */\n";
  print "#define _spinor_prod_re_${suff}(k,r,s) \\\n";
  print "   do { \\\n";
  print "      _vector_prod_re_${suff}((k),(r).$cname\[0\],(s).$cname\[0\]);\\\n";
  for (my $k=1; $k<4; $k++){
    print "      _vector_prod_add_assign_re_${suff}((k),(r).$cname\[$k\],(s).$cname\[$k\]); \\\n";
  }
	print "   } while(0) \n\n";
}

sub write_spinor_prod_im {
  print "/* k=Im part of the scalar product r*s (r,s spinors) */\n";
  print "#define _spinor_prod_im_${suff}(k,r,s) \\\n";
  print "   do { \\\n";
  print "      _vector_prod_im_${suff}((k),(r).$cname\[0\],(s).$cname\[0\]);\\\n";
  for (my $k=1; $k<4; $k++){
    print "      _vector_prod_add_assign_im_${suff}((k),(r).$cname\[$k\],(s).$cname\[$k\]); \\\n";
  }
	print "   } while(0) \n\n";
}

sub write_spinor_prod {
  print "/* z=r*s (r,s spinors, z complex) */\n";
  print "#define _spinor_prod_${suff}(z,r,s) \\\n";	
  	print "   do { \\\n";
	print "      (z).re=0.;(z).im=0.; \\\n";
  for (my $k=0; $k<4; $k++){
    print "      _vector_prod_assign_${suff}((z),(r).$cname\[$k\],(s).$cname\[$k\])";
    if($k==3) {print"; \\\n   } while(0) \n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_prod_assign {
  print "/* z+=r*s (r,s spinors, z complex) */\n";
  print "#define _spinor_prod_assign_${suff}(z,r,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_prod_assign_${suff}((z),(r).$cname\[$k\],(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_g5_prod_re {
  print "/* k=Real part of the scalar product (g5*r)*s (r,s spinors) */\n";
  print "#define _spinor_g5_prod_re_${suff}(k,r,s) \\\n";
  print "   do { \\\n";
  print "      _vector_prod_re_${suff}((k),(r).$cname\[0\],(s).$cname\[0\]);\\\n";
  print "      _vector_prod_add_assign_re_${suff}((k),(r).$cname\[1\],(s).$cname\[1\]);\\\n";
  for (my $k=2; $k<4; $k++){
    print "      _vector_prod_sub_assign_re_${suff}((k),(r).$cname\[$k\],(s).$cname\[$k\]);\\\n";
  }
  print "   } while(0) \n\n";
}

sub write_spinor_g5_prod_im {
  print "/* k=Imaginary part of the scalar product (g5*r)*s (r,s spinors) */\n";
  print "#define _spinor_g5_prod_im_${suff}(k,r,s) \\\n";
  print "   do { \\\n";
  print "      _vector_prod_im_${suff}((k),(r).$cname\[0\],(s).$cname\[0\]);\\\n";
  print "      _vector_prod_add_assign_im_${suff}((k),(r).$cname\[1\],(s).$cname\[1\]);\\\n";
  for (my $k=2; $k<4; $k++){
    print "      _vector_prod_sub_assign_im_${suff}((k),(r).$cname\[$k\],(s).$cname\[$k\]);\\\n";
  }
  print "   } while(0) \n\n";
}

sub write_spinor_project {
  print "/* r-=z*s (z complex; r,s spinors) */\n";
  print "#define _spinor_project_${suff}(r,z,s) \\\n";
  for (my $k=0; $k<4; $k++){
    print "  _vector_project_${suff}((r).$cname\[$k\],z,(s).$cname\[$k\])";
    if($k==3) {print"\n\n";} else {print "; \\\n"}
  }
}


sub write_spinor_pminus {
  print "/* r=(1-g0)/2 * s (r,s spinors) */\n";
  print "#define _spinor_pminus_${suff}(r,s) \\\n";
  print "  _vector_add_${suff}((r).$cname\[0\],(s).$cname\[0\],(s).$cname\[2\]); \\\n";
  print "  _vector_add_${suff}((r).$cname\[1\],(s).$cname\[1\],(s).$cname\[3\]); \\\n";
  print "  _vector_mul_${suff}((r).$cname\[0\],0.5,(r).$cname\[0\]); \\\n";
  print "  _vector_mul_${suff}((r).$cname\[1\],0.5,(r).$cname\[1\]); \\\n";
  print "  (r).$cname\[2\] = (r).$cname\[0\]; \\\n";
  print "  (r).$cname\[3\] = (r).$cname\[1\]\n\n";
}

sub write_spinor_pplus {
  print "/* r=(1+g0)/2 * s (r,s spinors) */\n";
  print "#define _spinor_pplus_${suff}(r,s) \\\n";
  print "  _vector_sub_${suff}((r).$cname\[0\],(s).$cname\[0\],(s).$cname\[2\]); \\\n";
  print "  _vector_sub_${suff}((r).$cname\[1\],(s).$cname\[1\],(s).$cname\[3\]); \\\n";
  print "  _vector_mul_${suff}((r).$cname\[0\],0.5,(r).$cname\[0\]); \\\n";
  print "  _vector_mul_${suff}((r).$cname\[1\],0.5,(r).$cname\[1\]); \\\n";
  print "  _vector_mul_${suff}((r).$cname\[2\],-1.,(r).$cname\[0\]); \\\n";
  print "  _vector_mul_${suff}((r).$cname\[3\],-1.,(r).$cname\[1\])\n\n";
}


sub write_vector_myrandom {
  print "/* random vector */\n";
  print "#define _vector_myrand(r) \\\n";
  for(my $i=1;$i<$N;$i++){
    print "   (r).$cname\[$i\].re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
    print "   (r).$cname\[$i\].im=1.0*((double)rand())/(double)RAND_MAX; \\\n";
  }
  print "   (r).$cname\[$N\].re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
  print "   (r).$cname\[$N\].im=1.0*((double)rand())/(double)RAND_MAX\n\n";

}

sub write_vector_iszero {
  print "/* check if zero vector */\n";
  print "#define _vector_iszero(r,e) \\\n";
  for(my $i=1;$i<$N;$i++){
    print "   fabs((r).$cname\[$i\].re)<(e) && fabs((r).$cname\[$i\].im)<(e) && \\\n";
  }
  print "   fabs((r).$cname\[$N\].re)<(e) && fabs((r).$cname\[$N\].im)<(e) \n\n";

}

sub write_su3_myrandom {
  print "/* random matrix */\n";
  print "#define _${dataname}_myrand(r) \\\n";
  for(my $i=0;$i<$N*$N;$i++){
	  print "   (r).$cname\[$i\].re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
	  print "   (r).$cname\[$i\].im=1.0*((double)rand())/(double)RAND_MAX";
		if($i==$N*$N-1) { print "\n\n" } else { print ";\\\n"; } 
	}
}


#
# SU2  OPERATIONS
#


sub old2_write_su2_multiply {
    print "/* SU(2) matrix u times SU(2) vector s */\n";
    print "/* r=u*s */\n";
    print "/* using quaternionic representations for SU(2) matrices */\n";
        
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_multiply(r,u,s) \\\n";
        print "   (r).$cname\[0\].re=(u).$cname\[0\]*(s).$cname\[0\].re-(u).$cname\[1\]*(s).$cname\[1\].im-(u).$cname\[2\]*(s).$cname\[1\].re-(u).$cname\[3\]*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[0\].im=(u).$cname\[0\]*(s).$cname\[0\].im+(u).$cname\[1\]*(s).$cname\[1\].re-(u).$cname\[2\]*(s).$cname\[1\].im+(u).$cname\[3\]*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[1\].re=(u).$cname\[0\]*(s).$cname\[1\].re-(u).$cname\[1\]*(s).$cname\[0\].im+(u).$cname\[2\]*(s).$cname\[0\].re+(u).$cname\[3\]*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[1\].im=(u).$cname\[0\]*(s).$cname\[1\].im+(u).$cname\[1\]*(s).$cname\[0\].re+(u).$cname\[2\]*(s).$cname\[0\].im-(u).$cname\[3\]*(s).$cname\[1\].re \n\n";
        
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_multiply(r,u,s) \\\n";
        print "   (r).$cname\[0\].re=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[0\].re+=2.*((u).$cname\[0\]*(u).$cname\[3\]+(u).$cname\[1\]*(u).$cname\[2\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[0\].re+=2.*((u).$cname\[2\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[1\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[0\].im=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[0\].im+=2.*((u).$cname\[0\]*(u).$cname\[3\]+(u).$cname\[1\]*(u).$cname\[2\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[0\].im+=2.*((u).$cname\[2\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[1\])*(s).$cname\[2\].im; \\\n";
        print "   (r).$cname\[1\].re=2.*((u).$cname\[1\]*(u).$cname\[2\]-(u).$cname\[0\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[1\].re+=((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[1\].re+=2.*((u).$cname\[0\]*(u).$cname\[2\]+(u).$cname\[1\]*(u).$cname\[3\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[1\].im=2.*((u).$cname\[1\]*(u).$cname\[2\]-(u).$cname\[0\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[1\].im+=((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[1\].im+=2.*((u).$cname\[0\]*(u).$cname\[2\]+(u).$cname\[1\]*(u).$cname\[3\])*(s).$cname\[2\].im; \\\n";
        print "   (r).$cname\[2\].re=2.*((u).$cname\[0\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[2\].re+=2.*((u).$cname\[1\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[2\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[2\].re+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[2\].im=2.*((u).$cname\[0\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[2\].im+=2.*((u).$cname\[1\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[2\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[2\].im+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[2\].im \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }
    
}

sub write_su2_multiply {
    print "/* SU(2) matrix u times SU(2) vector s */\n";
    print "/* r=u*s */\n";
    print "/* using quaternionic representations for SU(2) matrices */\n";
        
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_multiply(r,u,s) \\\n";
        print "   (r).$cname\[0\].re=(u).$cname\[0\]*(s).$cname\[0\].re-(u).$cname\[1\]*(s).$cname\[1\].im-(u).$cname\[2\]*(s).$cname\[1\].re-(u).$cname\[3\]*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[0\].im=(u).$cname\[0\]*(s).$cname\[0\].im+(u).$cname\[1\]*(s).$cname\[1\].re-(u).$cname\[2\]*(s).$cname\[1\].im+(u).$cname\[3\]*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[1\].re=(u).$cname\[0\]*(s).$cname\[1\].re-(u).$cname\[1\]*(s).$cname\[0\].im+(u).$cname\[2\]*(s).$cname\[0\].re+(u).$cname\[3\]*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[1\].im=(u).$cname\[0\]*(s).$cname\[1\].im+(u).$cname\[1\]*(s).$cname\[0\].re+(u).$cname\[2\]*(s).$cname\[0\].im-(u).$cname\[3\]*(s).$cname\[1\].re \n\n";
        
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_multiply(r,u,s) \\\n";
        print "   (r).$cname\[0\].re=((u).$cname\[0\]*(u).$cname\[3\]+(u).$cname\[1\]*(u).$cname\[2\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[0\].re+=((u).$cname\[2\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[1\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[0\].re+=(r).$cname\[0\].re; \\\n";
        print "   (r).$cname\[0\].re+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";

        print "   (r).$cname\[0\].im=((u).$cname\[0\]*(u).$cname\[3\]+(u).$cname\[1\]*(u).$cname\[2\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[0\].im+=((u).$cname\[2\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[1\])*(s).$cname\[2\].im; \\\n";
        print "   (r).$cname\[0\].im+=(r).$cname\[0\].im; \\\n";
        print "   (r).$cname\[0\].im+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";

        print "   (r).$cname\[1\].re=((u).$cname\[1\]*(u).$cname\[2\]-(u).$cname\[0\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[1\].re+=((u).$cname\[0\]*(u).$cname\[2\]+(u).$cname\[1\]*(u).$cname\[3\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[1\].re+=(r).$cname\[1\].re; \\\n";
        print "   (r).$cname\[1\].re+=((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[1\].re; \\\n";

        print "   (r).$cname\[1\].im=((u).$cname\[1\]*(u).$cname\[2\]-(u).$cname\[0\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[1\].im+=((u).$cname\[0\]*(u).$cname\[2\]+(u).$cname\[1\]*(u).$cname\[3\])*(s).$cname\[2\].im; \\\n";
        print "   (r).$cname\[1\].im+=(r).$cname\[1\].im; \\\n";
        print "   (r).$cname\[1\].im+=((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[1\].im; \\\n";

        print "   (r).$cname\[2\].re=((u).$cname\[0\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[2\].re+=((u).$cname\[1\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[2\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[2\].re+=(r).$cname\[2\].re; \\\n";
        print "   (r).$cname\[2\].re+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[2\].re; \\\n";

        print "   (r).$cname\[2\].im=((u).$cname\[0\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[2\].im+=((u).$cname\[1\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[2\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[2\].im+=(r).$cname\[2\].im; \\\n";
        print "   (r).$cname\[2\].im+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[2\].im \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }
    
}

sub write_su2_inverse_multiply {
    print "/* SU(2) matrix u^dagger times SU(2) vector s */\n";
    print "/* r=u^(dagger)*s */\n";
    print "/* using quaternionic representations for SU(2) matrices */\n";
    
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_inverse_multiply(r,u,s) \\\n";
        print "   (r).$cname\[0\].re=(u).$cname\[0\]*(s).$cname\[0\].re+(u).$cname\[1\]*(s).$cname\[1\].im+(u).$cname\[2\]*(s).$cname\[1\].re+(u).$cname\[3\]*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[0\].im=(u).$cname\[0\]*(s).$cname\[0\].im-(u).$cname\[1\]*(s).$cname\[1\].re+(u).$cname\[2\]*(s).$cname\[1\].im-(u).$cname\[3\]*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[1\].re=(u).$cname\[0\]*(s).$cname\[1\].re+(u).$cname\[1\]*(s).$cname\[0\].im-(u).$cname\[2\]*(s).$cname\[0\].re-(u).$cname\[3\]*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[1\].im=(u).$cname\[0\]*(s).$cname\[1\].im-(u).$cname\[1\]*(s).$cname\[0\].re-(u).$cname\[2\]*(s).$cname\[0\].im+(u).$cname\[3\]*(s).$cname\[1\].re \n\n";
        
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_inverse_multiply(r,u,s) \\\n";
        print "   (r).$cname\[0\].re=((u).$cname\[1\]*(u).$cname\[2\]-(u).$cname\[0\]*(u).$cname\[3\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[0\].re+=((u).$cname\[0\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[3\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[0\].re+=(r).$cname\[0\].re; \\\n";
        print "   (r).$cname\[0\].re+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[0\].re; \\\n";

        print "   (r).$cname\[0\].im=((u).$cname\[1\]*(u).$cname\[2\]-(u).$cname\[0\]*(u).$cname\[3\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[0\].im+=((u).$cname\[0\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[3\])*(s).$cname\[2\].im; \\\n";
        print "   (r).$cname\[0\].im+=(r).$cname\[0\].im; \\\n";
        print "   (r).$cname\[0\].im+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[0\].im; \\\n";

        print "   (r).$cname\[1\].re=((u).$cname\[0\]*(u).$cname\[3\]+(u).$cname\[1\]*(u).$cname\[2\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[1\].re+=((u).$cname\[1\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[2\])*(s).$cname\[2\].re; \\\n";
        print "   (r).$cname\[1\].re+=(r).$cname\[1\].re; \\\n";
        print "   (r).$cname\[1\].re+=((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[1\].re; \\\n";

        print "   (r).$cname\[1\].im=((u).$cname\[0\]*(u).$cname\[3\]+(u).$cname\[1\]*(u).$cname\[2\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[1\].im+=((u).$cname\[1\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[2\])*(s).$cname\[2\].im; \\\n";
        print "   (r).$cname\[1\].im+=(r).$cname\[1\].im; \\\n";
        print "   (r).$cname\[1\].im+=((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[1\].im; \\\n";

        print "   (r).$cname\[2\].re=((u).$cname\[2\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[1\])*(s).$cname\[0\].re; \\\n";
        print "   (r).$cname\[2\].re+=((u).$cname\[0\]*(u).$cname\[2\]+(u).$cname\[1\]*(u).$cname\[3\])*(s).$cname\[1\].re; \\\n";
        print "   (r).$cname\[2\].re+=(r).$cname\[2\].re; \\\n";
        print "   (r).$cname\[2\].re+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[2\].re; \\\n";

        print "   (r).$cname\[2\].im=((u).$cname\[2\]*(u).$cname\[3\]-(u).$cname\[0\]*(u).$cname\[1\])*(s).$cname\[0\].im; \\\n";
        print "   (r).$cname\[2\].im+=((u).$cname\[0\]*(u).$cname\[2\]+(u).$cname\[1\]*(u).$cname\[3\])*(s).$cname\[1\].im; \\\n";
        print "   (r).$cname\[2\].im+=(r).$cname\[2\].im; \\\n";
        print "   (r).$cname\[2\].im+=((u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\])*(s).$cname\[2\].im \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }
    
}



sub write_su2_dagger {
    print "/* u=v^dagger */\n";
    if ($N==2) {
    	print "#define _${dataname}_dagger(u,v) \\\n";
        print "   (u).$cname\[0\]=(v).$cname\[0\]; \\\n";
        print "   (u).$cname\[1\]=-(v).$cname\[1\]; \\\n";
        print "   (u).$cname\[2\]=-(v).$cname\[2\]; \\\n";
        print "   (u).$cname\[3\]=-(v).$cname\[3\] \n\n";
    } else {
    	print "#define _${basename}${repsuff}_dagger(u,v) _${basename}${fundsuff}_dagger((u),(v))\n\n";
    }
}    

sub write_su2_zero {
    print "/* u=0 */\n";
    if ($N==2) {
    	print "#define _${dataname}_zero(u) \\\n";
        print "   (u).$cname\[0\]=0.; \\\n";
        print "   (u).$cname\[1\]=0.; \\\n";
        print "   (u).$cname\[2\]=0.; \\\n";
        print "   (u).$cname\[3\]=0.\n\n";
    } else {
    	print "#define _${basename}${repsuff}_zero(u) _${basename}${fundsuff}_zero((u))\n\n";
    }
}   

sub write_su2_unit {
    print "/* u=1 */\n";
    if ($N==2) {
    	print "#define _${dataname}_unit(u) \\\n";
        print "   (u).$cname\[0\]=1.; \\\n";
        print "   (u).$cname\[1\]=0.; \\\n";
        print "   (u).$cname\[2\]=0.; \\\n";
        print "   (u).$cname\[3\]=0.\n\n";
    } else {
    	print "#define _${basename}${repsuff}_unit(u) _${basename}${fundsuff}_unit((u))\n\n";
    }
}  

sub write_su2_minus {
    print "/* u=-v */\n";
    if ($N==2) {
    	print "#define _${dataname}_minus(u,v) \\\n";
        print "   (u).$cname\[0\]=-(v).$cname\[0\]; \\\n";
        print "   (u).$cname\[1\]=-(v).$cname\[1\]; \\\n";
        print "   (u).$cname\[2\]=-(v).$cname\[2\]; \\\n";
        print "   (u).$cname\[3\]=-(v).$cname\[3\] \n\n";
    } else {
    	print "#define _${basename}${repsuff}_minus(u,v) _${basename}${fundsuff}_minus((u),(v))\n\n";
    }
}    

sub write_su2_add_assign{
    print "/* u+=v */\n";
    if ($N==2) {
    	print "#define _${dataname}_add_assign(u,v) \\\n";
        print "   (u).$cname\[0\]+=(v).$cname\[0\]; \\\n";
        print "   (u).$cname\[1\]+=(v).$cname\[1\]; \\\n";
        print "   (u).$cname\[2\]+=(v).$cname\[2\]; \\\n";
        print "   (u).$cname\[3\]+=(v).$cname\[3\] \n\n";
    } else {
    	print "#define _${basename}${repsuff}_add_assign(u,v) _${basename}${fundsuff}_add_assign((u),(v))\n\n";
    }
}    

sub write_su2_sub_assign{
    print "/* u-=v */\n";
    if ($N==2) {
    	print "#define _${dataname}_sub_assign(u,v) \\\n";
        print "   (u).$cname\[0\]-=(v).$cname\[0\]; \\\n";
        print "   (u).$cname\[1\]-=(v).$cname\[1\]; \\\n";
        print "   (u).$cname\[2\]-=(v).$cname\[2\]; \\\n";
        print "   (u).$cname\[3\]-=(v).$cname\[3\] \n\n";
    } else {
    	print "#define _${basename}${repsuff}_sub_assign(u,v) _${basename}${fundsuff}_sub_assign((u),(v))\n\n";
    }
}    

sub write_su2_mul{
    print "/* u=r*v (u,v mat, r real) */\n";
    if ($N==2) {
    	print "#define _${dataname}_mul(u,r,v) \\\n";
        print "   (u).$cname\[0\]=(r)*(v).$cname\[0\]; \\\n";
        print "   (u).$cname\[1\]=(r)*(v).$cname\[1\]; \\\n";
        print "   (u).$cname\[2\]=(r)*(v).$cname\[2\]; \\\n";
        print "   (u).$cname\[3\]=(r)*(v).$cname\[3\] \n\n";
    } else {
    	print "#define _${basename}${repsuff}_mul(u,r,v) _${basename}${fundsuff}_mul((u),(r),(v))\n\n";
    }
}    

sub write_su2_trace_re{
    print "/* k=Re Tr (u) */\n";
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_trace_re(k,u) \\\n";
        print "   (k)=2.*(u).$cname\[0\] \n\n";
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_trace_re(k,u) \\\n";
        print "   (k)=3.*(u).$cname\[0\]*(u).$cname\[0\]-(u).$cname\[1\]*(u).$cname\[1\]-(u).$cname\[2\]*(u).$cname\[2\]-(u).$cname\[3\]*(u).$cname\[3\] \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }
}    

sub write_su2_trace_im{
    print "/* k=Im Tr (u) */\n";
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_trace_im(k,u) \\\n";
        print "   (k)=0. \n\n";
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_trace_im(k,u) \\\n";
        print "   (k)=0. \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }    
}    

sub write_su2_times_su2{
    print "/* u=v*w */\n";
    if ($N==2) {
        print "#define _${dataname}_times_${dataname}(u,v,w) \\\n";
        print "   (u).$cname\[0\]=(v).$cname\[0\]*(w).$cname\[0\]-(v).$cname\[1\]*(w).$cname\[1\]-(v).$cname\[2\]*(w).$cname\[2\]-(v).$cname\[3\]*(w).$cname\[3\]; \\\n";
        print "   (u).$cname\[1\]=(v).$cname\[0\]*(w).$cname\[1\]+(v).$cname\[1\]*(w).$cname\[0\]+(v).$cname\[2\]*(w).$cname\[3\]-(v).$cname\[3\]*(w).$cname\[2\]; \\\n";
        print "   (u).$cname\[2\]=(v).$cname\[0\]*(w).$cname\[2\]-(v).$cname\[1\]*(w).$cname\[3\]+(v).$cname\[2\]*(w).$cname\[0\]+(v).$cname\[3\]*(w).$cname\[1\]; \\\n";
        print "   (u).$cname\[3\]=(v).$cname\[0\]*(w).$cname\[3\]+(v).$cname\[1\]*(w).$cname\[2\]-(v).$cname\[2\]*(w).$cname\[1\]+(v).$cname\[3\]*(w).$cname\[0\]\n\n";
    } else {
    	print "#define _${basename}${repsuff}_times_${basename}${repsuff}(u,v,w) _${basename}${fundsuff}_times_${basename}${fundsuff}((u),(v),(w))\n\n";
    }
}

sub write_su2_times_su2_dagger{
    print "/* u=v*w^+ */\n";
    if ($N==2) {
        print "#define _${dataname}_times_${dataname}_dagger(u,v,w) \\\n";
        print "   (u).$cname\[0\]=(v).$cname\[0\]*(w).$cname\[0\]+(v).$cname\[1\]*(w).$cname\[1\]+(v).$cname\[2\]*(w).$cname\[2\]+(v).$cname\[3\]*(w).$cname\[3\]; \\\n";
        print "   (u).$cname\[1\]=(v).$cname\[1\]*(w).$cname\[0\]-(v).$cname\[2\]*(w).$cname\[3\]+(v).$cname\[3\]*(w).$cname\[2\]-(v).$cname\[0\]*(w).$cname\[1\]; \\\n";
        print "   (u).$cname\[2\]=(v).$cname\[1\]*(w).$cname\[3\]+(v).$cname\[2\]*(w).$cname\[0\]-(v).$cname\[3\]*(w).$cname\[1\]-(v).$cname\[0\]*(w).$cname\[2\]; \\\n";
        print "   (u).$cname\[3\]=(v).$cname\[3\]*(w).$cname\[0\]-(v).$cname\[0\]*(w).$cname\[3\]-(v).$cname\[1\]*(w).$cname\[2\]+(v).$cname\[2\]*(w).$cname\[1\]\n\n";
    } else {
    	print "#define _${basename}${repsuff}_times_${basename}${repsuff}_dagger(u,v,w) _${basename}${fundsuff}_times_${basename}${fundsuff}_dagger((u),(v),(w))\n\n";
    }
}

sub write_su2_dagger_times_su2{
    print "/* u=v^+*w */\n";
    if ($N==2) {
        print "#define _${dataname}_dagger_times_${dataname}(u,v,w) \\\n";
        print "   (u).$cname\[0\]=(v).$cname\[0\]*(w).$cname\[0\]+(v).$cname\[1\]*(w).$cname\[1\]+(v).$cname\[2\]*(w).$cname\[2\]+(v).$cname\[3\]*(w).$cname\[3\]; \\\n";
        print "   (u).$cname\[1\]=(v).$cname\[0\]*(w).$cname\[1\]-(v).$cname\[1\]*(w).$cname\[0\]-(v).$cname\[2\]*(w).$cname\[3\]+(v).$cname\[3\]*(w).$cname\[2\]; \\\n";
        print "   (u).$cname\[2\]=(v).$cname\[0\]*(w).$cname\[2\]+(v).$cname\[1\]*(w).$cname\[3\]-(v).$cname\[2\]*(w).$cname\[0\]-(v).$cname\[3\]*(w).$cname\[1\]; \\\n";
        print "   (u).$cname\[3\]=(v).$cname\[0\]*(w).$cname\[3\]-(v).$cname\[1\]*(w).$cname\[2\]+(v).$cname\[2\]*(w).$cname\[1\]-(v).$cname\[3\]*(w).$cname\[0\]\n\n";
    } else {
    	print "#define _${basename}${repsuff}_dagger_times_${basename}${repsuff}(u,v,w) _${basename}${fundsuff}_dagger_times_${basename}${fundsuff}((u),(v),(w))\n\n";
    }
}

sub write_su2_sqnorm {
    print "/* k=| u |2 ) */\n";
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_sqnorm(k,u) \\\n";
        print "   (k)=2.*((u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\]) \n\n";
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_sqnorm(k,u) \\\n";
        print "   (k)=(u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\]; \\\n";
        print "   (k)*=3.*(k) \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }    
}

sub write_su2_sqnorm_m1 {
    print "/* k=| 1 - u |2 ) */\n";
    if ($N==2) { #fundamental representation
        print "#define _${dataname}_sqnorm_m1(k,u) \\\n";
        print "   (k)=2.*(1.+(u).$cname\[0\]*((u).$cname\[0\]-2.)+(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\]) \n\n";
    } elsif ($N==3) { #adjoint representation
        print "#define _${rdataname}_sqnorm_m1(k,u) \\\n";
        print "   (k)=(u).$cname\[0\]*(u).$cname\[0\]+(u).$cname\[1\]*(u).$cname\[1\]+(u).$cname\[2\]*(u).$cname\[2\]+(u).$cname\[3\]*(u).$cname\[3\]; \\\n";
        print "   (k)*=3.*(k)+2.; \\\n";
        print "   (k)+=3.-8.*(u).$cname\[0\]*(u).$cname\[0\] \n\n";
    } else {
        die("Undefined fermion representation in quaternionic code. Exiting...\n");
    }    
}

sub write_su2_exp {
    print "/* u=Exp(dt*iT[n]h[n]) */\n";
    print "/* dt real; T[n] are the generators; h[n] is an algebra vector */\n";
    print "/* tmp is a temporary real */\n";
    if ($N==2) {
        print "#define _${dataname}_exp(dt,h,u) \\\n";
        print "   _algebra_vector_sqnorm_g((u).$cname\[0\],(h)); \\\n";
        print "   (u).$cname\[0\]=sqrt((u).$cname\[0\]); \\\n";
        print "   (u).$cname\[1\]=sin((dt)*0.5*(u).$cname\[0\])/((u).$cname\[0\]); \\\n";
        print "   (u).$cname\[2\]=(h).$cname\[0\]*(u).$cname\[1\]; \\\n";
        print "   (u).$cname\[3\]=(h).$cname\[2\]*(u).$cname\[1\]; \\\n";
        print "   (u).$cname\[1\]*=(h).$cname\[1\]; \\\n";
        print "   (u).$cname\[0\]=cos((dt)*0.5*(u).$cname\[0\]) \n\n";
    } else {
    	print "#define _${basename}${repsuff}_exp(dt,h,u) _${basename}${fundsuff}_exp((dt),(h),(u))\n\n";
    }
}

sub  write_read_spinor_gpu {
    my $i;
    print "/* Read spinor field component from GPU memory */\n";
    print "/* (output) v = ${dataname}_vector ; (input) in = ${dataname}_spinor* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 spinor component; */\n"; 
    print "#define _${rdataname}_read_spinor_flt_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$N)*(stride); \\\n";
    for($i=0; $i<$N-1; $i++) {
        print "      (v).c\[$i\]=((complex_flt*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\]=((complex_flt*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";

    print "#define _${rdataname}_read_spinor_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$N)*(stride); \\\n";
    for($i=0; $i<$N-1; $i++) {
        print "      (v).c\[$i\]=((complex*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\]=((complex*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";

}

sub  write_write_spinor_gpu {
    my $i;
    print "/* Write spinor field component to GPU memory */\n";
    print "/* (input) v = ${dataname}_vector ; (output) out = ${dataname}_spinor* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 spinor component; */\n"; 
    print "#define _${rdataname}_write_spinor_flt_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$N)*(stride); \\\n";
    for($i=0; $i<$N-1; $i++) {
        print "      ((complex_flt*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
    }
    print "      ((complex_flt*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${rdataname}_write_spinor_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$N)*(stride); \\\n";
    for($i=0; $i<$N-1; $i++) {
        print "      ((complex*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
    }
    print "      ((complex*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
    print "   } while (0) \n\n";
    
}

sub write_su2_read_gpu {
    print "/* Read an suN matrix from GPU memory */\n";
    print "/* (output) v = suN ; (input) in = suN* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    if ($N==2) { #fundamental representation
        my $i; 
        my $dim=4; #real components
        print "#define _${dataname}_flt_read_gpu(stride,v,in,iy,x) \\\n";
        print "   do {  \\\n";
        print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
        for($i=0; $i<$dim-1; $i++) {
            print "      (v).c\[$i\]=((float*)(in))\[__iz\]; __iz+=(stride); \\\n";
        }
        print "      (v).c\[$i\]=((float*)(in))\[__iz\]; \\\n";
        print "   } while (0) \n\n";

        print "#define _${dataname}_read_gpu(stride,v,in,iy,x) \\\n";
        print "   do {  \\\n";
        print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
        for($i=0; $i<$dim-1; $i++) {
            print "      (v).c\[$i\]=((double*)(in))\[__iz\]; __iz+=(stride); \\\n";
        }
        print "      (v).c\[$i\]=((double*)(in))\[__iz\]; \\\n";
        print "   } while (0) \n\n";
    
    } else {
        print "#define _${basename}${repsuff}_flt_read_gpu(stride,v,in,iy,x) _${basename}${fundsuff}_flt_read_gpu(stride,v,in,iy,x)\n\n";
        print "#define _${basename}${repsuff}_read_gpu(stride,v,in,iy,x) _${basename}${fundsuff}_read_gpu(stride,v,in,iy,x)\n\n";
    }

}

sub write_su2_write_gpu {
    print "/* Write an suN matrix to GPU memory */\n";
    print "/* (input) v = suN ; (output) in = suN* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    if ($N==2) { #fundamental representation
        my $i; 
        my $dim=4; #real components
        print "#define _${dataname}_flt_write_gpu(stride,v,out,iy,x) \\\n";
        print "   do {  \\\n";
        print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
        for($i=0; $i<$dim-1; $i++) {
            print "      ((float*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
        }
        print "      ((float*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
        print "   } while (0) \n\n";
        
        print "#define _${dataname}_write_gpu(stride,v,out,iy,x) \\\n";
        print "   do {  \\\n";
        print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
        for($i=0; $i<$dim-1; $i++) {
            print "      ((double*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
        }
        print "      ((double*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
        print "   } while (0) \n\n";
        
    } else {
        print "#define _${basename}${repsuff}_flt_write_gpu(stride,v,in,iy,x) _${basename}${fundsuff}_flt_write_gpu(stride,v,in,iy,x)\n\n";
        print "#define _${basename}${repsuff}_write_gpu(stride,v,in,iy,x) _${basename}${fundsuff}_write_gpu(stride,v,in,iy,x)\n\n";
    }
    
}

sub write_suN_read_gpu {
    my $i; 
    my $dim=$N*$N; #complex components
    my $rdim=2*$dim; #real components

    print "/* Read an suN matrix from GPU memory */\n";
    print "/* (output) v = suN ; (input) in = suN* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    
    print "#define _${dataname}_flt_read_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$rdim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      (v).c\[$i\].re=((float*)(in))\[__iz\]; __iz+=(stride); \\\n";
        print "      (v).c\[$i\].im=((float*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\].re=((float*)(in))\[__iz\]; __iz+=(stride); \\\n";
    print "      (v).c\[$i\].im=((float*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${dataname}_read_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$rdim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      (v).c\[$i\].re=((double*)(in))\[__iz\]; __iz+=(stride); \\\n";
        print "      (v).c\[$i\].im=((double*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\].re=((double*)(in))\[__iz\]; __iz+=(stride); \\\n";
    print "      (v).c\[$i\].im=((double*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";
    
}

sub write_suN_write_gpu {
    my $i; 
    my $dim=$N*$N; #complex components
    my $rdim=2*$dim; #real components
    
    print "/* Write an suN matrix to GPU memory */\n";
    print "/* (input) v = suN ; (output) out = suN* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    
    print "#define _${dataname}_flt_write_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$rdim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((float*)(out))\[__iz\]=(v).c\[$i\].re; __iz+=(stride); \\\n";
        print "      ((float*)(out))\[__iz\]=(v).c\[$i\].im; __iz+=(stride); \\\n";
    }
    print "      ((float*)(out))\[__iz\]=(v).c\[$i\].re; __iz+=(stride); \\\n";
    print "      ((float*)(out))\[__iz\]=(v).c\[$i\].im; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${dataname}_write_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$rdim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((double*)(out))\[__iz\]=(v).c\[$i\].re; __iz+=(stride); \\\n";
        print "      ((double*)(out))\[__iz\]=(v).c\[$i\].im; __iz+=(stride); \\\n";
    }
    print "      ((double*)(out))\[__iz\]=(v).c\[$i\].re; __iz+=(stride); \\\n";
    print "      ((double*)(out))\[__iz\]=(v).c\[$i\].im; \\\n";
    print "   } while (0) \n\n";
    
}


sub write_suNr_read_gpu {
    my $i; 
    my $dim=$N*$N; #real components
    
    print "/* Read an suN matrix from GPU memory */\n";
    print "/* (output) v = suN ; (input) in = suN* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    
    print "#define _${rdataname}_flt_read_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      (v).c\[$i\]=((float*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\]=((float*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${rdataname}_read_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      (v).c\[$i\]=((double*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\]=((double*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";
      
}

sub write_suNr_write_gpu {
    my $i; 
    my $dim=$N*$N; #real components
    
    print "/* Write an suN matrix to GPU memory */\n";
    print "/* (input) v = suN ; (output) out = suN* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    
    print "#define _${rdataname}_flt_write_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((float*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
    }
    print "      ((float*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${rdataname}_write_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((double*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
    }
    print "      ((double*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
    print "   } while (0) \n\n";
    
}

sub write_suN_av_read_gpu {
    my $i; 
    my $dim=$N*$N-1; #real components
    
    print "/* Read an suN algebra vector from GPU memory */\n";
    print "/* (output) v = suN_algebra_vector ; (input) in = suN_algebra_vector* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    
    print "#define _${rdataname}_av_flt_read_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      (v).c\[$i\]=((float*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\]=((float*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${rdataname}_av_read_gpu(stride,v,in,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      (v).c\[$i\]=((double*)(in))\[__iz\]; __iz+=(stride); \\\n";
    }
    print "      (v).c\[$i\]=((double*)(in))\[__iz\]; \\\n";
    print "   } while (0) \n\n";
    
}

sub write_suN_av_write_gpu {
    my $i; 
    my $dim=$N*$N-1; #real components
    
    print "/* Write an suN algebra vector to GPU memory */\n";
    print "/* (input) v = suN_algebra_vector ; (output) out = suN_algebra_vector* */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; */\n"; 
    
    print "#define _${rdataname}_av_flt_write_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((float*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
    }
    print "      ((float*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
    print "   } while (0) \n\n";
    
    print "#define _${rdataname}_av_write_gpu(stride,v,out,iy,x) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((double*)(out))\[__iz\]=(v).c\[$i\]; __iz+=(stride); \\\n";
    }
    print "      ((double*)(out))\[__iz\]=(v).c\[$i\]; \\\n";
    print "   } while (0) \n\n";
    
}

sub write_suN_av_mul_add_assign_gpu {
    my $i; 
    my $dim=$N*$N-1; #real components
    
    print "/* Mul_add_assign on a suN algebra vector on GPU  */\n";
    print "/* (in/out) v = suN_algebra_vector* ; (input) in = suN_algebra_vector */\n";
    print "/* (input) iy = site ; (input) x = 0..3 direction; (input) r = real */\n"; 
    
    print "#define _algebra_vector_mul_add_assign_gpu_${suff}_flt(stride,v,iy,x,r,in) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((float*)(v))\[__iz\]+=(in).c\[$i\]*(r); __iz+=(stride); \\\n";
    }
    print "      ((float*)(v))\[__iz\]+=(in).c\[$i\]*(r); \\\n";
    print "   } while (0) \n\n";

    print "#define _algebra_vector_mul_add_assign_gpu_${suff}(stride,v,iy,x,r,in) \\\n";
    print "   do {  \\\n";
    print "      int __iz=(iy)+((x)*$dim)*(stride); \\\n";
    for($i=0; $i<$dim-1; $i++) {
        print "      ((double*)(v))\[__iz\]+=(in).c\[$i\]*(r); __iz+=(stride); \\\n";
    }
    print "      ((double*)(v))\[__iz\]+=(in).c\[$i\]*(r); \\\n";
    print "   } while (0) \n\n";

    
}


