#!/usr/bin/perl -w
use strict;

(@ARGV==2) or die("Usage: $0 Ng rep\n");

my ($Ng,$rep)=@ARGV;
my ($Nf,$c1,$c2);
$c1="C"; #gauge field always complex
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

my ($N,$suff,$complex,$to);

my $dataname;
my $rdataname;
my $cname="c";
my $structdef="typedef struct\n{\n";

open STDOUT, ">suN_types.h";

write_prolog_suN_types();

print "#define NG $Ng\n";
#system("./write_suN_def.pl $Ng g $c1 T");
write_suN_h($Ng,"g",$c1,"T");

print "#define NF $Nf\n";
#system("./write_suN_def.pl $Nf f $c2 T");
write_suN_h($Nf,"f",$c2,"T");

write_epilog();

open STDOUT, ">suN.h";

write_prolog_suN();

#system("./write_suN_def.pl $Ng g $c1 O");
write_suN_h($Ng,"g",$c1,"O");

#system("./write_suN_def.pl $Nf f $c2 O");
write_suN_h($Nf,"f",$c2,"O");

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

if ((not ($complex eq "C")) and (not ($complex eq "R"))) {
    die("Error: type must be C or R!\n");
}
if ((not ($to eq "T")) and (not ($to eq "O"))) {
    die("Error: Specify T for data types or O for operations!\n");
}

$dataname="suN$suff";
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
write_suN();
if ($complex eq "R") {
    write_suNr();
} else {
    print "typedef $dataname ${dataname}c;\n\n";
    print "typedef ${dataname}_dble ${dataname}c_dble;\n\n";
}
write_spinor();
write_suN_algebra_vector();

} else {

print <<END
/*******************************************************************************
*
* The following macros are the same for single and double precision types
*
* Depending on the macro, arguments are variables of type suN_vector and suN
* (or suN_vector_dble and suN_dble)
*
*******************************************************************************/

END
;
#write_vector_copy();
write_vector_zero();
write_vector_minus();
write_vector_mul();
write_vector_mulc();
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
write_vector_project();
if ($complex eq "R") {
    write_suNr_multiply();
    write_suNr_inverse_multiply();
}
write_suN_multiply();
write_suN_inverse_multiply();
write_algebra_vector_mul_add_assign();
write_algebra_vector_mul();
write_algebra_vector_zero();
write_algebra_vector_sqnorm();

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
write_suN_dagger();
write_suN_times_suN();
write_suN_times_suN_dagger();
write_suN_zero();
write_suN_unit();
write_suN_minus();
#write_suN_copy();
write_suN_mul();
write_suN_add_assign();
write_suN_sub_assign();
write_suN_sqnorm();
write_suN_sqnorm_m1();
write_suN_trace_re();
write_suN_trace_im();
write_suN_2TA();
write_suN_TA();

write_suN_FMAT();

if ($complex eq "R") { # we only need these functions at the moment...
    write_suNr_zero();
    write_suNr_FMAT();
    write_suNr_unit();
    write_suNr_times_suNr();
    write_suNr_add_assign();
    write_suNr_sub_assign();
    write_suNr_mul();
    write_suNr_trace_re();
    write_suNr_sqnorm();
		write_suNr_minus();
}

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
write_spinor_prod_re();
write_spinor_prod_im();
write_spinor_prod_assign();
write_spinor_g5_prod_re();
write_spinor_project();

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

sub write_suN_vector {
  print $structdef;
  print "   complex ";
  for(my $i=1;$i<$N;$i++){
    print "$cname$i,";
  }
  print "$cname$N;\n";
  print "} ${rdataname}_vector;\n\n";
  print $structdef;
  print "   complex_dble ";
  for(my $i=1;$i<$N;$i++){
    print "$cname$i,";
  }
  print "$cname$N;\n";
  print "} ${rdataname}_vector_dble;\n\n";
}

sub write_suN_algebra_vector {
  print $structdef;
  print "   float ";
  for(my $i=1;$i<($N*$N)-1;$i++){
    print "$cname$i,";
  }
  my $last=$N*$N-1;
  print "$cname$last;\n";
  print "} ${rdataname}_algebra_vector;\n\n";
  print $structdef;
  print "   double ";
  for(my $i=1;$i<($N*$N)-1;$i++){
    print "$cname$i,";
  }
  print "$cname$last;\n";
  print "} ${rdataname}_algebra_vector_dble;\n\n";
}

sub write_suN {
  print $structdef;
  print "   complex ";
  for(my $i=1;$i<$N;$i++){
    for(my $j=1;$j<$N;$j++){
      print "$cname$i\_$j, ";
    }
    print "$cname$i\_$N,\n";
    print "           ";
  }
  for(my $i=1;$i<$N;$i++){
    print "$cname$N\_$i, ";
  }
  print "$cname$N\_$N;\n";
  print "} $dataname;\n\n";
  print $structdef;
  print "   complex_dble ";
  for(my $i=1;$i<$N;$i++){
    for(my $j=1;$j<$N;$j++){
      print "$cname$i\_$j, ";
    }
    print "$cname$i\_$N,\n";
    print "                ";
  }
  for(my $i=1;$i<$N;$i++){
    print "$cname$N\_$i, ";
  }
  print "$cname$N\_$N;\n";
  print "} ${dataname}_dble;\n\n";
}

sub write_suNr {
  print $structdef;
  print "   float ";
  for(my $i=1;$i<$N;$i++){
    for(my $j=1;$j<$N;$j++){
      print "$cname$i\_$j, ";
    }
    print "$cname$i\_$N,\n";
    print "         ";
  }
  for(my $i=1;$i<$N;$i++){
    print "$cname$N\_$i, ";
  }
  print "$cname$N\_$N;\n";
  print "} $rdataname;\n\n";
  print $structdef;
  print "   double ";
  for(my $i=1;$i<$N;$i++){
    for(my $j=1;$j<$N;$j++){
      print "$cname$i\_$j, ";
    }
    print "$cname$i\_$N,\n";
    print "          ";
  }
  for(my $i=1;$i<$N;$i++){
    print "$cname$N\_$i, ";
  }
  print "$cname$N\_$N;\n";
  print "} ${rdataname}_dble;\n\n";
}

sub write_spinor {
  print $structdef;
  print "   ${rdataname}_vector ";
  my $slen=4;
  for(my $i=1;$i<$slen;$i++){
    print "$cname$i, ";
  }
  print "$cname$slen;\n";
  print "} ${rdataname}_spinor;\n\n";
  print $structdef;
  print "   ${rdataname}_vector_dble ";
  for(my $i=1;$i<$slen;$i++){
    print "$cname$i, ";
  }
  print "$cname$slen;\n";
  print "} ${rdataname}_spinor_dble;\n\n";
}

sub write_vector_copy {
  print "/* r=s */\n";
  print "#define _vector_copy_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   (r).$cname$i.re=(s).$cname$i.re; \\\n";
    print "   (r).$cname$i.im=(s).$cname$i.im";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_zero {
  print "/* r=0 */\n";
  print "#define _vector_zero_${suff}(r) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   (r).$cname$i.re=0.; \\\n";
    print "   (r).$cname$i.im=0.";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_minus {
  print "/* r=-s */\n";
  print "#define _vector_minus_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_minus((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_mul {
  print "/* r=c*s (c real) */\n";
  print "#define _vector_mul_${suff}(r,c,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_mulr((r).$cname$i,(c),(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_mulc {
  print "/* r=c*s (c complex) */\n";
  print "#define _vector_mulc_${suff}(r,c,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_mul((r).$cname$i,(c),(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_add {
  print "/* r=s1+s2 */\n";
  print "#define _vector_add_${suff}(r,s1,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_add((r).$cname$i,(s1).$cname$i,(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_sub {
  print "/* r=s1-s2 */\n";
  print "#define _vector_sub_${suff}(r,s1,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_sub((r).$cname$i,(s1).$cname$i,(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_i_add {
  print "/* r=s1+i*s2 */\n";
  print "#define _vector_i_add_${suff}(r,s1,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_i_add((r).$cname$i,(s1).$cname$i,(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_i_sub {
  print "/* r=s1-i*s2 */\n";
  print "#define _vector_i_sub_${suff}(r,s1,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_i_sub((r).$cname$i,(s1).$cname$i,(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_add_assign {
  print "/* r+=s */\n";
  print "#define _vector_add_assign_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_add_assign((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_sub_assign {
  print "/* r-=s */\n";
  print "#define _vector_sub_assign_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_sub_assign((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_i_add_assign {
  print "/* r+=i*s */\n";
  print "#define _vector_i_add_assign_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_i_add_assign((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_i_sub_assign {
  print "/* r-=i*s */\n";
  print "#define _vector_i_sub_assign_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_i_sub_assign((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_prod_re {
  print "/* Re(r*s) */\n";
  print "#define _vector_prod_re_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    my $par=" ";
    if ($i==1) { $par ="("; }
    print "   ${par}_complex_prod_re((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print ")\n\n"; } else { print "+ \\\n"; }
  }
}

sub write_vector_prod_im {
  print "/* Im(r*s) */\n";
  print "#define _vector_prod_im_${suff}(r,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    my $par=" ";
    if ($i==1) { $par ="("; }
    print "   ${par}_complex_prod_im((r).$cname$i,(s).$cname$i)";
    if($i==$N) { print ")\n\n"; } else { print "+ \\\n"; }
  }
}

sub write_vector_mulc_add_assign {
  print "/* r+=z*s (z complex) */\n";
  print "#define _vector_mulc_add_assign_${suff}(r,z,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_mul_assign((r).$cname$i,(z),(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_mul_add_assign {
  print "/* r+=c*s (c real) */\n";
  print "#define _vector_mul_add_assign_${suff}(r,c,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_mulr_assign((r).$cname$i,(c),(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_algebra_vector_mul_add_assign {
  print "/* r+=c*s (c real) */\n";
  print "#define _algebra_vector_mul_add_assign_${suff}(r,c,s) \\\n";
  my $last=$N*$N-1;
  for(my $i=1;$i<=$last;$i++){
    print "   (r).$cname$i+=(c)*(s).$cname$i";
    if($i==$last) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_algebra_vector_mul {
  print "/* r=c*s (c real) */\n";
  print "#define _algebra_vector_mul_${suff}(r,c,s) \\\n";
  my $last=$N*$N-1;
  for(my $i=1;$i<=$last;$i++){
    print "   (r).$cname$i=(c)*(s).$cname$i";
    if($i==$last) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_algebra_vector_zero {
  print "/* r=0  */\n";
  print "#define _algebra_vector_zero_${suff}(r) \\\n";
  my $last=$N*$N-1;
  for(my $i=1;$i<=$last;$i++){
    print "   (r).$cname$i=0.";
    if($i==$last) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_algebra_vector_sqnorm {
  print "/* |v|^2  */\n";
  print "#define _algebra_vector_sqnorm_${suff}(r) \\\n";
  my $last=$N*$N-1;
  print "    (";
  for(my $i=1;$i<=$last;$i++){
    print "((r).$cname$i*(r).$cname$i)";
    if($i==$last) { print ")\n\n"; } else { print "+"; }
  }
}

sub write_vector_lc {
  print "/* r=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */\n";
  print "#define _vector_lc_${suff}(r,k1,s1,k2,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_rlc((r).$cname$i,(k1),(s1).$cname$i,(k2),(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_lc_add_assign {
  print "/* r+=k1*s1+k2*s2 (k1,k2 real, s1,s2 vectors) */\n";
  print "#define _vector_lc_add_assign_${suff}(r,k1,s1,k2,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_rlc_assign((r).$cname$i,(k1),(s1).$cname$i,(k2),(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_clc {
  print "/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */\n";
  print "#define _vector_clc_${suff}(r,z1,s1,z2,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_clc((r).$cname$i,(z1),(s1).$cname$i,(z2),(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_clc_add_assign {
  print "/* r=z1*s1+z2*s2 (z1,z2 complex, s1,s2 vectors) */\n";
  print "#define _vector_clc_add_assign_${suff}(r,z1,s1,z2,s2) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_clc_assign((r).$cname$i,(z1),(s1).$cname$i,(z2),(s2).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_prod_assign {
  print "/* c+=r*s (c complex) */\n";
  print "#define _vector_prod_assign_${suff}(c,r,s) \\\n";
  my $var="   (c).re+=_vector_prod_re_${suff}(r,s); \\\n";
  print $var;
  $var="   (c).im+=_vector_prod_im_${suff}(r,s)\n\n";
  print $var;
}

sub write_vector_project {
  print "/* r-=z*s (z complex) */\n";
  print "#define _vector_project_${suff}(r,z,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    print "   _complex_mul_sub_assign((r).$cname$i,(z),(s).$cname$i)";
    if($i==$N) { print "\n\n"; } else { print "; \\\n"; }
  }
}

sub write_vector_cross_prod {
  print <<END

/*
* v.c1=(w.c2*z.c3-w.c3*z.c2)^*
* v.c2=(w.c3*z.c1-w.c1*z.c3)^*
* v.c3=(w.c1*z.c2-w.c2*z.c1)^*
* SU(3) only! (used in su3_fcts.c)
*/

#define _vector_cross_prod(v,w,z) \\
   (v).c1.re= (w).c2.re*(z).c3.re-(w).c2.im*(z).c3.im  \\
             -(w).c3.re*(z).c2.re+(w).c3.im*(z).c2.im; \\
   (v).c1.im= (w).c3.re*(z).c2.im+(w).c3.im*(z).c2.re  \\
             -(w).c2.re*(z).c3.im-(w).c2.im*(z).c3.re; \\
   (v).c2.re= (w).c3.re*(z).c1.re-(w).c3.im*(z).c1.im  \\
             -(w).c1.re*(z).c3.re+(w).c1.im*(z).c3.im; \\
   (v).c2.im= (w).c1.re*(z).c3.im+(w).c1.im*(z).c3.re  \\
             -(w).c3.re*(z).c1.im-(w).c3.im*(z).c1.re; \\
   (v).c3.re= (w).c1.re*(z).c2.re-(w).c1.im*(z).c2.im  \\
             -(w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \\
   (v).c3.im= (w).c2.re*(z).c1.im+(w).c2.im*(z).c1.re  \\
             -(w).c1.re*(z).c2.im-(w).c1.im*(z).c2.re


END
}

sub write_suN_multiply {
  print "/* SU(N) matrix u times SU(N) vector s */\n";
  print "/* r=u*s */\n";
  print "#define _${dataname}_multiply(r,u,s) \\\n";
  my $var;
  my $i=1;
  for(;$i<=$N;$i++){
    $var="   (r).$cname$i.re=";
    print $var;
    $var=~s/./ /g;
    my $j;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$i\_$j.re*(s).$cname$j.re-(u).$cname$i\_$j.im*(s).$cname$j.im";
      if($j==$N) {print "; \\\n";} else {print "  \\\n";}
    }
    $var="   (r).$cname$i.im=";
    print $var;
    $var=~s/./ /g;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$i\_$j.re*(s).$cname$j.im+(u).$cname$i\_$j.im*(s).$cname$j.re";
      if($i!=$N) {
	if($j==$N) {print "; \\\n";} else {print "  \\\n";}
      } else {
	if($j!=$N) {print "  \\\n";}	
      }
    }
  }
  print "\n\n";
}

sub write_suNr_multiply {
  print "/* SU(N) matrix u times SU(N) vector s */\n";
  print "/* r=u*s */\n";
  print "#define _${rdataname}_multiply(r,u,s) \\\n";
  my $var;
  my $i=1;
  for(;$i<=$N;$i++){
    $var="   (r).$cname$i.re=";
    print $var;
    $var=~s/./ /g;
    my $j;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$i\_$j*(s).$cname$j.re";
      if($j==$N) {print "; \\\n";} else {print "  \\\n";}
    }
    $var="   (r).$cname$i.im=";
    print $var;
    $var=~s/./ /g;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$i\_$j*(s).$cname$j.im";
      if($i!=$N) {
	if($j==$N) {print "; \\\n";} else {print "  \\\n";}
      } else {
	if($j!=$N) {print "  \\\n";}	
      }
    }
  }
  print "\n\n";
}

sub write_suN_inverse_multiply {
  print "/* SU(N) matrix u^dagger times SU(N) vector s */\n";
  print "/* r=(u^dagger)*s */\n";
  print "#define _${dataname}_inverse_multiply(r,u,s) \\\n";
  my $var;
  my $i=1;
  for(;$i<=$N;$i++){
    $var="   (r).$cname$i.re=";
    print $var;
    $var=~s/./ /g;
    my $j;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$j\_$i.re*(s).$cname$j.re+(u).$cname$j\_$i.im*(s).$cname$j.im";
      if($j==$N) {print "; \\\n";} else {print "  \\\n";}
    }
    $var="   (r).$cname$i.im=";
    print $var;
    $var=~s/./ /g;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$j\_$i.re*(s).$cname$j.im-(u).$cname$j\_$i.im*(s).$cname$j.re";
      if($i!=$N) {
	if($j==$N) {print "; \\\n";} else {print "  \\\n";}
      } else {
	if($j!=$N) {print "  \\\n";}	
      }
    }
  }
  print "\n\n";
}

sub write_suNr_inverse_multiply {
  print "/* SU(N) matrix u^dagger times SU(N) vector s */\n";
  print "/* r=(u^dagger)*s */\n";
  print "#define _${rdataname}_inverse_multiply(r,u,s) \\\n";
  my $var;
  my $i=1;
  for(;$i<=$N;$i++){
    $var="   (r).$cname$i.re=";
    print $var;
    $var=~s/./ /g;
    my $j;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$j\_$i*(s).$cname$j.re";
      if($j==$N) {print "; \\\n";} else {print "  \\\n";}
    }
    $var="   (r).$cname$i.im=";
    print $var;
    $var=~s/./ /g;
    for($j=1;$j<=$N;$j++){
      if($j==1) {print " ";} else {print $var,"+";}
      print "(u).$cname$j\_$i*(s).$cname$j.im";
      if($i!=$N) {
	if($j==$N) {print "; \\\n";} else {print "  \\\n";}
      } else {
	if($j!=$N) {print "  \\\n";}	
      }
    }
  }
  print "\n\n";
}

sub write_suN_dagger {
  print "/* u=v^dagger */\n";
  print "#define _${dataname}_dagger(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re= (v).$cname$j\_$i.re; \\\n";
      print "   (u).$cname$i\_$j.im=-(v).$cname$j\_$i.im";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_times_suN {
  print "/* u=v*w */\n";
  print "#define _${dataname}_times_${dataname}(u,v,w) \\\n";
  my $var;
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      $var="   (u).$cname$i\_$j.re=";
      print $var;
      $var=~s/./ /g;
      my $k;
      for ($k=1;$k<=$N; $k++){
	if($k==1) {print " ";} else {print $var,"+";}
	print "(v).$cname$i\_$k.re*(w).$cname$k\_$j.re-(v).$cname$i\_$k.im*(w).$cname$k\_$j.im";
	if($k==$N) {print"; \\\n";} else {print "  \\\n";}
      }
      $var="   (u).$cname$i\_$j.im=";
      print $var;
      $var=~s/./ /g;
      for ($k=1;$k<=$N; $k++){
	if($k==1) {print " ";} else {print $var,"+";}
	print "(v).$cname$i\_$k.re*(w).$cname$k\_$j.im+(v).$cname$i\_$k.im*(w).$cname$k\_$j.re";
	if($i==$N and $j==$N and $k==$N) {
	  print "\n\n";
	} else {
	  if($k==$N) {print"; \\\n";} else {print "  \\\n"}
	}
      }
    }
  }
}

sub write_suNr_times_suNr {
  print "/* u=v*w */\n";
  print "#define _${rdataname}_times_${rdataname}(u,v,w) \\\n";
  my $var;
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      $var="   (u).$cname$i\_$j=";
      print $var;
      $var=~s/./ /g;
      my $k;
      for ($k=1;$k<=$N; $k++){
	if($k==1) {print " ";} else {print $var,"+";}
	print "(v).$cname$i\_$k*(w).$cname$k\_$j";
	if($i==$N and $j==$N and $k==$N) {
	  print "\n\n";
	} else {
	  if($k==$N) {print"; \\\n";} else {print "  \\\n";}
	}
      }
    }
  }
}

sub write_suN_times_suN_dagger {
  print "/* u=v*w^+ */\n";
  print "#define _${dataname}_times_${dataname}_dagger(u,v,w) \\\n";
  my $var;
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      $var="   (u).$cname$i\_$j.re=";
      print $var;
      $var=~s/./ /g;
      my $k;
      for ($k=1;$k<=$N; $k++){
	if($k==1) {print " ";} else {print $var,"+";}
	print "(v).$cname$i\_$k.re*(w).$cname$j\_$k.re+(v).$cname$i\_$k.im*(w).$cname$j\_$k.im";
	if($k==$N) {print"; \\\n";} else {print "  \\\n";}
      }
      $var="   (u).$cname$i\_$j.im=";
      print $var;
      $var=~s/./ /g;
      for ($k=1;$k<=$N; $k++){
	if($k==1) {print " ";} else {print $var,"+";}
	print "(v).$cname$i\_$k.im*(w).$cname$j\_$k.re-(v).$cname$i\_$k.re*(w).$cname$j\_$k.im";
	if($i==$N and $j==$N and $k==$N) {
	  print "\n\n";
	} else {
	  if($k==$N) {print"; \\\n";} else {print "  \\\n"}
	}
      }
    }
  }
}

sub write_suN_zero {
  print "/* u=0 */\n";
  print "#define _${dataname}_zero(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re= 0.; \\\n";
      print "   (u).$cname$i\_$j.im= 0.";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_zero {
  print "/* u=0 */\n";
  print "#define _${rdataname}_zero(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j= 0.";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_unit {
  print "/* u=1 */\n";
  print "#define _${dataname}_unit(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      my $val = "0.";
      if ($i==$j) { $val="1."; }
      print "   (u).$cname$i\_$j.re= $val; \\\n";
      print "   (u).$cname$i\_$j.im= 0.";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_unit {
  print "/* u=1 */\n";
  print "#define _${rdataname}_unit(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      my $val = "0.";
      if ($i==$j) { $val="1."; }
      print "   (u).$cname$i\_$j= $val";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_copy {
  print "/* u=v */\n";
  print "#define _${dataname}_copy(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re = (v).$cname$i\_$j.re; \\\n";
      print "   (u).$cname$i\_$j.im = (v).$cname$i\_$j.im";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_minus {
  print "/* u=-v */\n";
  print "#define _${dataname}_minus(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re = -(v).$cname$i\_$j.re; \\\n";
      print "   (u).$cname$i\_$j.im = -(v).$cname$i\_$j.im";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_minus {
  print "/* u=-v */\n";
  print "#define _${rdataname}_minus(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j = -(v).$cname$i\_$j";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_mul {
  print "/* u=r*v (u,v mat, r real) */\n";
  print "#define _${dataname}_mul(u,r,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re = (r)*(v).$cname$i\_$j.re; \\\n";
      print "   (u).$cname$i\_$j.im = (r)*(v).$cname$i\_$j.im";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_mul {
  print "/* u=r*v (u,v mat, r real) */\n";
  print "#define _${rdataname}_mul(u,r,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j = (r)*(v).$cname$i\_$j";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_add_assign {
  print "/* u+=v */\n";
  print "#define _${dataname}_add_assign(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re += (v).$cname$i\_$j.re; \\\n";
      print "   (u).$cname$i\_$j.im += (v).$cname$i\_$j.im";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_add_assign {
  print "/* u+=v */\n";
  print "#define _${rdataname}_add_assign(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j += (v).$cname$i\_$j";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_sub_assign {
  print "/* u-=v */\n";
  print "#define _${dataname}_sub_assign(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j.re -= (v).$cname$i\_$j.re; \\\n";
      print "   (u).$cname$i\_$j.im -= (v).$cname$i\_$j.im";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_sub_assign {
  print "/* u-=v */\n";
  print "#define _${rdataname}_sub_assign(u,v) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   (u).$cname$i\_$j -= (v).$cname$i\_$j";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suN_sqnorm {
  print "/* returns: | u |2 ) */\n";
  print "#define _${dataname}_sqnorm(u) \\\n";
  print "   (";
  for(my $i=1;$i<=$N;$i++){
      for(my $j=1;$j<=$N;$j++){
	  print "   ((u).$cname$i\_$j.re)*((u).$cname$i\_$j.re) + \\\n";
	  print "   ((u).$cname$i\_$j.im)*((u).$cname$i\_$j.im) ";
	  if($j==$N and $i==$N) {print ")\n\n";} else {print "+ \\\n";}
      }
  }
}

sub write_suNr_sqnorm {
  print "/* returns: | u |2 ) */\n";
  print "#define _${rdataname}_sqnorm(u) \\\n";
  print "   (";
  for(my $i=1;$i<=$N;$i++){
      for(my $j=1;$j<=$N;$j++){
	  print "   ((u).$cname$i\_$j)*((u).$cname$i\_$j) ";
	  if($j==$N and $i==$N) {print ")\n\n";} else {print "+ \\\n";}
      }
  }
}

sub write_suN_sqnorm_m1 {
  print "/* returns: | 1 - u |2 ) */\n";
  print "#define _${dataname}_sqnorm_m1(u) \\\n";
  print "   (";
  for(my $i=1;$i<=$N;$i++){
      for(my $j=1;$j<=$N;$j++){
	  if ($i==$j) {
	      print "   (1.-(u).$cname$i\_$j.re)*(1.-(u).$cname$i\_$j.re) + \\\n";
	  } else {
	      print "   ((u).$cname$i\_$j.re)*((u).$cname$i\_$j.re) + \\\n";
	  }
	  print "   ((u).$cname$i\_$j.im)*((u).$cname$i\_$j.im) ";
	  if($j==$N and $i==$N) {print ")\n\n";} else {print "+ \\\n";}
      }
  }
}

sub write_suN_trace_re {
  print "/* returns: Re Tr (u) */\n";
  print "#define _${dataname}_trace_re(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
      my $par=" ";
      if ($i==1) { $par="("; }
      print "   $par(u).$cname$i\_$i.re";
      if($i==$N) {print ")\n\n";} else {print "+ \\\n";}
  }
}

sub write_suNr_trace_re {
  print "/* returns: Re Tr (u) */\n";
  print "#define _${rdataname}_trace_re(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
      my $par=" ";
      if ($i==1) { $par="("; }
      print "   $par(u).$cname$i\_$i";
      if($i==$N) {print ")\n\n";} else {print "+ \\\n";}
  }
}

sub write_suN_trace_im {
  print "/* returns: Re Im (u) */\n";
  print "#define _${dataname}_trace_im(u) \\\n";
  for(my $i=1;$i<=$N;$i++){
      my $par=" ";
      if ($i==1) { $par="("; }
      print "   $par(u).$cname$i\_$i.im";
      if($i==$N) {print ")\n\n";} else {print "+ \\\n";}
  }
}

sub write_suN_2TA {
  print "/* u=v - v^+ -1/N Tr(v - v^+)*I */\n";
  print "#define _${dataname}_2TA(u,v) \\\n";
  print "  {float _trim = _suNg_trace_im(v)*(2./$N.);\\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=$i;$j<=$N;$j++){
	if($i==$j) {
	    print "   (u).$cname$i\_$j.re= 0.; \\\n";
	    print "   (u).$cname$i\_$j.im= 2.*(v).$cname$j\_$i.im-_trim; \\\n";
	} else {
	    print "   (u).$cname$i\_$j.re= (v).$cname$i\_$j.re-(v).$cname$j\_$i.re; \\\n";
	    print "   (u).$cname$i\_$j.im= (v).$cname$j\_$i.im+(v).$cname$i\_$j.im; \\\n";
	    print "   (u).$cname$j\_$i.re= -(u).$cname$i\_$j.re; \\\n";
	    print "   (u).$cname$j\_$i.im= (u).$cname$i\_$j.im; \\\n";
	}
    }
  }
  print "  } \n\n";
}

sub write_suN_TA {
  print "/* u=0.5(v - v^+) -1/(2N) Tr(v - v^+)*I */\n";
  print "#define _${dataname}_TA(u,v) \\\n";
  print "  {float _trim = _suNg_trace_im(v)*(1./$N.);\\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=$i;$j<=$N;$j++){
	if($i==$j) {
	    print "   (u).$cname$i\_$j.re= 0.; \\\n";
	    print "   (u).$cname$i\_$j.im= (v).$cname$j\_$i.im-_trim; \\\n";
	} else {
	    print "   (u).$cname$i\_$j.re= 0.5*((v).$cname$i\_$j.re-(v).$cname$j\_$i.re); \\\n";
	    print "   (u).$cname$i\_$j.im= 0.5*((v).$cname$j\_$i.im+(v).$cname$i\_$j.im); \\\n";
	    print "   (u).$cname$j\_$i.re= -(u).$cname$i\_$j.re; \\\n";
	    print "   (u).$cname$j\_$i.im= (u).$cname$i\_$j.im; \\\n";
	}
    }
  }
  print "  } \n\n";
}

sub write_suN_FMAT {
  print "/* This fuction compute the hmc force matrix */\n";
  print "/* this fuction accumulates on the original matrix u */\n";
  print "#define _${dataname}_FMAT(u,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   _complex_mul_star_assign((u).$cname$i\_$j,(s).${cname}1.$cname$i,(s).${cname}3.$cname$j); \\\n";
      print "   _complex_mul_star_assign((u).$cname$i\_$j,(s).${cname}2.$cname$i,(s).${cname}4.$cname$j)";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_suNr_FMAT {
  print "/* This fuction compute the hmc force matrix */\n";
  print "/* this fuction accumulates on the original matrix u */\n";
  print "#define _${rdataname}_FMAT(u,s) \\\n";
  for(my $i=1;$i<=$N;$i++){
    for(my $j=1;$j<=$N;$j++){
      print "   _complex_mul_star_assign_re((u).$cname$i\_$j,(s).${cname}1.$cname$i,(s).${cname}3.$cname$j); \\\n";
      print "   _complex_mul_star_assign_re((u).$cname$i\_$j,(s).${cname}2.$cname$i,(s).${cname}4.$cname$j)";
      if($i==$N and $j==$N) {print "\n\n";} else {print "; \\\n";}
    }
  }
}

sub write_spinor_copy {
  print "/*  r=s (r,s spinors) */\n";
  print "#define _spinor_copy_${suff}(r,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_copy_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_zero {
  print "/*  r=0  (r spinor) */\n";
  print "#define _spinor_zero_${suff}(r) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_zero_${suff}((r).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_g5 {
  print "/*  s=g5*r (r,s spinors, g5 matrix) */\n";
  print "#define _spinor_g5_${suff}(s,r) \\\n";
  for (my $k=1; $k<3; $k++){ # g5 acts only on 3,4 components
    print "  (s).$cname$k=(r).$cname$k; \\\n";
  }
  for (my $k=3; $k<5; $k++){ # g5 acts only on 3,4 components
    print "  _vector_minus_${suff}((s).$cname$k,(r).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
  print "/*  r=g5*r (r,s spinors, g5 matrix) */\n";
  print "#define _spinor_g5_assign_${suff}(r) \\\n";
  for (my $k=3; $k<5; $k++){ # g5 acts only on 3,4 components
    print "  _vector_minus_${suff}((r).$cname$k,(r).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_minus {
  print "/*  s=-r (r,s spinors) */\n";
  print "#define _spinor_minus_${suff}(s,r) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_minus_${suff}((s).$cname$k,(r).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mul {
  print "/*  r=c*s (c real; r,s spinors) */\n";
  print "#define _spinor_mul_${suff}(r,c,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_mul_${suff}((r).$cname$k,c,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mulc {
  print "/*  r=c*s (c complex; r,s spinors) */\n";
  print "#define _spinor_mulc_${suff}(r,c,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_mulc_${suff}((r).$cname$k,c,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mulc_add_assign {
  print "/*  r+=c*s (c complex; r,s spinors) */\n";
  print "#define _spinor_mulc_add_assign_${suff}(r,c,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_mulc_add_assign_${suff}((r).$cname$k,(c),(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_mul_add_assign {
  print "/*  r+=k*s (k real; r,s spinors) */\n";
  print "#define _spinor_mul_add_assign_${suff}(r,k,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_mul_add_assign_${suff}((r).$cname$k,(k),(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_lc {
  print "/*  r=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */\n";
  print "#define _spinor_lc_${suff}(r,k1,s1,k2,s2) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_lc_${suff}((r).$cname$k,(k1),(s1).$cname$k,(k2),(s2).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_lc_add_assign {
  print "/*  r+=k1*s1+k2*s2 (k1,k2 real; r,s1,s2 spinors) */\n";
  print "#define _spinor_lc_add_assign_${suff}(r,k1,s1,k2,s2) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_lc_add_assign_${suff}((r).$cname$k,(k1),(s1).$cname$k,(k2),(s2).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_clc {
  print "/*  r=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */\n";
  print "#define _spinor_clc_${suff}(r,z1,s1,z2,s2) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_clc_${suff}((r).$cname$k,(z1),(s1).$cname$k,(z2),(s2).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_clc_add_assign {
  print "/*  r+=z1*s1+z2*s2 (z1,z2 complex; r,s1,s2 spinors) */\n";
  print "#define _spinor_clc_add_assign_${suff}(r,z1,s1,z2,s2) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_clc_add_assign_${suff}((r).$cname$k,(z1),(s1).$cname$k,(z2),(s2).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_add {
  print "/*  r=s1+s2 (r,s1,s2 spinors) */\n";
  print "#define _spinor_add_${suff}(r,s1,s2) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_add_${suff}((r).$cname$k,(s1).$cname$k,(s2).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_sub {
  print "/*  r=s1-s2 (r,s1,s2 spinors) */\n";
  print "#define _spinor_sub_${suff}(r,s1,s2) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_sub_${suff}((r).$cname$k,(s1).$cname$k,(s2).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_add_assign {
  print "/*  r+=s (r,s spinors) */\n";
  print "#define _spinor_add_assign_${suff}(r,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_add_assign_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_sub_assign {
  print "/*  r-=s (r,s spinors) */\n";
  print "#define _spinor_sub_assign_${suff}(r,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_sub_assign_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_prod_re {
  print "/* Real part of the scalar product r*s (r,s spinors) */\n";
  print "#define _spinor_prod_re_${suff}(r,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_prod_re_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "+ \\\n"}
  }
}

sub write_spinor_prod_im {
  print "/* Im part of the scalar product r*s (r,s spinors) */\n";
  print "#define _spinor_prod_im_${suff}(r,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_prod_im_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "+ \\\n"}
  }
}

sub write_spinor_prod_assign {
  print "/* c+=r*s (r,s spinors, c complex) */\n";
  print "#define _spinor_prod_assign_${suff}(c,r,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_prod_assign_${suff}((c),(r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}

sub write_spinor_g5_prod_re {
  print "/* Real part of the scalar product (g5*r)*s (r,s spinors) */\n";
  print "#define _spinor_g5_prod_re_${suff}(r,s) \\\n";
  for (my $k=1; $k<3; $k++){
    print "  _vector_prod_re_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==2) {print"- \\\n";} else {print "+ \\\n"}
  }
  for (my $k=3; $k<5; $k++){
    print "  _vector_prod_re_${suff}((r).$cname$k,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "- \\\n"}
  }
}

sub write_spinor_project {
  print "/* r-=z*s (z complex; r,s spinors) */\n";
  print "#define _spinor_project_${suff}(r,z,s) \\\n";
  for (my $k=1; $k<5; $k++){
    print "  _vector_project_${suff}((r).$cname$k,z,(s).$cname$k)";
    if($k==4) {print"\n\n";} else {print "; \\\n"}
  }
}


sub write_vector_myrandom {
  print "/* random vector */\n";
  print "#define _vector_myrand(r) \\\n";
  for(my $i=1;$i<$N;$i++){
    print "   (r).$cname$i.re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
    print "   (r).$cname$i.im=1.0*((double)rand())/(double)RAND_MAX; \\\n";
  }
  print "   (r).$cname$N.re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
  print "   (r).$cname$N.im=1.0*((double)rand())/(double)RAND_MAX\n\n";

}

sub write_vector_iszero {
  print "/* check if zero vector */\n";
  print "#define _vector_iszero(r,e) \\\n";
  for(my $i=1;$i<$N;$i++){
    print "   fabs((r).$cname$i.re)<(e) && fabs((r).$cname$i.im)<(e) && \\\n";
  }
  print "   fabs((r).$cname$N.re)<(e) && fabs((r).$cname$N.im)<(e) \n\n";

}

sub write_su3_myrandom {
  print "/* random matrix */\n";
  print "#define _${dataname}_myrand(r) \\\n";
  for(my $i=1;$i<$N;$i++){
      for(my $j=1; $j<=$N; $j++){ 
	  print "   (r).$cname$i$j.re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
	  print "   (r).$cname$i$j.im=1.0*((double)rand())/(double)RAND_MAX; \\\n";
      }
  }
  for(my $j=1; $j<$N; $j++){ 
      my $i=$N;
      print "   (r).$cname$i$j.re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
      print "   (r).$cname$i$j.im=1.0*((double)rand())/(double)RAND_MAX; \\\n";
  }
  print "   (r).$cname$N$N.re=1.0*((double)rand())/(double)RAND_MAX; \\\n";
  print "   (r).$cname$N$N.im=1.0*((double)rand())/(double)RAND_MAX\n\n";

}

