#!/usr/bin/env perl
use warnings;
use strict;

my $rootdir = ".";
require "$rootdir/Make/NinjaBuild.pl";

build_rules($rootdir);

###############################################################################
# LibHR
###############################################################################
{
  my $topdir = "LibHR";
  my @subdirs = ("Observables", "Error", "Geometry", "Inverters", "IO","Makefile", "Memory", "Random","Statistics", "Update", "Utils");

  my %exclude = ( Inverters => [ "HBiCGstab.c", "HBiCGstab_mshift.c", "BiCGstab_mshift.c", "dirac_eva.c", ], 
                  Observables => [ 
                                  #  "mesons.c", "measure_mesons.c", "measure_baryons.c", "z2semwall_new.c", "measure_mesons_ff.c", "z2semwall.c", "measure_scattering.c", "measure_renormalization.c",
                                  #  "loop_tools.c", "glueballs.c", "glueballs_op.c", "wilsonloops.c", "torellons.c",
                                  "trunc_hairpin.c",
                                ],
  );

  my @lib_objs =();
  for ( @subdirs ) {
      my $dir = "$topdir/$_";
      my @c_sources = glob($dir . "/*.c");
      exclude_files(\@c_sources, $exclude{$_});
      my @c_objs = obj_rules($dir, @c_sources);
      push(@lib_objs, @c_objs);
  }

  lib_rules("libhr.a", @lib_objs);
  print "build $topdir: phony libhr.a\n";
  print "default libhr.a\n";
}

my @libs = ("libhr.a"); #this is for later programs to use

###############################################################################
# HMC
###############################################################################
{
  my $topdir = "HMC";
  my %exes = (
    "hmc" => [ "hmc.c", "hmc_utils.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# Pure Gauge
###############################################################################
{
  my $topdir = "PureGauge";
  my %exes = (
    "suN" => [ "suN.c", "suN_utils.c", ],
    "suN_multilevel" => [ "suN_multilevel.c", "suN_utils_multilevel.c", ],
    "suN_multilevel_measure" => [ "suN_multilevel_measure.c", "suN_utils_multilevel.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# Spectrum
###############################################################################
{
  my $topdir = "Spectrum";
  my %exes = (
    "random_cnfg" => [ "random_cnfg.c", ],
    "random_spinor" => [ "random_spinor.c", ],
    "mk_mesons_with_z2semwall" => [ "mk_mesons_with_z2semwall.c", ],
    "mk_mesons_with_z2semwall_new" => [ "mk_mesons_with_z2semwall_new.c", ],
    "measure_formfactor" => [ "measure_formfactor.c", ],
    # "mk_sfcoupling" => [ "mk_sfcoupling.c", ],
    # "trunc_mesons" => [ "trunc_mesons.c", ],
    # "mk_mesons" => [ "mk_mesons.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  # print "default $topdir\n";
}

###############################################################################
# Tests
###############################################################################
{
  my $topdir = "TestProgram";
  my @subdirs = ( "PureGauge", "GaugeFix", "Inverters", "Random", 
                  "DiracOperator", "Algebra", "Disconnected", "Scattering", 
                  "Memory", "Sources", "RIMOM", "Deflate", "Utils", "Propagator", 
                  "Integrators", "StoredConfs", "SpinorField", "WilsonLoops", 
                  "Update", "Mesons", "Geometry", "LinearAlgebra",
              #    "RotatedSF", # this is broken 
  );

  my %exclude = ( Integrators => [ "check_integrator_utils_1.c", ], 
                  PureGauge => [ "check_puregauge_3.c", ], # this test is broken
                  Utils => [ "check_utils_3_gb_functions.c", "check_utils_3_tor_functions.c", # these 2 files are included directly in the main c test file
                            "check_complex.c", #this is broken
                  ], 
  );

  my %extra_sources = ( "check_integrator_1" => [ "check_integrator_utils_1.c"],
                        "check_update_1" => [ "../../HMC/hmc_utils.c" ],
                        "check_update_2" => [ "../../HMC/hmc_utils.c" ],
  );

  for my $dir ( @subdirs ) {
      my %exes = ();
      my @tests = glob("$topdir/$dir/*.c");
      @tests = map { local $_ = $_; s/^$topdir\/$dir\///; $_ } @tests;
      exclude_files(\@tests, $exclude{$dir});
      my @testnames = ();
      for ( @tests ) {
          my @c_sources = ( $_ );
          my $exe_name = $_; $exe_name =~ s{\.[^.]+$}{}; # remove extension
          if (exists $extra_sources{$exe_name}) { push(@c_sources,@{$extra_sources{$exe_name}}); }
          $exes{$exe_name} = [ @c_sources ] ;
          push(@testnames, "$topdir/$exe_name");
      }
      add_exes("$topdir/$dir", \%exes, \@libs);
  }

  @subdirs = map { "$topdir/$_" } @subdirs;
  print "build $topdir: phony @subdirs\n";
  print "default $topdir\n";
}

###############################################################################
# GaugeFix
###############################################################################
{
  my $topdir = "GaugeFix";
  my %exes = (
    "gaugefix_measure" => [ "gaugefix_measure.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# Converter
###############################################################################
{
  my $topdir = "Converter";
  my %exes = (
    "converter" => [ "converter.c", "archive_eolexi.c", "archive_milc.c", "archive_more_mpieo.c", 
                     "archive_ascii.c", "archive_fortran.c", "archive_openQCD.c", ],
    "converter_openQCD" => [ "converter_openQCD.c", "archive_eolexi.c", "archive_milc.c", "archive_more_mpieo.c", 
                             "archive_ascii.c", "archive_fortran.c", "archive_openQCD.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# Disconnected
###############################################################################
{
  my $topdir = "Disconnected";
  my %exes = (
    "compute_loops" => [ "compute_loops.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# ModeNumber
###############################################################################
{
  my $topdir = "ModeNumber";
  my %exes = (
    "mk_modenumber" => [ "mk_modenumber.c", ],
    "mk_eigvals" => [ "mk_eigvals.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  # print "default $topdir\n";
}

###############################################################################
# RenormalizationFactors
###############################################################################
{
  my $topdir = "RenormalizationFactors";
  my %exes = (
    "measure_Z_mom" => [ "measure_Z_mom.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# Reweight
###############################################################################
{
  my $topdir = "Reweight";
  my %exes = (
    "reweight_theta" => [ "reweight_theta.c", ],
    "reweight_mass" => [ "reweight_mass.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# Scattering
###############################################################################
{
  my $topdir = "Scattering";
  my %exes = (
    "scatter" => [ "scatter.c", ],
    "scattering_lengths" => [ "scattering_lengths.c", ],
    "sigma_triangle" => [ "sigma_triangle.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# StaticPotential
###############################################################################
{
  my $topdir = "StaticPotential";
  my %exes = (
    "mk_static_observables" => [ "mk_static_observables.c", ],
    # "tune_HYP_smearing" => [ "tune_HYP_smearing.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

###############################################################################
# WilsonFlow
###############################################################################
{
  my $topdir = "WilsonFlow";
  my %exes = (
    "WF_measure" => [ "WF_measure.c", ],
  );

  add_exes($topdir, \%exes, \@libs);
  print "default $topdir\n";
}

