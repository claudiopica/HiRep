TOPDIR = ../..
MKDIR=$(TOPDIR)/Make

EXES = check_rotatedsf_1 check_rotatedsf_2 check_rotatedsf_3 check_rotatedsf_4 check_rotatedsf_5 check_rotatedsf_6 check_rotatedsf_7 check_rotatedsf_8 check_rotatedsf_9 check_rotatedsf_10 check_rotatedsf_11 check_rotatedsf_12

default: $(EXES)


MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo communications process_init
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := logger read_input archive_nocheck archive archive_spinors_sint endian
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := suN_utils suN_exp single_double_utils boundary_conditions time max_eig

MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields ranz2
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := Dphi staples random_su2 representation random_momenta random_spinor_field #pdate_rhmc update_hmc leapfrog local_action force0 copy_gfield find_spec_limits rational_functions force_rhmc O2MN_multistep cabmar update_hb
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra cg_mshift eva jacobi MINRES
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette sf_pcac
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

COMMON_OBJS = $(MOD_OBJS)

include $(MKDIR)/MkRules

.PHONY: default

