TOPDIR = ../..

SRCS = check_diracoperator_3.c

MAIN = $(patsubst %.c, %, $(SRCS))

default:

MKDIR=$(TOPDIR)/Make
include $(MKDIR)/MkRules

MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo test_geometry_mpi blocked eo_lexi lexi communications
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := boundary_conditions
#MODS := suN_utils suN_exp single_double_utils boundary_conditions
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := 
#MODS := avr_plaquette
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := Dphi random_spinor_field random_su2 representation
#MODS := cabmar staples update_hb random_su2 Dphi representation random_momenta random_spinor_field update_rhmc leapfrog local_action force0 copy_gfield find_spec_limits rational_functions force_rhmc O2MN_multistep
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Statistics
MODS := 
#MODS := base_stat jacknife
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra cg_mshift
#MODS := linear_algebra cg_mshift MINRES_mshift g5QMR_mshift MINRES eva jacobi
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

$(MAIN): $(OBJS) $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(MOD_OBJS) &&\
	echo " done." || echo " failed."


.PHONY: default

