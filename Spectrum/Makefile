TOPDIR = ..

SRCS = mk_mesons.c random_cnfg.c

default: mk_mesons random_cnfg

MKDIR=$(TOPDIR)/Make
include $(MKDIR)/MkRules

MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo communications
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := boundary_conditions time #suN_utils suN_exp single_double_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette mesons pta_qprop
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss ranlxs ranlxd random_fields random_suNg 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := Dphi representation random_spinor_field random_su2
#staples update_hb  random_momenta update_rhmc leapfrog local_action force0 copy_gfield find_spec_limits rational_functions force_rhmc O2MN_multistep cabmar
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Statistics
MODS := 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra cg_mshift MINRES_mshift g5QMR_mshift MINRES #eva jacobi
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

mk_mesons: mk_mesons.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ mk_mesons.o $(MOD_OBJS) &&\
	echo " done." || echo " failed."

random_cnfg: random_cnfg.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ random_cnfg.o $(MOD_OBJS) &&\
	echo " done." || echo " failed."


.PHONY: default

