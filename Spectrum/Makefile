TOPDIR = ..

SRCS = qspectrum.c

default: qspectrum mk_mesons mk_prop

MKDIR=$(TOPDIR)/Make
include $(MKDIR)/MkRules

MODDIR := Geometry
MODS := geometry_init blocked eo_lexi test_geometry
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := suN_utils suN_exp single_double_utils boundary_conditions
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette quark_propagator mesons 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := cabmar staples update_hb random_su2 Dphi representation random_momenta random_spinor_field 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Statistics
MODS := base_stat jacknife
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra cg_mshift MINRES_mshift g5QMR_mshift MINRES
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

mk_prop: mk_prop.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ mk_prop.o $(MOD_OBJS) &&\
	echo " done." || echo " failed."

qspectrum: $(OBJS) $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(MOD_OBJS) &&\
	echo " done." || echo " failed."

mk_mesons: mk_mesons.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ mk_mesons.o $(MOD_OBJS) &&\
	echo " done." || echo " failed."

mk_dmesons: mk_dmesons.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ mk_dmesons.o $(MOD_OBJS) &&\
	echo " done." || echo " failed."

.PHONY: default

