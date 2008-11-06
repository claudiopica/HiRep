TOPDIR = ../..

MAIN := check_puregauge_1

SRCS = $(addsuffix .c,$(MAIN))

default: $(MAIN)

MKDIR=$(TOPDIR)/Make
include $(MKDIR)/MkRules

MODDIR := Geometry
MODS := geometry_init  geometry_mpi_eo test_geometry_mpi geometry_mpi
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := logger
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := boundary_conditions #suN_utils suN_exp single_double_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Observables
# MODS := avr_plaquette mesons quark_propagator
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := random_su2 representation # cabmar staples update_hb Dphi Dphi_old random_momenta random_spinor_field find_spec_limits
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

$(MAIN): %: %.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(MOD_OBJS) &&\
	echo " done." || echo " failed."


.PHONY: default

cleanall: clean
	rm -rf $(MAIN) *~
