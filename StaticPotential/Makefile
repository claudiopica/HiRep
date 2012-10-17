TOPDIR = ..
MKDIR=$(TOPDIR)/Make

EXES = mk_wilsonloops monopoles poly_local poly_in_time polyakov_loops verify_group #tune_HYP_smearing

default: $(EXES)


MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo communications
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input endian
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := suN_utils time boundary_conditions HYP_smearing diag_hmat det_suNg
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette wilsonloops polyakov
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := cabmar staples update_hb random_su2 representation copy_gfield
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

#COMMON_OBJS = $(MOD_OBJS)

mk_wilsonloops_OBJS = $(MOD_OBJS)
monopoles_OBJS = $(MOD_OBJS)
polyakov_loops_OBJS = $(MOD_OBJS)
poly_in_time_OBJS = $(MOD_OBJS)
poly_local_OBJS = $(MOD_OBJS)
tune_HYP_smearing_OBJS = $(MOD_OBJS)
verify_group_OBJS = $(MOD_OBJS)


include $(MKDIR)/MkRules

.PHONY: default

