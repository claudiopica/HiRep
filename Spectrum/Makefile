TOPDIR = ..
MKDIR=$(TOPDIR)/Make

EXES = mk_mesons random_cnfg random_spinor mk_mesons_with_z2semwall mk_sfcoupling trunc_mesons

default: $(EXES)


MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo communications geometry_mpi_eo_patch
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input archive_spinors endian
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := boundary_conditions time max_eig single_double_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette trunc_hairpin mesons g5QMR_trunc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss ranlxs ranlxd ranz2 random_fields random_suNg 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := representation random_spinor_field random_su2 Dphi
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Statistics
MODS := 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra cg_mshift MINRES_mshift g5QMR_mshift MINRES eva jacobi
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Converter
MODS := archive_eolexi archive_henty
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

COMMON_OBJS = $(MOD_OBJS)

mk_mesons_with_z2semwall_OBJS = $(TOPDIR)/Observables/mesons.o $(TOPDIR)/Observables/pta_qprop.o $(TOPDIR)/Observables/z2semwall.o

mk_mesons_OBJS = $(TOPDIR)/Observables/mesons.o $(TOPDIR)/Observables/pta_qprop.o

mk_sfcoupling_OBJS = $(TOPDIR)/Observables/sf_action.o

include $(MKDIR)/MkRules

.PHONY: default

