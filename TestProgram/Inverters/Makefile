TOPDIR = ../..
MKDIR=$(TOPDIR)/Make

EXES = check_inverters_DDalphaAMG check_inverters_SAP # check_inverters_1 check_inverters_2 check_inverters_3 check_inverters_4 check_inverters_5 check_inverters_6 check_inverters_8

default: $(EXES)


MODDIR := Geometry
MODS := geometry_SAP geometry_init geometry_mpi_eo test_geometry_mpi communications communications_flt process_init
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := logger read_input
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := boundary_conditions single_double_utils suN_utils time
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields ranz2
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := Dphi random_spinor_field random_su2 representation Dphi_flt
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := DDalphaAMG FGMRES GMRES SAP_precondition linear_algebra cg_mshift MINRES_mshift MINRES g5QMR_mshift eva jacobi cg_mshift_flt
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := g5QMR_trunc trunc_hairpin avr_plaquette
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

COMMON_OBJS = $(MOD_OBJS)


include $(MKDIR)/MkRules

.PHONY: default

