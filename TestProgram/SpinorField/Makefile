TOPDIR = ../..
MKDIR=$(TOPDIR)/Make

EXES = check_spinorfield_gpu_1 global_sum_gpu global_sum_gpu_flt #speed_flt_gpu speed_dbl_gpu check_spinorfield_gpu_1_flt speed_spinor_gpu_flt

default: $(EXES)


MODDIR := Geometry
MODS := geometry_init  geometry_mpi_eo test_geometry_mpi communications communications_flt
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := logger read_input
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := timer_gpu boundary_conditions # suN_utils suN_exp single_double_utils 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc field_pack_gpu field_copy_gpu 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette #mesons quark_propagator
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields ranz2
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := random_su2 representation random_spinor_field random_spinor_field_gpu # cabmar staples update_hb  Dphi Dphi_old  random_momenta  find_spec_limits
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error error_gpu
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Statistics
# MODS := base_stat jacknife
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra # cg_mshift BiCGstab_mshift HBiCGstab_mshift MINRES_mshift g5QMR_mshift MINRES eva jacobi dirac_eva
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

COMMON_OBJS = $(MOD_OBJS)

include $(MKDIR)/MkRules

.PHONY: default

