TOPDIR = ..
MKDIR=$(TOPDIR)/Make

SUBDIRS = Algebra Geometry PureGauge SpinorField DiracOperator Inverters Update Random TripletMesons

# MODDIR := Geometry
# MODS := geometry_init  geometry_mpi_eo test_geometry_mpi geometry_mpi
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := IO
# MODS := archive logger
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Utils
# MODS := suN_utils suN_exp single_double_utils boundary_conditions
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Memory
# MODS := amalloc #field_alloc spinor_field_alloc
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Observables
# MODS := avr_plaquette mesons quark_propagator
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Random
# MODS := gauss random_suNg ranlxs ranlxd random_fields
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Update
# MODS := cabmar staples update_hb random_su2 Dphi Dphi_old representation random_momenta random_spinor_field find_spec_limits
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Error
# MODS := error
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Statistics
# MODS := base_stat jacknife
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

# MODDIR := Inverters
# MODS := linear_algebra cg_mshift BiCGstab_mshift HBiCGstab_mshift MINRES_mshift g5QMR_mshift MINRES eva jacobi dirac_eva
# MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))


include $(MKDIR)/MkRules

.PHONY: default

