TOPDIR = ..
MKDIR=$(TOPDIR)/Make

EXES = hmc

default: $(EXES)

MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo test_geometry_mpi communications communications_flt
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input endian
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := suN_utils suN_exp single_double_utils boundary_conditions time max_eig
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette polyakov z2semwall mesons # sf_action sf_pcac
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields ranz2
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := Dphi Dphi_flt staples update_hb random_su2 representation random_momenta random_spinor_field update_hmc local_action force0 copy_gfield force_hmc O2MN_multistep cabmar
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Statistics
MODS := 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := linear_algebra cg_mshift MINRES_mshift g5QMR_mshift MINRES eva jacobi HBiCGstab
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := HMC
MODS := hmc_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

hmc_OBJS = $(MOD_OBJS)

#rhmc: $(OBJS) $(MOD_OBJS)
#	echo -n "Linking [$@]..." 
#	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(MOD_OBJS) &&\
#	echo " done." || ( echo " failed." && false )


include $(MKDIR)/MkRules

.PHONY: default

