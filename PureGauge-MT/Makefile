TOPDIR = ..
MKDIR=$(TOPDIR)/Make

EXES = puregauge

default: $(EXES)

MODDIR := Geometry
MODS := process_init geometry_init geometry_mpi_eo test_geometry_mpi communications communications_flt
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := archive logger read_input endian
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := suN_utils suN_exp single_double_utils boundary_conditions time
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc spinor_field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette polyakov
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss random_suNg ranlxs ranlxd random_fields ran_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := staples update_hb random_su2 random_momenta update_puregauge local_action force0 copy_gfield O2MN_multistep cabmar representation
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Statistics
MODS := 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Inverters
MODS := #linear_algebra
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := PureGauge-MP
MODS := puregauge_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

puregauge_OBJS = $(MOD_OBJS)


include $(MKDIR)/MkRules

.PHONY: default

