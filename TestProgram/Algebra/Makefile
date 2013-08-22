TOPDIR = ../..
MKDIR=$(TOPDIR)/Make

EXES = check_algebra_1 check_algebra_2

default: $(EXES)


MODDIR := Geometry
MODS := communications geometry_init 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := logger
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := boundary_conditions single_double_utils
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Random
MODS := gauss ranlxs ranlxd 
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Update
MODS := representation
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))


MODDIR := Observables
MODS := avr_plaquette
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))


COMMON_OBJS = $(MOD_OBJS)

include $(MKDIR)/MkRules

.PHONY: default

