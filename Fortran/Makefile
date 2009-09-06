FC = ifort -cxxlib -O2 -align -Zp16 -r8
#FC = gfortran -O2

TOPDIR = ..
MKDIR=$(TOPDIR)/Make

EXES = su4dbl

default: $(EXES)


MODDIR := Random
MODS := ranlxd
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Geometry
MODS := geometry_init geometry_mpi_eo communications
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := IO
MODS := logger endian archive
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Memory
MODS := amalloc field_alloc
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Error
MODS := error
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Converter
MODS := archive_eolexi archive_fortran archive_more_mpieo
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Utils
MODS := time
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

MODDIR := Observables
MODS := avr_plaquette
MOD_OBJS := $(MOD_OBJS) $(patsubst %, $(TOPDIR)/$(MODDIR)/%.o, $(MODS))

$(EXES): %: %.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(FC) -o $@ luesher.f $< $(MOD_OBJS) &&\
	echo " done." || echo " failed."


include $(MKDIR)/MkRules

.PHONY: default

