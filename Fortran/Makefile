FC = ifort 
FOPT = -cxxlib -O2 -align -Zp16 -r8
#Fc = gfortran -O2
FSOURCES=su4dbl.f luesher.f
FOBJECTS=$(FSOURCES:.f=.o)
EXECUTABLE=su4dbl 

TOPDIR = ..
MKDIR=$(TOPDIR)/Make

default: all

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

all: $(FSOURCES) $(EXECUTABLE)

tidy: 
	rm $(EXECUTABLE) *.o

%.o: %.f 
	echo -n Compiling [$@]... 
	$(FC) $(FOPT) -c $< -o $@ &&\
	echo " done." || echo " failed."

$(EXECUTABLE): $(FOBJECTS) $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(FC) $(FOBJECTS) $(MOD_OBJS) -o su4dbl &&\
	echo " done." || echo " failed."

include $(MKDIR)/MkRules

.PHONY: default

