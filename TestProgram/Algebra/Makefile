TOPDIR = ../../

MAIN := check_algebra_1

SRCS = $(addsuffix .c,$(MAIN))

default: $(MAIN)

MKDIR=$(TOPDIR)/Make
include $(MKDIR)/MkRules

MODDIR := IO
MODS := logger
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

$(MAIN): %: %.o $(MOD_OBJS)
	echo -n Linking [$@]... 
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(MOD_OBJS) &&\
	echo " done." || echo " failed."


.PHONY: default

cleanall: clean
	rm -rf $(MAIN) *~
