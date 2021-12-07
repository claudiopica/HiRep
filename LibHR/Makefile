TOPDIR = ..
MKDIR = $(TOPDIR)/Make
SUBDIRS1 = $(shell find * -maxdepth 0 -type d)
SUBDIRS = $(filter-out Observables, $(SUBDIRS1))

LIBDIRS = $(SUBDIRS)
LIBNAME = ./libhr.a

include $(MKDIR)/MkRules

