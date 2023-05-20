TOPDIR = .
SUBDIRS := $(shell find * -maxdepth 0 -type d)

#Exclude these directories
EXCLUDEDIR := Analysis Doc Fortran HiRep.xcodeproj Run Benchmarks build example_omp
SUBDIRS := $(filter-out $(EXCLUDEDIR), $(SUBDIRS))

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules

