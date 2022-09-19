TOPDIR = .
SUBDIRS := $(shell find * -maxdepth 0 -type d)

#Exclude these directories
EXCLUDEDIR := Analysis Doc Fortran HiRep.xcodeproj Run
SUBDIRS := $(filter-out $(EXCLUDEDIR), $(SUBDIRS))

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules

