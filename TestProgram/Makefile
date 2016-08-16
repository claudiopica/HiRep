TOPDIR = ..
SUBDIRS := $(shell find * -maxdepth 0 -type d)

#Exclude these directories
EXCLUDEDIR := TestConfigurations TestConfigurations-test
SUBDIRS := $(filter-out $(EXCLUDEDIR), $(SUBDIRS))

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules

tests: 
	for i in $(SUBDIRS) ; \
	do \
	  $(MAKE) -C $$i ;\
	done

.PHONY: tests
