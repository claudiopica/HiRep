TOPDIR = ..
SUBDIRS := $(shell find * -maxdepth 0 -type d)

#Exclude these directories
EXCLUDEDIR :=  Deflate RotatedSF Scattering
SUBDIRS := $(filter-out $(EXCLUDEDIR), $(SUBDIRS))

LIBS += $(TOPDIR)/LibHR/libhr.a

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules


tests: $(LIBS)
	touch $(TOPDIR)/LibHR/.lib.updated
	for i in $(SUBDIRS) ; \
	do \
	  $(MAKE) -C $$i ;\
	done
	rm -f touch $(TOPDIR)/LibHR/.lib.updated

runalltests: 
	for i in $(SUBDIRS) ; \
	do \
	  echo "-----------------------------------"; echo $$i ; \
	  $(MAKE) -C $$i runtests;\
	done
	echo ; echo "-----------------------------------"; echo "Summary"
	rm -f .test_failed
	for i in $(SUBDIRS) ; \
	do \
	  if [ -f $$i/.test_failed ];\
	  then \
	  	echo "$(_col_failed) $$i FAIL$(_col_norm)" ;\
		touch .test_failed ; \
	  else\
	  	echo "$(_col_done) $$i OK$(_col_norm)" ;\
	  fi ;\
	done
	if [ -f .test_failed ];\
	then exit 1 ; \
	else exit 0 ; \
	fi

.PHONY: tests
