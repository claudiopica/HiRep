TOPDIR = ..
MKDIR=$(TOPDIR)/Make

EXES = gaugefix_measure gaugefix_unittest

default: $(EXES)

gaugefix_measure_OBJS = gaugefix.o
gaugefix_unittest_OBJS = gaugefix.o


LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules

.PHONY: default

