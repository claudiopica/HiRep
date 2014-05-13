TOPDIR = ..
MKDIR = $(TOPDIR)/Make

EXES = mk_modenumber mk_eigvals

mk_modenumber_OBJS = modenumber.o


LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules


# THIS HAS TO BE ADAPTED TO THE MACHINE
approx_for_modenumber: approx_for_modenumber.c
	gcc -o approx_for_modenumber approx_for_modenumber.c -Wall -std=c99 -lm -lgsl -lgslcblas -I/sw/include -L/sw/lib -m32

