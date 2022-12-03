TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

SRCS = $(wildcard *.c)
LIBNAME = ../libhr.a

#Exclude this source files
EXCLUDE = 
EXCLUDE += measure_mesons.c measure_baryons.c z2semwall_new.c measure_mesons_ff.c z2semwall.c measure_ff.a
#EXCLUDE += loop_tools.c mesons.c polyakov.c
EXCLUDE += measure_scattering.c
#EXCLUDE += glueballs.c glueballs_op.c
EXCLUDE += measure_renormalization.c
#EXCLUDE += wilsonloops.c
EXCLUDE += trunc_hairpin.c

include $(MKDIR)/MkRules

