TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

SRCS = $(wildcard *.c)
GPUFILES = linear_algebra.c global_sum_gpu.c
LIBNAME = ../libhr.a

#Exclude this source files
EXCLUDE = HBiCGstab.c BiCGstab_mshift.c HBiCGstab_mshift.c dirac_eva.c

include $(MKDIR)/MkRules

