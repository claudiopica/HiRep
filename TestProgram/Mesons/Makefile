TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

TESTS = check_triplets_1 check_spinmatrix check_triplets_2 
FAIL = check_triplets_3 check_triplets_4

check_triplets_3_OBJS = ../../Spectrum/meson_measurements.o ../../GaugeFix/gaugefix.o
check_triplets_4_OBJS = ../../Spectrum/meson_measurements.o ../../GaugeFix/gaugefix.o

LIBS += $(TOPDIR)/LibHR/libhr.a

include $(MKDIR)/MkRules
