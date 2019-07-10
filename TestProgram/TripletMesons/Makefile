TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

TESTS = check_spinmatrix check_triplets_1 check_triplets_2 check_triplets_3
#TESTS = check_triplets_3

check_triplets_3_OBJS = ../../Spectrum/meson_measurements.o ../../GaugeFix/gaugefix.o

LIBS += $(TOPDIR)/LibHR/libhr.a

include $(MKDIR)/MkRules
