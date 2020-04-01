TOPDIR = ..
MKDIR = $(TOPDIR)/Make

#EXES =  random_cnfg random_spinor mk_mesons_with_z2semwall mk_mesons_with_z2semwall_new  measure_spectrum measure_formfactor #mk_sfcoupling #trunc_mesons mk_mesons
EXES=measure_spectrum mk_sfcorrelators

#COMMON_OBJS 

measure_spectrum_OBJS = meson_measurements.o baryon_measurements.o meson_measurements_ff.o
#measure_formfactor_OBJS = meson_measurements.o 


LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules

