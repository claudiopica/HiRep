TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

EXES = check_algebra_1 
#check_algebra_2 only compiles for NG=3 and REPR=AS
EXES += $(if $(findstring 3,$(NG)),$(if $(findstring REPR_ANTISYMMETRIC,$(REPR)),check_algebra_2,),)



LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules

