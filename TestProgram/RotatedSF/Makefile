TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

EXES = check_rotatedsf_1 check_rotatedsf_2 check_rotatedsf_3 check_rotatedsf_4

EXES += $(if $(findstring 3,$(NG)),$(if $(findstring REPR_FUNDAMENTAL,$(REPR)),check_rotatedsf_5 check_rotatedsf_6 check_rotatedsf_7 check_rotatedsf_8 check_rotatedsf_9 check_rotatedsf_10,),)

LIBS += $(TOPDIR)/LibHR/libhr.a

include $(MKDIR)/MkRules

