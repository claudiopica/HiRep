TOPDIR = .


SUBDIRS = Include Geometry IO Memory Observables Random Statistics Update Utils Error Inverters Make Converter TestProgram Disconnected Spectrum GaugeFix Scattering

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules

