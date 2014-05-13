TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

SRCS = $(wildcard *.c) suN_exp.c
LIBNAME = ../libhr.a

#Exclude this source files
EXCLUDE += $(if $(findstring 2,$(NG)),,HYP_smearing.c)

include $(MKDIR)/MkRules

suN_exp.c: $(MKDIR)/MkFlags $(WRREPR) TMPL/suN_exp.c.tmpl
	echo -n "$(_col_write)Writing representation code$(_col_norm)..."
	$(WRREPR) $(NG) TMPL/suN_exp.c.tmpl > $@ &&\
	  echo " $(_cdone)" || ( echo " $(_cfailed)" && false )

clean::
	rm -f suN_exp.c

