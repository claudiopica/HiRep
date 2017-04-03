TOPDIR = ../..
MKDIR = $(TOPDIR)/Make

SRCS = $(wildcard *.c)
LIBNAME = ../libhr.a

#Exclude these source files
EXCLUDE = 

include $(MKDIR)/MkRules

approx_data.db: remez/approx_*
	echo -n "$(_col_write)Writing rational approximation database$(_col_norm)..."
	( cd remez && cat approx_* | ./mkappdata.pl >../$@ && cd .. \
		&& echo " $(_cdone)" ) || ( echo " $(_cfailed)" && false )

