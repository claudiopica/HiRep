TOPDIR = ..

SRCS = 

all: suN.h suN_types.h suN_repr_func.h

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules


suN.h suN_types.h: $(MKDIR)/MkFlags $(WRHEAD)
	echo -n "Writing suN headers..."
	$(WRHEAD) $(NG) $(REPR) &&\
		echo " done." || echo " failed."

suN_repr_func.h: $(MKDIR)/MkFlags $(WRREPR) TMPL/suN_repr_func.h.tmpl
	echo -n "Writing representation header..."
	$(WRREPR) $(NG) TMPL/suN_repr_func.h.tmpl > $@ &&\
		echo " done." || ( echo " failed." && false )


clean::
	rm -f *~ suN.h suN_types.h suN_repr_func.h HEADERS.DONE.h

.PHONY: all
