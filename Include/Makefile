TOPDIR = ..
MKDIR = $(TOPDIR)/Make

#the default target cannot start with .
default: .autoheaders.h
 
.autoheaders.h: suN.h suN_types.h suN_repr_func.h
	echo done > $@

include $(MKDIR)/MkRules

suN.h suN_types.h: $(MKDIR)/MkFlags $(WRHEAD)
	echo -n "Writing suN headers..."
	$(WRHEAD) $(NG) $(REPR) $(WQUAT) &&\
		echo " done." || ( echo " failed." && false )

suN_repr_func.h: $(MKDIR)/MkFlags $(WRREPR) TMPL/suN_repr_func.h.tmpl
	echo -n "Writing representation header..."
	$(WRREPR) $(NG) TMPL/suN_repr_func.h.tmpl > $@ &&\
		echo " done." || ( echo " failed." && false )


clean::
	rm -f *~ suN.h suN_types.h suN_repr_func.h .autoheaders.h


.PHONY: default

