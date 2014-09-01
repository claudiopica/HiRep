TOPDIR = ..
MKDIR = $(TOPDIR)/Make

#the default target cannot start with .
default: .autoheaders.h
 
.autoheaders.h: suN.h suN_types.h suN_repr_func.h macro_opt.h
	echo done > $@


#this must be a stem rule otherwise it will be called twice in parallel make
su%.h su%_types.h: $(MKDIR)/MkFlags $(WRHEAD)
	echo -n "$(_col_write)Writing suN headers$(_col_norm)..."
	$(WRHEAD) $(NG) $(REPR) $(WQUAT) $(GAUGE_GROUP) &&\
		echo " $(_cdone)" || ( echo " $(_cfailed)" && false )

suN_repr_func.h: $(MKDIR)/MkFlags TMPL/suN_repr_func.h.tmpl
	echo -n "$(_col_write)Writing representation header$(_col_norm)..."
	$(WRREPR) $(NG) TMPL/suN_repr_func.h.tmpl > $@ &&\
		echo " $(_cdone)" || ( echo " (_cfailed)" && false )


clean::
	rm -f *~ suN.h suN_types.h suN_repr_func.h macro_opt.h .autoheaders.h

include $(MKDIR)/MkRules

#this must be after MkRules which defines WRREPR
suN_repr_func.h: $(WRREPR)

suN.h suN_types.h: $(WRHEAD)

macro_opt.h: $(MKDIR)/MkFlags
	perl -e '$$_="$(CPPFLAGS)"."\n"; s/\s+/\n/g; s/-I\s*.*?\n//g; s/-D/#define /g; s/=/ /g; print; print "\n";' > $@

