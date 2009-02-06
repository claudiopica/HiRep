TOPDIR = ..

SRCS = suN_utils.c suN_exp.c boundary_conditions.c time.c single_double_utils.c

# these are not yet modified 
# single_double_utils.c

MKDIR = $(TOPDIR)/Make
include $(MKDIR)/MkRules

suN_exp.c: $(MKDIR)/MkFlags $(WRREPR) TMPL/suN_exp.c.tmpl
	echo -n "Writing representation code..."
	$(WRREPR) $(NG) TMPL/suN_exp.c.tmpl > $@ &&\
	  echo " done." || echo " failed."

clean::
	rm -f suN_exp.c
