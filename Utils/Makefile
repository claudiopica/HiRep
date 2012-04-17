TOPDIR = ..
MKDIR = $(TOPDIR)/Make

SRCS = suN_utils.c suN_exp.c boundary_conditions.c time.c single_double_utils.c HYP_smearing.c single_double_utils_gpu.c update_hb_gpu.c
GPUFILES = single_double_utils_gpu.c update_hb_gpu.c

GPUFILES = single_double_utils_gpu.c

include $(MKDIR)/MkRules


suN_exp.c: $(MKDIR)/MkFlags $(WRREPR) TMPL/suN_exp.c.tmpl
	echo -n "Writing representation code..."
	$(WRREPR) $(NG) TMPL/suN_exp.c.tmpl > $@ &&\
	  echo " done." || ( echo " failed." && false )

clean::
	rm -f suN_exp.c


