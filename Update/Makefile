TOPDIR = ..
MKDIR = $(TOPDIR)/Make

SRCS = cabmar.c staples.c update_hb.c random_su2.c Dphi.c Dphi_flt.c representation.c leapfrog.c force0.c random_momenta.c local_action.c update_rhmc.c update_hmc.c rational_functions.c copy_gfield.c random_momenta.c random_spinor_field.c find_spec_limits.c force_hmc.c force_rhmc.c O2MN_multistep.c

include $(MKDIR)/MkRules

approx_data.db: remez/approx_*
	echo -n "Writing rational approximation database..."
	( cd remez && cat approx_* | ./mkappdata.pl >../$@ && cd .. \
		&& echo " done." ) || ( echo " failed." && false )
