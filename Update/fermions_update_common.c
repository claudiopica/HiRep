#include "global.h"
#include "suN.h"
#include "update.h"
#include "memory.h"
#include "rational_functions.h"
#include "dirac.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"

/* State quantities for HMC */
suNg_algebra_vector *momenta=NULL;
suNf_spinor **pf=NULL;
rhmc_par _update_par;
rational_app r_S;  /* used for computing the action S in the metropolis test */
rational_app r_MD; /* used in the action MD evolution */
rational_app r_HB;  /* used in pseudofermions heatbath */
double minev, maxev; /* min and max eigenvalue of H^2 */
/* END of State */

static int init=1;

/* this is the basic operator used in the update */
static suNf_spinor *h2tmp=NULL;
#ifdef UPDATE_EO
void H2(suNf_spinor *out, suNf_spinor *in){
  g5Dphi_eopre(_update_par.mass, h2tmp, in);
  g5Dphi_eopre(_update_par.mass, out, h2tmp);
}
void H(suNf_spinor *out, suNf_spinor *in){
  g5Dphi_eopre(_update_par.mass, h2tmp, in);
	spinor_field_copy_f(out,h2tmp); /* in this way out can be == in */
}
#else
void H2(suNf_spinor *out, suNf_spinor *in){
  g5Dphi(_update_par.mass, h2tmp, in);
  g5Dphi(_update_par.mass, out, h2tmp);
}
void H(suNf_spinor *out, suNf_spinor *in){
  g5Dphi(_update_par.mass, h2tmp, in);
	spinor_field_copy_f(out,h2tmp); /* in this way out can be == in */
}
#endif

static void cleanup() {
	if (h2tmp!=NULL) {
		free_field(h2tmp);
		h2tmp=NULL;
	}
}

void init_fermions_common() {
	unsigned int slen;
	if (init) {
		get_spinor_len(&slen);
#ifdef UPDATE_EO
		set_spinor_len(VOLUME/2);
#else
		set_spinor_len(VOLUME);
#endif
		h2tmp=alloc_spinor_field_f(1);
		set_spinor_len(slen);
		atexit(&cleanup);
		init=0;
	}
}

