#include "global.h"
#include "observables.h"
#include "update.h"
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "representation.h"

extern rhmc_par _update_par;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      double *loc_action,
                      suNg_algebra_vector *momenta,
                      suNf_spinor **phi1,
                      suNf_spinor **phi2) {

	int i,j;
	double a,tmp;
	for(i=0; i<VOLUME; ++i){

		/* Momenta */
		a=0.;
		for (j=0;j<4;++j) {
			suNg_algebra_vector *cmom=gfield_ordering(momenta,i,j);
			_algebra_vector_sqnorm_g(tmp,*cmom); 
			a+=tmp; /* this must be positive */
		}
		a*=0.5*_FUND_NORM2;

		/* Gauge action */
		a -= (_update_par.beta/((double)NG))*local_plaq(i);

		/* Fermions */
		for (j=0;j<_update_par.n_pf;++j) {
			_spinor_prod_re_f(tmp,phi1[j][i],phi2[j][i]);
			a += tmp;
		}

		switch(type) {
			case NEW:
				*(loc_action+i)=a;
				break;
			case DELTA:
				a-=*(loc_action+i);
				*(loc_action+i)=a;
				break;
			default:
				error(1,1,"local_hmc_action","Invalid type");
		}

	}
   
}
