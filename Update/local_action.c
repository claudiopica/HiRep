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
 * | momenta |^2 = -6c1^2 - 6c2^2 + 6c1c2 - 2c3^2 - ... - 2c8^2
 */
void local_hmc_action(local_action_type type,
                      float *loc_action,
                      suNg_algebra_vector *momenta,
                      suNf_spinor *phi1,
                      suNf_spinor *phi2) {

  int i,j;
   double a;
   for(i=0; i<VOLUME; ++i){

     /* Momenta */
     a=0.;
     for (j=0;j<4;++j) {
       suNg_algebra_vector *cmom=gfield_ordering(momenta,i,j);
/*  
     a += ((cmom->c1)*(cmom->c1))+((cmom->c2)*(cmom->c2))
	 +((cmom->c3)*(cmom->c3))+((cmom->c4)*(cmom->c4))
	 +((cmom->c5)*(cmom->c5))+((cmom->c6)*(cmom->c6))
	 +((cmom->c7)*(cmom->c7))+((cmom->c8)*(cmom->c8));
*/
       a+=_algebra_vector_sqnorm_g(*cmom);
     }

     a*=0.5*_FUND_NORM2;
     /* a/=0.5*_FUND_NORM2; */
     
     /* Gauge action */
     a -= (_update_par.beta/((double)NG))*local_plaq(i);

     /* Fermions */
     a += _spinor_prod_re_f(phi1[i],phi2[i]);
     
     switch(type) {
         case NEW:
            *(loc_action+i)=a;
            break;
         case DELTA:
            *(loc_action+i)-=a;
            break;
         default:
            error(1,1,"local_hmc_action","Invalid type");
      }
      
   }
   
}
