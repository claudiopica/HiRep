#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"

#include <stdio.h>

extern rhmc_par _update_par;

#define _print_avect(a) printf("(%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f)\n",(a).c1,(a).c2,(a).c3,(a).c4,(a).c5,(a).c6,(a).c7,(a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.re,(a).c1_2.re,(a).c1_3.re,(a).c2_1.re,(a).c2_2.re,(a).c2_3.re,(a).c3_1.re,(a).c3_2.re,(a).c3_3.re);printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.im,(a).c1_2.im,(a).c1_3.im,(a).c2_1.im,(a).c2_2.im,(a).c2_3.im,(a).c3_1.im,(a).c3_2.im,(a).c3_3.im)


void Force0(float dt, suNg_algebra_vector *force){
  static suNg s1,s2;
  static suNg_algebra_vector f;
  int i;
  for (i=0; i<4*VOLUME; ++i){
    int x, mu;
    index_to_coord(i,x,mu);
    staples(x,mu,&s1);
    _suNg_times_suNg_dagger(s2,*(u_gauge+i),s1);
    
    /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
    _fund_algebra_project(f,s2);
    
    _algebra_vector_mul_add_assign_g(force[i], dt*(-_update_par.beta/((double)(NG))), f);
  }
}

void Force(float dt, suNg_algebra_vector *force){
  Force0(dt, force);
  Force_rhmc_f(dt, force);
}
