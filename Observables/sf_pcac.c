#include "global.h"
#include "communications.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "error.h"
#include "geometry.h"
#include "spinor_field.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"

static double hmass;
static void H_sf(spinor_field *out, spinor_field *in){
  g5Dphi(hmass,out,in);
}

void SF_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc) {
  static MINRES_par MINRESpar;
  int cgiter;
  hmass = mass;
  
  MINRESpar.err2 = acc;
  MINRESpar.max_iter = 0;
  cgiter=0;
  cgiter+=MINRES(&MINRESpar, &H_sf, in, out, 0);
	  lprintf("PROPAGATOR",10,"MINRES MVM = %d",cgiter);
}



/*
f1 can be computed inserting the sources in the correlator in the lower or higher border. The two procedures give identical results. We keep only the computation with sources on the lower border.
*/
double SF_PCAC_wall_mass(double mass)
{

#ifdef BASIC_SF

  int i,j,ix0,ix1,ix2,ix3;
  double f_P[GLB_T], f_A[GLB_T], f_Pt[GLB_T], f_At[GLB_T], f_1=0, temp;
  double acc = 1.e-10;
  spinor_field *prop;	
  spinor_field *source;
  suNf_spinor* sptr;
  suNf_spinor* stmp;
  stmp=malloc(2*sizeof(suNf_spinor));
  suNf_spinor* sbord;
  sbord=malloc(4*NF*sizeof(suNf_spinor));
  suNf *uptr;
  prop=alloc_spinor_field_f(4*NF,&glattice);
  source=alloc_spinor_field_f(4*NF,&glattice);    

  for(ix0=0;ix0<GLB_T;ix0++)
    {
      f_A[ix0]=0;
      f_At[ix0]=0;
      f_P[ix0]=0;
      f_Pt[ix0]=0;
    }


  
  /*U' and P+ on source (actually P- since there is a g5 that needs to be commuted through)*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    if(COORD[0]==0){
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s%4].c[s/4].re=1.;
      stmp[0].c[s%4].c[s/4].im=0.;
      _spinor_g5_assign_f(stmp[0]);

      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(idn(i,0),0);
      	for(j=0;j<4;j++) {
      	  _suNf_inverse_multiply(stmp[1].c[j],*uptr,stmp[0].c[j]);
    	  }
      	sptr = _FIELD_AT(&source[s],i);
      	_spinor_pplus_f(*sptr,stmp[1]);
      }
    }
  }
  
  /*get propagator to all points*/
  for(int s=0; s<4*NF; s++){
    spinor_field_zero_f(&prop[s]);
    SF_quark_propagator(&source[s], mass, &prop[s], acc); 
  }
  
  /*get time averaged correlators for each timeslice*/
  /*f_P*/
  for(int s=0;s<4*NF;s++){
    for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      i=ipt(ix0,ix1,ix2,ix3);
      sptr = _FIELD_AT(&prop[s],i);
      /*f_P*/
      _spinor_prod_re_f(temp,*sptr,*sptr);
      f_P[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=temp;
      /*f_A*/
      /*gamma_0*/
      stmp[0].c[0]=sptr->c[2];
      stmp[0].c[1]=sptr->c[3];
      stmp[0].c[2]=sptr->c[0];
      stmp[0].c[3]=sptr->c[1];
      _spinor_prod_re_f(temp,*sptr,stmp[0]);
      f_A[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=temp;
    }
  }
  
/* f_P = prop^dag prop */
/* f_A = prop^dag gamma0 prop */

  global_sum((double*)f_P,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Ppost%d = %.15f\n",ix0,f_P[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	

  global_sum((double*)f_A,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Apost%d = %.15f\n",ix0,f_A[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
  

  /*f_1 - NEED TO DO EACH color/dirac component separately, then combine at the end*/
  /*U' and P+ on prop at T-2 (actually P- since there is a g5 that needs to be commuted through)*/

  for(int s=0;s<4*NF;s++){
    _spinor_zero_f(sbord[s]);
    if(COORD[0]==NP_T-1){
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(T-2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(i,0);
      	sptr = _FIELD_AT(&prop[s],i);
      	for(j=0;j<4;j++) {
      	  _suNf_inverse_multiply(stmp[0].c[j],*uptr,sptr->c[j]);
      	}
      	_spinor_pminus_f(stmp[1],stmp[0]);
      	_spinor_add_assign_f(sbord[s],stmp[1]);
      }
    }
  }
  
  global_sum((double*)sbord,sizeof(suNf_spinor)/sizeof(double)*4*NF);

  if(PID==0){
    f_1=0;
    for(int s=0;s<4*NF;s++){
      _spinor_prod_re_f(temp,sbord[s],sbord[s]);
      f_1+=temp;
    }
  }
  lprintf("PC_wall_AC",0,"f1 = %.15f\n",f_1/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z)));
  lprintf("PC_wall_AC",0,"ZP_pos = %.15f\n",(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));
  
  for (ix0=2;ix0<GLB_T-3;ix0++)
    lprintf("PC_wall_AC",0,"PCACpost%d = %f\n",ix0,(double)(f_A[(int)(ix0)+1] - f_A[(int)(ix0)-1])/(4*f_P[(int)(ix0)]));
  
  /*Create wall source with g5 factor at t=T-2*/
  /*U and P- on source (again actually use P+ to account for commuting with g5 in source)*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    if(COORD[0]==NP_T-1){
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s%4].c[s/4].re=1.;
      stmp[0].c[s%4].c[s/4].im=0.;
      _spinor_g5_assign_f(stmp[0]);
    
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++){
      	i=ipt(T-2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(i,0);
      	for(j=0;j<4;j++) {
      	  _suNf_multiply(stmp[1].c[j],*uptr,stmp[0].c[j]);
    	  }
      	sptr = _FIELD_AT(&source[s],i);
      	_spinor_pminus_f(*sptr,stmp[1]);
      }
    }
  }

  /*get propagator to all points*/
  for(int s=0; s<4*NF; s++){
    spinor_field_zero_f(&prop[s]);
    SF_quark_propagator(&source[s], mass, &prop[s], acc); 
  }

  /*get time averaged correlators for each timeslice (going back from T in time)*/
  for(int s=0;s<4*NF;s++){
    for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) for(ix0=0;ix0<T;ix0++){
      i=ipt(ix0,ix1,ix2,ix3);
      sptr = _FIELD_AT(&prop[s],i);
      /*f_P*/
      _spinor_prod_re_f(temp,*sptr,*sptr);
      f_Pt[((GLB_T-1)-(COORD[0]*T+ix0))%GLB_T]+=temp;
      /*f_A*/
      /*gamma_0*/
      _vector_mul_f(stmp[0].c[0],-1,sptr->c[2]);
      _vector_mul_f(stmp[0].c[1],-1,sptr->c[3]);
      _vector_mul_f(stmp[0].c[2],-1,sptr->c[0]);
      _vector_mul_f(stmp[0].c[3],-1,sptr->c[1]);
      _spinor_prod_re_f(temp,*sptr,stmp[0]);
      f_At[((GLB_T-1)-(COORD[0]*T+ix0))%GLB_T]+=temp;
    }
  }
  
  global_sum((double*)f_Pt,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Pnegt%d = %.15f\n",ix0,f_Pt[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	

  global_sum((double*)f_At,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Anegt%d = %.15f\n",ix0,f_At[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	

  lprintf("PC_wall_AC",0,"ZP_neg = %.15f\n",(sqrt(f_1)/(f_Pt[(int)(GLB_X/2)])));
  
  lprintf("PC_wall_AC",0,"Z_P = %.15f\n",0.5*(sqrt(f_1)/(f_Pt[(int)(GLB_X/2)]))+0.5*(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));
  
  for (ix0=2;ix0<GLB_T-3;ix0++)
    lprintf("PC_wall_AC",0,"PCACnegt%d = %f\n",ix0,(double)(f_At[(int)(ix0)+1] - f_At[(int)(ix0)-1])/(4*f_Pt[(int)(ix0)]));
  
  free_spinor_field(source);
  free_spinor_field(prop);
  

  free(stmp);
  free(sbord);

  return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]);
  
  
#elif defined(ROTATED_SF)

  int i,j,ix0,ix1,ix2,ix3;
  double f_P[GLB_T], f_A[GLB_T], f_Pt[GLB_T], f_At[GLB_T], f_1=0, temp;
  double acc = 1.e-10;
  spinor_field* prop[2];	
  spinor_field *source;
  suNf_spinor* sptr[2];
  suNf_spinor* stmp;
  stmp=malloc(2*sizeof(suNf_spinor));
  suNf_spinor* sbord[2];
  sbord[0]=malloc(8*NF*sizeof(suNf_spinor));
  sbord[1]=sbord[0]+4*NF;
  suNf *uptr;
  prop[0]=alloc_spinor_field_f(4*NF,&glattice);
  prop[1]=alloc_spinor_field_f(4*NF,&glattice);
  source=alloc_spinor_field_f(4*NF,&glattice);    

  for(ix0=0;ix0<GLB_T;ix0++)
    {
      f_A[ix0]=0;
      f_At[ix0]=0;
      f_P[ix0]=0;
      f_Pt[ix0]=0;
    }

/*
[eta_(a0,alpha0)](t,x)_{a,alpha} = \delta_{a,a0} \delta_{alpha,alpha0) \delta_{t,1}
source_s = gamma_5 Q- U_0(1,y)^\dag \eta_s = Q+ U_0(1,y)^\dag gamma_5 \eta_s
Q+ = (1 + i gamma_0 gamma_5)/2
*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    if(COORD[0]==0){
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s%4].c[s/4].re=1.;
      stmp[0].c[s%4].c[s/4].im=0.;
      _spinor_g5_assign_f(stmp[0]);

      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(idn(i,0),0);
      	for(j=0;j<4;j++) {
      	  _suNf_inverse_multiply(stmp[1].c[j],*uptr,stmp[0].c[j]);
    	  }
      	sptr[0] = _FIELD_AT(&source[s],i);
      	*sptr[0]=stmp[1];
      	_vector_i_sub_assign_f(sptr[0]->c[0],stmp[1].c[2]);
      	_vector_i_sub_assign_f(sptr[0]->c[1],stmp[1].c[3]);
      	_vector_i_add_assign_f(sptr[0]->c[2],stmp[1].c[0]);
      	_vector_i_add_assign_f(sptr[0]->c[3],stmp[1].c[1]);
        _spinor_mul_f(*sptr[0],.5,*sptr[0]);
      }
    }
  }
  
  /*
  prop[0][s] = H^{-1} source[s]
  prop[1][s] = H^{-1} gamma_5 source[s]
  */
  for(int s=0; s<4*NF; s++){
    spinor_field_zero_f(&prop[0][s]);
    SF_quark_propagator(&source[s], mass, &prop[0][s], acc);
    spinor_field_zero_f(&prop[1][s]);
    spinor_field_g5_assign_f(&source[s]);
    SF_quark_propagator(&source[s], mass, &prop[1][s], acc);
    spinor_field_g5_assign_f(&source[s]);
  }
  
  /*
  f_P(x) \propto \sum_s prop[1][s]^dag(x) gamma_5 prop[0][s](x)
  f_A(x) \propto \sum_s prop[1][s]^dag(x) gamma_0 prop[0][s](x)
  CONTROLLARE IN PARTICOLARE SE SERVE LA PARTE REALE O IMMAGINARIA
  */
  for(int s=0;s<4*NF;s++){
    for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      i=ipt(ix0,ix1,ix2,ix3);
      sptr[0] = _FIELD_AT(&prop[0][s],i);
      sptr[1] = _FIELD_AT(&prop[1][s],i);
      /*f_P*/
      _spinor_g5_prod_re_f(temp,*sptr[1],*sptr[0]);
      f_P[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=temp;
      /*f_A*/
      /*gamma_0*/
      stmp[0].c[0]=sptr[0]->c[2];
      stmp[0].c[1]=sptr[0]->c[3];
      stmp[0].c[2]=sptr[0]->c[0];
      stmp[0].c[3]=sptr[0]->c[1];
      _spinor_prod_re_f(temp,*sptr[1],stmp[0]);
      f_A[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=temp;
    }
  }
  
  global_sum((double*)f_P,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Ppost%d = %.15f\n",ix0,f_P[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	

  global_sum((double*)f_A,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Apost%d = %.15f\n",ix0,f_A[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
  

  /*
  f_1 \propto  \sum_s sum_u prop[1][s]^dag(T-2,u) U_0(T-2,u) gamma_5 Q+ \sum_v U_0(T-2,v)^dag prop[0][s](T-2,v)
  sbord[0][s] = \sum_v gamma_5 Q+ U_0(T-2,v)^dag prop[0][s](T-2,v)
  sbord[1][s] = \sum_v U_0(T-2,v)^dag prop[1][s](T-2,v)
  */

  for(int s=0;s<4*NF;s++){
    _spinor_zero_f(sbord[0][s]);
    _spinor_zero_f(sbord[1][s]);
    if(COORD[0]==NP_T-1){
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(T-2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(i,0);
      	sptr[0] = _FIELD_AT(&prop[0][s],i);
      	sptr[1] = _FIELD_AT(&prop[1][s],i);
      	for(j=0;j<4;j++) {
      	  _suNf_inverse_multiply(stmp[0].c[j],*uptr,sptr[0]->c[j]);
      	  _suNf_inverse_multiply(stmp[1].c[j],*uptr,sptr[1]->c[j]);
      	}

      	_spinor_add_assign_f(sbord[1][s],stmp[1]);

      	stmp[1].c[0]=stmp[0].c[0];
      	stmp[1].c[1]=stmp[0].c[1];
      	_vector_minus_f(stmp[1].c[2],stmp[0].c[2]);
      	_vector_minus_f(stmp[1].c[3],stmp[0].c[3]);
      	_vector_i_sub_assign_f(sptr[1]->c[0],stmp[0].c[2]);
      	_vector_i_sub_assign_f(sptr[1]->c[1],stmp[0].c[3]);
      	_vector_i_sub_assign_f(sptr[1]->c[2],stmp[0].c[0]);
      	_vector_i_sub_assign_f(sptr[1]->c[3],stmp[0].c[1]);
        _spinor_mul_f(*sptr[1],.5,*sptr[1]);
      	_spinor_add_assign_f(sbord[0][s],stmp[1]);
      }
    }
  }
  
  global_sum((double*)sbord[0],sizeof(suNf_spinor)/sizeof(double)*8*NF);

  if(PID==0){
    f_1=0;
    for(int s=0;s<4*NF;s++){
      _spinor_prod_re_f(temp,sbord[1][s],sbord[0][s]);
      f_1+=temp;
    }
  }
  lprintf("PC_wall_AC",0,"f1 = %.15f\n",f_1/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z)));
  lprintf("PC_wall_AC",0,"ZP_pos = %.15f\n",(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));
  
  for (ix0=2;ix0<GLB_T-3;ix0++)
    lprintf("PC_wall_AC",0,"PCACpost%d = %f\n",ix0,(double)(f_A[(int)(ix0)+1] - f_A[(int)(ix0)-1])/(4*f_P[(int)(ix0)]));
  



/*
[eta_(a0,alpha0)](t,x)_{a,alpha} = \delta_{a,a0} \delta_{alpha,alpha0) \delta_{t,T-1}
source_s = gamma_5 Q+ U_0(T-2,y) \eta_s = Q- U_0(T-2,y) gamma_5 \eta_s
Q- = (1 - i gamma_0 gamma_5)/2
*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    if(COORD[0]==NP_T-1){
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s%4].c[s/4].re=1.;
      stmp[0].c[s%4].c[s/4].im=0.;
      _spinor_g5_assign_f(stmp[0]);

      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(T-2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(i,0);
      	for(j=0;j<4;j++) {
      	  _suNf_multiply(stmp[1].c[j],*uptr,stmp[0].c[j]);
    	  }
      	sptr[0] = _FIELD_AT(&source[s],i);
      	*sptr[0]=stmp[1];
      	_vector_i_add_assign_f(sptr[0]->c[0],stmp[1].c[2]);
      	_vector_i_add_assign_f(sptr[0]->c[1],stmp[1].c[3]);
      	_vector_i_sub_assign_f(sptr[0]->c[2],stmp[1].c[0]);
      	_vector_i_sub_assign_f(sptr[0]->c[3],stmp[1].c[1]);
        _spinor_mul_f(*sptr[0],.5,*sptr[0]);
      }
    }
  }

  
  /*
  prop[0][s] = H^{-1} source[s]
  prop[1][s] = H^{-1} gamma_5 source[s]
  */
  for(int s=0; s<4*NF; s++){
    spinor_field_zero_f(&prop[0][s]);
    SF_quark_propagator(&source[s], mass, &prop[0][s], acc);
    spinor_field_zero_f(&prop[1][s]);
    spinor_field_g5_assign_f(&source[s]);
    SF_quark_propagator(&source[s], mass, &prop[1][s], acc);
    spinor_field_g5_assign_f(&source[s]);
  }
  
  /*
  f_P(x) \propto \sum_s prop[1][s]^dag(x) gamma_5 prop[0][s](x)
  f_A(x) \propto \sum_s prop[1][s]^dag(x) gamma_0 prop[0][s](x)
  CONTROLLARE IN PARTICOLARE SE SERVE LA PARTE REALE O IMMAGINARIA
  */
  for(int s=0;s<4*NF;s++){
    for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      i=ipt(ix0,ix1,ix2,ix3);
      sptr[0] = _FIELD_AT(&prop[0][s],i);
      sptr[1] = _FIELD_AT(&prop[1][s],i);
      /*f_P*/
      _spinor_g5_prod_re_f(temp,*sptr[1],*sptr[0]);
      f_P[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=temp;
      /*f_A*/
      /*gamma_0*/
      stmp[0].c[0]=sptr[0]->c[2];
      stmp[0].c[1]=sptr[0]->c[3];
      stmp[0].c[2]=sptr[0]->c[0];
      stmp[0].c[3]=sptr[0]->c[1];
      _spinor_prod_re_f(temp,*sptr[1],stmp[0]);
      f_A[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=temp;
    }
  }
   
  global_sum((double*)f_Pt,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Pnegt%d = %.15f\n",ix0,f_Pt[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	

  global_sum((double*)f_At,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Anegt%d = %.15f\n",ix0,f_At[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	

  lprintf("PC_wall_AC",0,"ZP_neg = %.15f\n",(sqrt(f_1)/(f_Pt[(int)(GLB_X/2)])));
  
  lprintf("PC_wall_AC",0,"Z_P = %.15f\n",0.5*(sqrt(f_1)/(f_Pt[(int)(GLB_X/2)]))+0.5*(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));
  
  for (ix0=2;ix0<GLB_T-3;ix0++)
    lprintf("PC_wall_AC",0,"PCACnegt%d = %f\n",ix0,(double)(f_At[(int)(ix0)+1] - f_At[(int)(ix0)-1])/(4*f_Pt[(int)(ix0)]));
  
 
  return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]);
  
  free(stmp);
  free(sbord[0]);
  free_spinor_field(prop[0]);
  free_spinor_field(prop[1]);
  free_spinor_field(source);
 
#else

  return 0.;

#endif
  
}
