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


#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */



static double hmass;
#ifdef BASIC_SF
static void H_sf(spinor_field *out, spinor_field *in){
  g5Dphi(hmass,out,in);
}
#else
static void H2_sf(spinor_field *out, spinor_field *in){
  g5Dphi_sq(hmass,out,in);
}
#endif

/*
this computes Hu^{-1}
*/
void SF_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc) {
#ifdef BASIC_SF

  static MINRES_par MINRESpar;
  int cgiter;
  hmass = mass;
  
  MINRESpar.err2 = acc;
  MINRESpar.max_iter = 0;
  cgiter=MINRES(&MINRESpar, &H_sf, in, out, 0);
	lprintf("PROPAGATOR",10,"MINRES MVM = %d",cgiter);
#else
	  
  int cgiter;
  static mshift_par inv_par;
  static spinor_field *chi=NULL;
  if(chi==NULL) chi=alloc_spinor_field_f(1,&glattice);
  hmass = mass;
                                        
  inv_par.n = 1;
  inv_par.shift = malloc(sizeof(double));
  inv_par.shift[0] = 0.;
  inv_par.err2= acc;
  inv_par.max_iter=0;
  /*this change of sign is equivalent to evaluate H^\dagger*/
  _update_par.SF_sign = -_update_par.SF_sign;
   g5Dphi(mass, chi, in);
  _update_par.SF_sign = -_update_par.SF_sign;

  cgiter=cg_mshift(&inv_par, &H2_sf, chi, out);
	lprintf("PROPAGATOR",10,"cg_mshift MVM = %d",cgiter);
	free(inv_par.shift);

#endif
}



/*
f1 can be computed inserting the sources in the correlator in the lower or higher border. The two procedures give identical results. We keep only the computation with sources on the lower border.
*/


void SF_PCAC_wall_mass(double mass, double acc )
{

#ifdef BASIC_SF

  int i,j,ix0,ix1,ix2,ix3;
  double f_P[GLB_T], f_A[GLB_T], f_Pt[GLB_T], f_At[GLB_T], f_1=0, temp;
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
      	_spinor_pminus_f(*sptr,stmp[1]);
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
    lprintf("PC_wall_AC",10,"f_Ppost%d = %.10e\n",ix0,f_P[ix0]/(double)GLB_VOL3);	

  global_sum((double*)f_A,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Apost%d = %.10e\n",ix0,f_A[ix0]/(double)GLB_VOL3);	
  

  /*f_1 - NEED TO DO EACH color/dirac component separately, then combine at the end*/
  /*U' and P- on prop at T-2 (actually P+ since there is a g5 that needs to be commuted through)*/

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
      	_spinor_pplus_f(stmp[1],stmp[0]);
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
  lprintf("PC_wall_AC",0,"f1 = %.10e\n",f_1/(double)GLB_VOL3/(double)GLB_VOL3);
  lprintf("PC_wall_AC",0,"ZP_pos = %.10e\n",(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));
  
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
      	_spinor_pplus_f(*sptr,stmp[1]);
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
    lprintf("PC_wall_AC",10,"f_Pnegt%d = %.10e\n",ix0,f_Pt[ix0]/(double)GLB_VOL3);	

  global_sum((double*)f_At,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Anegt%d = %.10e\n",ix0,f_At[ix0]/(double)GLB_VOL3);	

  lprintf("PC_wall_AC",0,"ZP_neg = %.10e\n",(sqrt(f_1)/(f_Pt[(int)(GLB_X/2)])));
  
  lprintf("PC_wall_AC",0,"Z_P = %.10e\n",0.5*(sqrt(f_1)/(f_Pt[(int)(GLB_X/2)]))+0.5*(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));
  
  for (ix0=2;ix0<GLB_T-3;ix0++)
    lprintf("PC_wall_AC",0,"PCACnegt%d = %f\n",ix0,(double)(f_At[(int)(ix0)+1] - f_At[(int)(ix0)-1])/(4*f_Pt[(int)(ix0)]));
  



  free_spinor_field(source);
  free_spinor_field(prop);
  

  free(stmp);
  free(sbord);

/*   return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]); */
  
  
#elif defined(ROTATED_SF)

/*

g1ud+ = 1/2 sum tr U0(1,z) Ddd^(-1)(2,z;T-2,z') U0(T-2,z') g5 Q+ U0(T-2,y')^dag Duu^(-1)(T-2,y',2,y) U0(1,y)^dag g5 Q+
      = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') g5 U0(T-2,z') g5 Q+ U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) g5 U0(1,y)^dag g5 Q+
      = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') U0(T-2,z') Q+  U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) U0(1,y)^dag Q+    
      = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') U0(T-2,z') Q+  U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) U0(1,y)^dag Q+



g1du+ = 1/2 sum tr csi^dag Q-  U0(1,z) Huu^(-1)(2,z;T-2,z') U0(T-2,z') Q-  U0(T-2,y')^dag Hdd^(-1)(T-2,y',2,y) U0(1,y)^dag Q- csi 

gpud+ = 1/2 sum tr Hdd^(-1)(x:2,y) U0(1,y)^dag  Q- csi csi^dag  Q- U0(1,z) Huu^(-1)(2,z;x)

gaud+ = - 1/2 sum tr Hdd^(-1)(x:2,y) U0(1,y)^dag  Q- csi csi^dag  Q- U0(1,z) Huu^(-1)(2,z;x) g0 

*/

  int i,j,ix0,ix1,ix2,ix3;
  double ga[GLB_T], gp[GLB_T], g1=0, temp;
  spinor_field* prop;	
  spinor_field *source;
  suNf_spinor *sptr[2];
  suNf_spinor *stmp;
  stmp=malloc((4*NF+1)*sizeof(suNf_spinor));
  suNf_spinor* sbord[2];
  sbord[0]=malloc(8*NF*sizeof(suNf_spinor));
  sbord[1]=sbord[0]+4*NF;
  suNf *uptr;
  prop=alloc_spinor_field_f(4*NF,&glattice);
  source=alloc_spinor_field_f(4*NF,&glattice);    

  for(ix0=0;ix0<GLB_T;ix0++)
    {
      ga[ix0]=0;
      gp[ix0]=0;
    }

/*
[csi_(a0,alpha0)]_{a,alpha} = \delta_{a,a0} \delta_{alpha,alpha0) 
source_s(t,x) = Q- U_0(1,x)^\dag  \csi_s \delta_{t,2}
Q- = (1 - i gamma_0 gamma_5)/2 

     | Id    -i Id |
Q- = |             | * 1/2
     | i Id     Id |
   

*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    if(COORD[0]==0){
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s%4].c[s/4].re=1.;
      stmp[0].c[s%4].c[s/4].im=0.;

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
  prop[s] = Hdd^{-1} source[s]

  */
  for(int s=0; s<4*NF; s++){
    spinor_field_zero_f(&prop[s]);
    _update_par.SF_sign=-_update_par.SF_sign;
    SF_quark_propagator(&source[s], mass, &prop[s], acc);
    _update_par.SF_sign=-_update_par.SF_sign;
  }
  
/*

gpud+ = 1/2 sum tr Hdd^(-1)(x:2,y) U0(1,y)^dag  Q- csi csi^dag  Q- U0(1,z) Huu^(-1)(2,z;x)

gaud+ = - 1/2 sum tr  g0  Hdd^(-1)(x:2,y) U0(1,y)^dag  Q- csi csi^dag  Q- U0(1,z) Huu^(-1)(2,z;x)


 */
  for(int s=0;s<4*NF;s++){
    for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      i=ipt(ix0,ix1,ix2,ix3);
      sptr[0] = _FIELD_AT(&prop[s],i);
      _spinor_prod_re_f(temp,*sptr[0],*sptr[0]);
      /*gpud*/
      gp[(COORD[0]*T+ix0-1+GLB_T)%GLB_T] += .5*(temp)/GLB_VOL3;
 
      /*gaud*/
      /*-gamma_0*/
      stmp[0].c[0]=sptr[0]->c[2];
      stmp[0].c[1]=sptr[0]->c[3];
      stmp[0].c[2]=sptr[0]->c[0];
      stmp[0].c[3]=sptr[0]->c[1];

      _spinor_prod_re_f(temp,*sptr[0],stmp[0]);
      ga[(COORD[0]*T+ix0-1+GLB_T)%GLB_T] += .5*temp/GLB_VOL3;


    }
  }
  
  global_sum((double*)ga,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_ud_+[%d] = %.10e\n",ix0,ga[ix0]);	
  lprintf("PC_twisted_AC",10,"\n");
  global_sum((double*)gp,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_ud_+[%d] = %.10e\n",ix0,gp[ix0]);	

  for(ix0=1;ix0<GLB_T-2;ix0++)
    lprintf("PC_twisted_AC",10," mpcac [%d] = %.10e\n",ix0,(ga[ix0+1] - ga[ix0-1])/(4*gp[ix0]));
  
  
  
   lprintf("PC_twisted_AC",10,"\n");


   /*
     source_s(t,x) = g5 Q- U_0(1,x)^\dag  \csi_s \delta_{t,2}
     Q- = (1 - i gamma_0 gamma_5)/2 
     Q+ = (1+ i gamma_0 gamma_5)/2 
     
     | Id    +i Id |
Q+ = |             | * 1/2
     | - i Id   Id |

    g1ud+ = 1/2 sum tr U0(1,z) Ddd^(-1)(2,z;T-2,z') U0(T-2,z') g5 Q+ U0(T-2,y')^dag Duu^(-1)(T-2,y',2,y) U0(1,y)^dag g5 Q+
           = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') g5 U0(T-2,z') g5 Q+ U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) g5 U0(1,y)^dag g5 Q+
           = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') U0(T-2,z') Q+  U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) U0(1,y)^dag Q+    
           = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') U0(T-2,z') Q+  U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) U0(1,y)^dag Q+
           = 1/2 sum tr U0(1,z) Hdd^(-1)(2,z;T-2,z') U0(T-2,z') Q+ U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) U0(1,y)^dag g5 Q- g5
           = 1/2 sum csi^+ Q- g5 U0(1,z) Hdd^(-1)(2,z;T-2,z') U0(T-2,z') Q+ Q+  U0(T-2,y')^dag Huu^(-1)(T-2,y',2,y) g5 U0(1,y)^dag Q- csi

*/


  
  /*
  prop[s] = Huu^{-1} g5 source[s]

  */
  for(int s=0; s<4*NF; s++){
    spinor_field_zero_f(&prop[s]);
    spinor_field_g5_assign_f(&source[s]);
    SF_quark_propagator(&source[s], mass, &prop[s], acc);
    spinor_field_g5_assign_f(&source[s]);
  }
  

  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(stmp[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
        i=ipt((T-2),ix1,ix2,ix3);
        sptr[0] = _FIELD_AT(&prop[s],i);
        uptr=pu_gauge_f(i,0);
        for(j=0;j<4;j++) {
          _suNf_inverse_multiply(stmp[4*NF].c[j],*uptr,sptr[0]->c[j]);
        }
        
        _spinor_add_assign_f(stmp[s],stmp[4*NF]);
        
      }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(stmp[s]);}
  }
  global_sum((double*)stmp,4*NF*4*NF*2);    
  g1=0.;
  for(int s=0;s<4*NF;s++){
    stmp[4*NF]=stmp[s];
    _vector_i_add_assign_f(stmp[4*NF].c[0],stmp[s].c[2]);
    _vector_i_add_assign_f(stmp[4*NF].c[1],stmp[s].c[3]);
    _vector_i_sub_assign_f(stmp[4*NF].c[2],stmp[s].c[0]);
    _vector_i_sub_assign_f(stmp[4*NF].c[3],stmp[s].c[1]);
    
    
    _spinor_prod_re_f(temp,stmp[s],stmp[4*NF]);
    
    
    g1+=temp;
  }
  g1*=(.5/GLB_VOL3)/GLB_VOL3;
  
  global_sum((double*)(&g1),1);
  lprintf("PC_twisted_AC",10," g1_ud_+ = %.10e\n",g1);	
  lprintf("PC_twisted_AC",10,"\n");

  free_spinor_field(prop);
  free_spinor_field(source);
  free(stmp);
  free(sbord[0]);
  
#endif
  
}
