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
#include "update.h"
#include "matrix.h"

extern rhmc_par _update_par; /* Update/update_rhmc.c */

static double hmass;
static void H2_sf(spinor_field *out, spinor_field *in){
  g5Dphi_sq(hmass,out,in);
}

void SF_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc) {
  int cgiter;
  static mshift_par inv_par;
  static spinor_field *chi=NULL;
  if(chi==NULL) chi=alloc_spinor_field_f(1,&glattice);
  hmass = mass;

  inv_par.n = 1;
  inv_par.shift = malloc(sizeof(double));
  inv_par.shift[0] = 0.;
  inv_par.err2= acc;
  inv_par.max_iter=0; /* no limit */
  _update_par.SF_sign = -_update_par.SF_sign;
   g5Dphi(mass, chi, in);
  _update_par.SF_sign = -_update_par.SF_sign;
  cgiter=cg_mshift(&inv_par, &H2_sf, chi, out);
	lprintf("PROPAGATOR",10,"cg_mshift MVM = %d",cgiter);
	free(inv_par.shift);
}



/*
f1 can be computed inserting the sources in the correlator in the lower or higher border. The two procedures give identical results. We keep only the computation with sources on the lower border.
*/
double SF_PCAC_wall_mass(double mass)
{
#ifdef BASIC_SF

  int i,ix0,ix1,ix2,ix3;
  double f_P[GLB_T], f_A[GLB_T], f_Pt[GLB_T], f_At[GLB_T], f_1=0, k_1=0, temp;
  complex z1,z2;
  double acc = 1.e-10;
  spinor_field *prop_x0;	
  prop_x0=alloc_spinor_field_f(4*NF,&glattice);
  spinor_field *prop_xT;	
  prop_xT=alloc_spinor_field_f(4*NF,&glattice);
  spinor_field *source;
  source=alloc_spinor_field_f(4*NF,&glattice);    
  suNf_spinor* sptr;
  suNf_spinor* stmp;
  stmp=malloc(2*sizeof(suNf_spinor));
  suNf_spinor* H_T0;
  H_T0=malloc(4*NF*sizeof(suNf_spinor));
  suNf_spinor** H_x0;
	H_x0=(suNf_spinor **)malloc(4*NF*sizeof(suNf_spinor*));
  suNf_spinor** H_xT;
	H_xT=(suNf_spinor **)malloc(4*NF*sizeof(suNf_spinor*));
  suNf *uptr;
  suNf_spinor* ptmp1;
  ptmp1=malloc(4*NF*sizeof(suNf_spinor));
  suNf_spinor* ptmp2;
  ptmp2=malloc(4*NF*sizeof(suNf_spinor));
  suNf_spinor* A;
  A=malloc(4*NF*sizeof(suNf_spinor));
  suNf_spinor* B;
  B=malloc(4*NF*sizeof(suNf_spinor));

  for(ix0=0;ix0<GLB_T;ix0++)
    { 
      f_A[ix0]=0; 
      f_At[ix0]=0;
      f_P[ix0]=0;
      f_Pt[ix0]=0;
    }

  /*Get site-to-wall propagator from each site x to all points on t=0 boundary*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    spinor_field_zero_f(&prop_x0[s]);
    if(COORD[0]==0){
      //Create g5 spinor
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s/NF].c[s%NF].re=1.;
      stmp[0].c[s/NF].c[s%NF].im=0.;
      _spinor_g5_assign_f(stmp[0]);
      //Wall source at t=0 boundary - U' and P- on g5 spinor    
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(idn(i,0),0);
      	for(int j=0;j<4;j++) {
      	  _suNf_inverse_multiply(stmp[1].c[j],*uptr,stmp[0].c[j]);
	}
      	sptr = _FIELD_AT(&source[s],i);
      	_spinor_pminus_f(*sptr,stmp[1]);
      }
    }
    //Invert source
    SF_quark_propagator(&source[s], mass, &prop_x0[s], acc); 
  }
  
  /*Get wall-to-wall boundary to boundary propagator*/
  for(int s=0;s<4*NF;s++){
    _spinor_zero_f(H_T0[s]);
    if(COORD[0]==NP_T-1){
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      	i=ipt(T-2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(i,0);
      	sptr = _FIELD_AT(&prop_x0[s],i);
      	for(int j=0;j<4;j++) {
      	  _suNf_inverse_multiply(stmp[0].c[j],*uptr,sptr->c[j]);
      	}
      	_spinor_pplus_f(stmp[1],stmp[0]);
      	_spinor_add_assign_f(H_T0[s],stmp[1]);
      }
			_spinor_g5_assign_f(H_T0[s]); //added g5 prefactor to match alpha definition
    }
  }
  global_sum((double*)H_T0,sizeof(suNf_spinor)/sizeof(double)*4*NF);

  /*Get site-to-wall propagator from each site x to all points on t=T boundary*/
  for(int s=0;s<4*NF;s++){
    spinor_field_zero_f(&source[s]);
    spinor_field_zero_f(&prop_xT[s]);
    if(COORD[0]==NP_T-1){
      //Create g5 spinor
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s/NF].c[s%NF].re=1.;
      stmp[0].c[s/NF].c[s%NF].im=0.;
      _spinor_g5_assign_f(stmp[0]);
      //Wall source at t=T boundary - U' and P+ on g5 spinor    
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++){
      	i=ipt(T-2,ix1,ix2,ix3);
      	uptr=pu_gauge_f(i,0);
      	for(int j=0;j<4;j++) {
      	  _suNf_multiply(stmp[1].c[j],*uptr,stmp[0].c[j]);
    	  }
      	sptr = _FIELD_AT(&source[s],i);
      	_spinor_pplus_f(*sptr,stmp[1]);
        }
      }
    //Invert source
    SF_quark_propagator(&source[s], mass, &prop_xT[s], acc); 
    }

  /*Get time averaged forward correlators for each timeslice*/
  for(int s=0;s<4*NF;s++){
    for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      i=ipt(ix0,ix1,ix2,ix3);
      sptr = _FIELD_AT(&prop_x0[s],i);
      //f_P = prop_x0^dag prop_x0
      _spinor_prod_re_f(temp,*sptr,*sptr);
      f_P[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]+=0.5*temp;
      //f_A = -prop_x0^dag gamma0 prop_x0
      _spinor_g0_f(stmp[0],*sptr);
      _spinor_prod_re_f(temp,*sptr,stmp[0]);
      f_A[(COORD[0]*T+ix0-1+GLB_T)%GLB_T]-=0.5*temp;
    }
  }  
  global_sum((double*)f_P,GLB_T);
  global_sum((double*)f_A,GLB_T);
  
  /*Get boundary-to-boundary correlators*/
  //f_1 = |H_T0|^2
  if(PID==0){
    f_1=0;
    for(int s=0;s<4*NF;s++){
      _spinor_prod_re_f(temp,H_T0[s],H_T0[s]);
      f_1+=temp;
    }
		f_1/=2.0;
  }
  //k_1 = |gamma_5 gamma_k H_T0|^2
  if(PID==0){
    k_1=0;
    for(int s=0;s<4*NF;s++){
	//gamma_5(\gamma_1 + \gamma_2 + \gamma_3)
	z1.re=0.;
	z1.im=-1.;
	z2.re=-1.;
	z2.im=-1.;
	_vector_clc_f(stmp[0].c[0],z1,H_T0[s].c[2],z2,H_T0[s].c[3]);
	_vector_clc_f(stmp[0].c[1],z2,H_T0[s].c[2],z1,H_T0[s].c[3]);
	_vector_clc_f(stmp[0].c[2],z1,H_T0[s].c[0],z2,H_T0[s].c[1]);
	_vector_clc_f(stmp[0].c[3],z2,H_T0[s].c[0],z1,H_T0[s].c[1]);
      _spinor_prod_re_f(temp,stmp[0],stmp[0]);
      k_1+=temp;
    }
		k_1/=6.0;
  }

  /*Get time averaged backwards correlators for each timeslice*/
  for(int s=0;s<4*NF;s++){
    for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) for(ix0=0;ix0<T;ix0++){
      i=ipt(ix0,ix1,ix2,ix3);
      sptr = _FIELD_AT(&prop_xT[s],i);
      //f_P = prop_x0^dag prop_x0
      _spinor_prod_re_f(temp,*sptr,*sptr);
      f_Pt[((GLB_T-1)-(COORD[0]*T+ix0))%GLB_T]+=0.5*temp;
      //f_A = prop_x0^dag gamma0 prop_x0
      _spinor_g0_f(stmp[0],*sptr);
      _spinor_prod_re_f(temp,*sptr,stmp[0]);
      f_At[((GLB_T-1)-(COORD[0]*T+ix0))%GLB_T]+=0.5*temp;
    }
  }
  global_sum((double*)f_Pt,GLB_T);
  global_sum((double*)f_At,GLB_T);
  
  /*Output correlators*/
  lprintf("PC_wall_AC",0,"f1 = %.10e\n",f_1/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z)));
  lprintf("PC_wall_AC",0,"k1 = %.10e\n",k_1/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z)));
  lprintf("PC_wall_AC",0,"ZP_pos = %.10e\n",(sqrt(3*f_1)/(f_P[(int)(GLB_X/2)])));
  lprintf("PC_wall_AC",0,"ZP_neg = %.10e\n",(sqrt(3.0*f_1)/(f_Pt[(int)(GLB_X/2)])));  
  lprintf("PC_wall_AC",0,"Z_P = %.10e\n",(0.5*sqrt(3.0*f_1)/(f_Pt[(int)(GLB_X/2)]))+(0.5*sqrt(3.0*f_1)/(f_P[(int)(GLB_X/2)])));
  for(ix0=0;ix0<GLB_T-1;ix0++){
    lprintf("PC_wall_AC",10,"f_Apost%d = %.10e\n",ix0,f_A[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
    lprintf("PC_wall_AC",10,"f_Ppost%d = %.10e\n",ix0,f_P[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
    lprintf("PC_wall_AC",10,"f_Anegt%d = %.10e\n",ix0,f_At[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));
    lprintf("PC_wall_AC",10,"f_Pnegt%d = %.10e\n",ix0,f_Pt[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
  }
  for (ix0=2;ix0<GLB_T-3;ix0++){
    lprintf("PC_wall_AC",0,"PCACpost%d = %10e\n",ix0,(double)(f_A[(int)(ix0)+1] - f_A[(int)(ix0)-1])/(4.0*f_P[(int)(ix0)]));
    lprintf("PC_wall_AC",0,"PCACnegt%d = %10e\n",ix0,(double)(f_At[(int)(ix0)+1] - f_At[(int)(ix0)-1])/(4.0*f_Pt[(int)(ix0)]));
  }

//Four-fermion correlators, assumes ix0, H_x0, H_xT and H_T0 exist as inputs, A,B as outputs, and uses temporary variables ptmp1, ptmp2. Assigns
//A = [H gA g5 Hdag]
//B = [H gB g5 Hcaldag gC g5 H'dag]
//then calls _TRACE
#define _FCORR_minus(_Findex,_gA,_gB,_gC) \
{ \
			_smatrix_g5ndagp(ptmp1,H_x0); \
			_smatrix_ ## _gA ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(A,H_x0,ptmp2); \
			_smatrix_g5ndagp(ptmp1,H_xT); \
			_smatrix_ ## _gC ## _f(ptmp2,ptmp1); \
			_smatrix_mdag_n(ptmp1,H_T0,ptmp2); \
			_smatrix_g5_assign_f(ptmp1); \
			_smatrix_ ## _gB ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(B,H_x0,ptmp2); \
			_TRACE_minus(A,B,F_disc[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T],F_conn[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T]); \
}((void)0)
#define _FCORR_plus(_Findex,_gA,_gB,_gC) \
{ \
			_smatrix_g5ndagp(ptmp1,H_x0); \
			_smatrix_ ## _gA ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(A,H_x0,ptmp2); \
			_smatrix_g5ndagp(ptmp1,H_xT); \
			_smatrix_ ## _gC ## _f(ptmp2,ptmp1); \
			_smatrix_mdag_n(ptmp1,H_T0,ptmp2); \
			_smatrix_g5_assign_f(ptmp1); \
			_smatrix_ ## _gB ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(B,H_x0,ptmp2); \
			_TRACE_plus(A,B,F_disc[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T],F_conn[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T]); \
}((void)0)
//Four-fermion correlators: TIME REVERSAL (using the fact that H_0T = H_T0^dag)
#define _FCORR_minus_t(_Findex,_gA,_gB,_gC) \
{ \
			_smatrix_g5ndagp(ptmp1,H_xT); \
			_smatrix_ ## _gA ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(A,H_xT,ptmp2); \
			_smatrix_g5ndagp(ptmp1,H_x0); \
			_smatrix_ ## _gC ## _f(ptmp2,ptmp1); \
			_smatrix_m_n(ptmp1,H_T0,ptmp2); \
			_smatrix_g5_assign_f(ptmp1); \
			_smatrix_ ## _gB ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(B,H_xT,ptmp2); \
			_TRACE_minus(A,B,F_disc_t[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T],F_conn_t[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T]); \
}((void)0)
#define _FCORR_plus_t(_Findex,_gA,_gB,_gC) \
{ \
			_smatrix_g5ndagp(ptmp1,H_xT); \
			_smatrix_ ## _gA ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(A,H_xT,ptmp2); \
			_smatrix_g5ndagp(ptmp1,H_x0); \
			_smatrix_ ## _gC ## _f(ptmp2,ptmp1); \
			_smatrix_m_n(ptmp1,H_T0,ptmp2); \
			_smatrix_g5_assign_f(ptmp1); \
			_smatrix_ ## _gB ## _f(ptmp2,ptmp1); \
			_smatrix_mp_n(B,H_xT,ptmp2); \
			_TRACE_plus(A,B,F_disc_t[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T],F_conn_t[_Findex-1][(COORD[0]*T+ix0-1+GLB_T)%GLB_T]); \
}((void)0)


	complex F_disc[5][GLB_T];
	complex F_conn[5][GLB_T];
	complex F_disc_t[5][GLB_T];
	complex F_conn_t[5][GLB_T];
  for(int i=0;i<5;i++)
	{
	  for(ix0=0;ix0<GLB_T;ix0++)
		{
			_complex_0(F_disc[i][ix0]);
			_complex_0(F_conn[i][ix0]);
			_complex_0(F_disc_t[i][ix0]);
			_complex_0(F_conn_t[i][ix0]);
		}
	}

	for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++)
	{
		//get H(x) and H'(x)
		i=ipt(ix0,ix1,ix2,ix3);
		for(int s=0;s<4*NF;s++)
		{
			H_x0[s] = _FIELD_AT(&prop_x0[s],i);
			H_xT[s] = _FIELD_AT(&prop_xT[s],i);
		}

		//F1: g5g5g5
		_FCORR_plus(1,g5,g5,g5);
		_FCORR_plus_t(1,g5,g5,g5);

		//F2: g1g2g3+g2g3g1+g3g1g2-g3g2g1-g2g1g3-g1g3g2
		_FCORR_plus(2,g1,g2,g3);
		_FCORR_plus(2,g2,g3,g1);
		_FCORR_plus(2,g3,g1,g2);
		_FCORR_minus(2,g1,g3,g2);
		_FCORR_minus(2,g3,g2,g1);
		_FCORR_minus(2,g2,g1,g3);
		_FCORR_plus_t(2,g1,g2,g3);
		_FCORR_plus_t(2,g2,g3,g1);
		_FCORR_plus_t(2,g3,g1,g2);
		_FCORR_minus_t(2,g1,g3,g2);
		_FCORR_minus_t(2,g3,g2,g1);
		_FCORR_minus_t(2,g2,g1,g3);

		//F3: g5g1g1+g5g2g2+g5g3g3
		_FCORR_plus(3,g5,g1,g1);
		_FCORR_plus(3,g5,g2,g2);
		_FCORR_plus(3,g5,g3,g3);
		_FCORR_plus_t(3,g5,g1,g1);
		_FCORR_plus_t(3,g5,g2,g2);
		_FCORR_plus_t(3,g5,g3,g3);

		//F4: g1g5g1+g2g5g2+g3g5g3
		_FCORR_plus(4,g1,g5,g1);
		_FCORR_plus(4,g2,g5,g2);
		_FCORR_plus(4,g3,g5,g3);
		_FCORR_plus_t(4,g1,g5,g1);
		_FCORR_plus_t(4,g2,g5,g2);
		_FCORR_plus_t(4,g3,g5,g3);

		//F5: g1g1g5+g2g2g5+g3g3g5
		_FCORR_plus(5,g1,g1,g5);
		_FCORR_plus(5,g2,g2,g5);
		_FCORR_plus(5,g3,g3,g5);
		_FCORR_plus_t(5,g1,g1,g5);
		_FCORR_plus_t(5,g2,g2,g5);
		_FCORR_plus_t(5,g3,g3,g5);

	}

  global_sum((double*)F_disc,5*GLB_T*2);
  global_sum((double*)F_conn,5*GLB_T*2);
  global_sum((double*)F_disc_t,5*GLB_T*2);
  global_sum((double*)F_conn_t,5*GLB_T*2);

	//Normalisation of F
#define _FNORM(_Findex,_ix0,_norm) \
{ \
		F_disc[_Findex-1][_ix0].re /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_disc[_Findex-1][_ix0].im /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_conn[_Findex-1][_ix0].re /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_conn[_Findex-1][_ix0].im /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_disc_t[_Findex-1][_ix0].re /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_disc_t[_Findex-1][_ix0].im /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_conn_t[_Findex-1][_ix0].re /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
		F_conn_t[_Findex-1][_ix0].im /= (double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z*_norm); \
}((void)0)

  for(ix0=0;ix0<GLB_T;ix0++)
	{
		_FNORM(1,ix0,1.0);
		_FNORM(2,ix0,6.0);
		_FNORM(3,ix0,3.0);
		_FNORM(4,ix0,3.0);
		_FNORM(5,ix0,3.0);
	}

  for(int Findex=0;Findex<5;Findex++)
	{
	  for(ix0=0;ix0<GLB_T-1;ix0++)
		{
	    lprintf("PC_wall_AC",10,"F%dplus_t%dre = %.10e\n",Findex+1,ix0,0.5*(F_disc[Findex][ix0].re+F_conn[Findex][ix0].re));	
	    lprintf("PC_wall_AC",10,"F%dplus_t%dre_t = %.10e\n",Findex+1,ix0,0.5*(F_disc_t[Findex][GLB_T-ix0-2].re+F_conn_t[Findex][GLB_T-ix0-2].re));	
	    lprintf("PC_wall_AC",10,"F%dplus_t%dim = %.10e\n",Findex+1,ix0,0.5*(F_disc[Findex][ix0].im+F_conn[Findex][ix0].im));	
	    lprintf("PC_wall_AC",10,"F%dplus_t%dim_t = %.10e\n",Findex+1,ix0,0.5*(F_disc_t[Findex][GLB_T-ix0-2].im+F_conn_t[Findex][GLB_T-ix0-2].im));	
	    lprintf("PC_wall_AC",10,"F%dminus_t%dre = %.10e\n",Findex+1,ix0,0.5*(F_disc[Findex][ix0].re-F_conn[Findex][ix0].re));	
	    lprintf("PC_wall_AC",10,"F%dminus_t%dre_t = %.10e\n",Findex+1,ix0,0.5*(F_disc_t[Findex][GLB_T-ix0-2].re-F_conn_t[Findex][GLB_T-ix0-2].re));	
	    lprintf("PC_wall_AC",10,"F%dminus_t%dim = %.10e\n",Findex+1,ix0,0.5*(F_disc[Findex][ix0].im-F_conn[Findex][ix0].im));
	    lprintf("PC_wall_AC",10,"F%dminus_t%dim_t = %.10e\n",Findex+1,ix0,0.5*(F_disc_t[Findex][GLB_T-ix0-2].im-F_conn_t[Findex][GLB_T-ix0-2].im));
		}
	}

  free_spinor_field(source);
  free_spinor_field(prop_x0);
  free_spinor_field(prop_xT);

  free(stmp);
  free(ptmp1);
  free(ptmp2);
  free(H_T0);
  free(H_x0);
  free(H_xT);
  free(A);
  free(B);

  return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]);
  
  
#elif defined(ROTATED_SF)

  int i,j,ix0,ix1,ix2,ix3;
  double f_P[GLB_T], f_A[GLB_T], f_1=0, temp;
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
      f_P[ix0]=0;
    }

/*
[eta_(a0,alpha0)]_{a,alpha} = \delta_{a,a0} \delta_{alpha,alpha0) 
source_s(t,x) = Q- U_0(1,y)^\dag gamma_5 \eta_s \delta_{t,1}
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
  f_P(x) = -1/2 \sum_s prop[0][s]^dag(x) gamma_5 prop[1][s](x)
  f_A(x) = -1/2 \sum_s prop[0][s]^dag(x) gamma_0 prop[1][s](x)
  */
  for(int s=0;s<4*NF;s++){
    for(ix0=0;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
      i=ipt(ix0,ix1,ix2,ix3);
      sptr[0] = _FIELD_AT(&prop[0][s],i);
      sptr[1] = _FIELD_AT(&prop[1][s],i);
      /*f_P*/
      _spinor_g5_prod_re_f(temp,*sptr[0],*sptr[1]);
      f_P[(COORD[0]*T+ix0-1+GLB_T)%GLB_T] += -.5*temp;
      /*f_A*/
      /*gamma_0*/
      stmp[0].c[0]=sptr[1]->c[2];
      stmp[0].c[1]=sptr[1]->c[3];
      stmp[0].c[2]=sptr[1]->c[0];
      stmp[0].c[3]=sptr[1]->c[1];
      _spinor_prod_re_f(temp,*sptr[0],stmp[0]);
      f_A[(COORD[0]*T+ix0-1+GLB_T)%GLB_T] += .5*temp;
    }
  }
  
  global_sum((double*)f_P,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Ppost%d = %.10e\n",ix0,f_P[ix0]);	

  global_sum((double*)f_A,GLB_T);
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_wall_AC",10,"f_Apost%d = %.10e\n",ix0,f_A[ix0]);	
  

  /*
  sbord[0][s] = \sum_v gamma_5 Q- U_0(T-2,v)^dag prop[0][s](T-2,v)
  sbord[1][s] = \sum_v U_0(T-2,v)^dag prop[1][s](T-2,v)
  f1 = 1/(2*VOL^2) \sum_s sbord[0][s]^dag sbord[1][s]
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
      _spinor_prod_re_f(temp,sbord[0][s],sbord[1][s]);
      f_1+=.5*temp/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z));
    }
  }
  lprintf("PC_wall_AC",0,"f1 = %.10e\n",f_1);
  lprintf("PC_wall_AC",0,"ZP_pos = %.10e\n",(sqrt(f_1)/(f_P[(int)(GLB_T/2-1)])));
  
  for (ix0=2;ix0<GLB_T-3;ix0++)
    lprintf("PC_wall_AC",0,"PCACpost%d = %f\n",ix0,(double)(f_A[(int)(ix0)+1] - f_A[(int)(ix0)-1])/(4*f_P[(int)(ix0)]));
  


/*As question of principle it  could be added also the calculation of the correlator with the upper border, these would leade to a double measure of the mpcac, in this moment we are not implemting it, but in case of need should be easy to use the different definition of the upper border to add it*/  
 
  return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]);
  
  free(stmp);
  free(sbord[0]);
  free_spinor_field(prop[0]);
  free_spinor_field(prop[1]);
  free_spinor_field(source );
 
#else

  return 0.;

#endif
  
}
