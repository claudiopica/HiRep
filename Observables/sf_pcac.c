#include "global.h"
#include "communications.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "geometry.h"
#include "spinor_field.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"

static double hmass;
static void H(spinor_field *out, spinor_field *in){
  g5Dphi(hmass,out,in);
}

void sf_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc) {
  static MINRES_par MINRESpar;
  int cgiter;
  hmass = mass;

  MINRESpar.err2 = acc;
  MINRESpar.max_iter = 0;
  cgiter=0;
  cgiter+=MINRES(&MINRESpar, &H, in, out, 0);
	  lprintf("PROPAGATOR",10,"MINRES MVM = %d",cgiter);
}

double sf_PCAC_wall_mass(double mass)
{
	int j,ix0,ix1,ix2,ix3,source;
	double f_P[GLB_T], f_A[GLB_T], f_Pt[GLB_T], f_At[GLB_T], f_1=0, f_1t=0, temp;
	for(ix0=0;ix0<GLB_T;ix0++)
	{
		f_A[ix0]=0;
		f_At[ix0]=0;
		f_P[ix0]=0;
		f_Pt[ix0]=0;
	}
	double acc = 1.e-10;
	spinor_field *prop;	
	spinor_field *prop_source;
	suNf_spinor *S_y;
	suNf_spinor *tmp=0;
	suNf_spinor S1, S2, Stmp1P, **S_top, **S_bottom;
	S_top=(suNf_spinor **)malloc(sizeof(suNf_spinor*)*(4*NF));
	S_bottom=(suNf_spinor **)malloc(sizeof(suNf_spinor*)*(4*NF));
	for(j=0;j<4*NF;j++)
	{
	S_top[j]=malloc(sizeof(suNf_spinor));
	_spinor_zero_f(*S_top[j]);
	S_bottom[j]=malloc(sizeof(suNf_spinor));
	_spinor_zero_f(*S_bottom[j]);
	}
	suNf *U1;
	prop=alloc_spinor_field_f(4*NF,&glattice);
	prop_source=alloc_spinor_field_f(4*NF,&glattice);

   /*Create wall source with g5 factor at t=1*/
   _DECLARE_INT_ITERATOR(i);
   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop_source[source]);
	spinor_field_zero_f(&prop_source[source]);

	if(COORD[0]==0)
	{
	_PIECE_FOR((&prop_source[source])->type,i)
	{
		_SITE_FOR((&prop_source[source])->type,i)
		{

	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		if (ipt(2,ix1,ix2,ix3)==i)
		{
		  tmp = _FIELD_AT(&prop_source[source],i);
		  if (source<2*NF)
		  {
		  	*(((double *) tmp+2*source))=1.;
		  }
		  else
		  {
		  	*(((double *) tmp+2*source))=-1.;
		  }
		}
	      }
	    }
	  }
		}
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop_source[source]);
		}
	}
	}
   }
   
   /*U' and P+ on source (actually P- since there is a g5 that needs to be commuted through)*/
   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop_source[source]);
	if(COORD[0]==0)
	{

	_PIECE_FOR((&prop_source[source])->type,i)
	{
		_SITE_FOR((&prop_source[source])->type,i)
		{
	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		if (ipt(2,ix1,ix2,ix3)==i)
		{
		      U1 = pu_gauge_f(idn(i,0),0);
		      S_y = _FIELD_AT(&prop_source[source],i);
		      S1 = (*S_y);
		      for(j=0;j<4;j++)
		      {
        		 _suNf_inverse_multiply(Stmp1P.c[j],(*U1),(*S_y).c[j]);
		      }
		      _vector_lc_f((*S_y).c[0],0.5,Stmp1P.c[0],+0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[1],0.5,Stmp1P.c[1],+0.5,Stmp1P.c[3]);
		      _vector_lc_f((*S_y).c[2],+0.5,Stmp1P.c[0],0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[3],+0.5,Stmp1P.c[1],0.5,Stmp1P.c[3]);

		 }
	       }
	     }
	   }
	      }
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop_source[source]);
		}
	}
	}
    }

	/*get propagator to all points*/
	for(source=0; source<4*NF; source++)
	{
		spinor_field_zero_f(&prop[source]);
		sf_quark_propagator(&prop_source[source], mass, &prop[source], acc); 
	}

	/*get time averaged correlators for each timeslice*/
	/*f_P*/
   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop[source]);

	_PIECE_FOR((&prop[source])->type,i)
	{
		_SITE_FOR((&prop[source])->type,i)
		{
	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		for(ix0=0;ix0<T;ix0++)
		{
		  if (ipt(ix0,ix1,ix2,ix3)==i)
		  {
		      S_y = _FIELD_AT(&prop[source],i);
		      S1 = (*S_y);
		      /*f_P*/
		      _spinor_prod_re_f(temp,S1,S1);
		      f_P[(COORD[0]*T+ix0-1)%GLB_T]+=temp;
		      /*f_A*/
			/*gamma_0*/
			_vector_mul_f((S2).c[0],1.0,S1.c[2]);
			_vector_mul_f((S2).c[1],1.0,S1.c[3]);
			_vector_mul_f((S2).c[2],1.0,S1.c[0]);
			_vector_mul_f((S2).c[3],1.0,S1.c[1]);
			_spinor_prod_re_f(temp,S1,S2);
		      f_A[(COORD[0]*T+ix0-1)%GLB_T]+=temp;
		   }
		 }
	       }
	     }
	   }
	      }
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop[source]);
		}
	}

    }
   global_sum((double*)f_P,GLB_T);
   for(ix0=0;ix0<GLB_T-1;ix0++)
   {
	lprintf("PC_wall_AC",10,"f_Ppost%d = %.15f\n",ix0,f_P[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
   }
   global_sum((double*)f_A,GLB_T);
   for(ix0=0;ix0<GLB_T-1;ix0++)
   {
	lprintf("PC_wall_AC",10,"f_Apost%d = %.15f\n",ix0,f_A[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
   }

/*f_1 - NEED TO DO EACH color/dirac component separately, then combine at the end*/
   /*U' and P+ on prop at T-2 (actually P- since there is a g5 that needs to be commuted through)*/

   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop[source]);
	if(COORD[0]==NP_T-1)
	{

	_PIECE_FOR((&prop[source])->type,i)
	{
		_SITE_FOR((&prop[source])->type,i)
		{
	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		if (ipt(T-2,ix1,ix2,ix3)==i)
		{
		      U1 = pu_gauge_f(i,0);
		      S_y = _FIELD_AT(&prop[source],i);
		      S1 = (*S_y);
		      for(j=0;j<4;j++)
		      {
        		 _suNf_inverse_multiply(Stmp1P.c[j],(*U1),(*S_y).c[j]);
		      }
		      _vector_lc_f((*S_y).c[0],0.5,Stmp1P.c[0],-0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[1],0.5,Stmp1P.c[1],-0.5,Stmp1P.c[3]);
		      _vector_lc_f((*S_y).c[2],-0.5,Stmp1P.c[0],0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[3],-0.5,Stmp1P.c[1],0.5,Stmp1P.c[3]);
		      _spinor_mul_add_assign_f((*S_top[source]),1.0,(*S_y));
		 }
	       }
	     }
	   }
	      }
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop[source]);
		}
	}
	}
    }

   for(source=0;source<4*NF;source++)
   {
	global_sum((double*)S_top[source],(int)(sizeof(suNf_spinor)/sizeof(double)));
   }
	if(PID==0)
	{
		f_1=0;
		for(source=0;source<4*NF;source++)
		{
			_spinor_prod_re_f(temp,(*S_top[source]),(*S_top[source]));
			f_1+=temp;
		}
	}
	lprintf("PC_wall_AC",0,"f1_pos = %.15f\n",f_1/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z)));
	lprintf("PC_wall_AC",0,"ZP_pos = %.15f\n",(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));

	for (ix0=2;ix0<GLB_T-3;ix0++)
	{
	lprintf("PC_wall_AC",0,"PCACpost%d = %f\n",ix0,(double)(f_A[(int)(ix0)+1] - f_A[(int)(ix0)-1])/(4*f_P[(int)(ix0)]));
	}

   /*Create wall source with g5 factor at t=T-2*/
   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop_source[source]);
	spinor_field_zero_f(&prop_source[source]);

	if(COORD[0]==NP_T-1)
	{
	_PIECE_FOR((&prop_source[source])->type,i)
	{
		_SITE_FOR((&prop_source[source])->type,i)
		{

	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		if (ipt(T-2,ix1,ix2,ix3)==i)
		{
		  tmp = _FIELD_AT(&prop_source[source],i);
		  if (source<2*NF)
		  {
		  	*(((double *) tmp+2*source))=1.;
		  }
		  else
		  {
		  	*(((double *) tmp+2*source))=-1.;
		  }
		}
	      }
	    }
	  }
		}
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop_source[source]);
		}
	}
	}
   }

   /*U and P- on source (again actually use P+ to account for commuting with g5 in source)*/
   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop_source[source]);

	if(COORD[0]==NP_T-1)
	{
	_PIECE_FOR((&prop_source[source])->type,i)
	{
		_SITE_FOR((&prop_source[source])->type,i)
		{
	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		if (ipt(T-2,ix1,ix2,ix3)==i)
		{
		      U1 = pu_gauge_f(i,0);
		      S_y = _FIELD_AT(&prop_source[source],i);
		      S1 = (*S_y);
		      for(j=0;j<4;j++)
		      {
        		 _suNf_multiply(Stmp1P.c[j],(*U1),(*S_y).c[j]);
		      }
		      _vector_lc_f((*S_y).c[0],0.5,Stmp1P.c[0],-0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[1],0.5,Stmp1P.c[1],-0.5,Stmp1P.c[3]);
		      _vector_lc_f((*S_y).c[2],-0.5,Stmp1P.c[0],0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[3],-0.5,Stmp1P.c[1],0.5,Stmp1P.c[3]);

		 }
	       }
	     }
	   }
	      }
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop_source[source]);
		}
	}
	}
    }

	/*get propagator to all points*/
	for(i=0; i<4*NF; i++)
	{
		spinor_field_zero_f(&prop[i]);
		sf_quark_propagator(&prop_source[i], mass, &prop[i], acc); 
	}
   /*get time averaged correlators for each timeslice (going back from T in time)*/
   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop[source]);

	_PIECE_FOR((&prop[source])->type,i)
	{
		_SITE_FOR((&prop[source])->type,i)
		{
	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		for(ix0=0;ix0<T;ix0++)
		{
		  if (ipt(ix0,ix1,ix2,ix3)==i)
		  {
		      S_y = _FIELD_AT(&prop[source],i);
		      S1 = (*S_y);
		      /*f_P*/
		      _spinor_prod_re_f(temp,S1,S1);
		      f_Pt[((GLB_T-1)-(COORD[0]*T+ix0))%GLB_T]+=temp;
		      /*f_A*/
			/*gamma_0*/
			_vector_mul_f((S2).c[0],-1.0,S1.c[2]);
			_vector_mul_f((S2).c[1],-1.0,S1.c[3]);
			_vector_mul_f((S2).c[2],-1.0,S1.c[0]);
			_vector_mul_f((S2).c[3],-1.0,S1.c[1]);
			_spinor_prod_re_f(temp,S1,S2);
		      f_At[((GLB_T-1)-(COORD[0]*T+ix0))%GLB_T]+=temp;
		   }
		 }
	       }
	     }
	   }
	      }
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop[source]);
		}
	}
    }
   global_sum((double*)f_Pt,GLB_T);
   for(ix0=0;ix0<GLB_T-1;ix0++)
   {
	lprintf("PC_wall_AC",10,"f_Pnegt%d = %.15f\n",ix0,f_Pt[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
   }
   global_sum((double*)f_At,GLB_T);
   for(ix0=0;ix0<GLB_T-1;ix0++)
   {
	lprintf("PC_wall_AC",10,"f_Anegt%d = %.15f\n",ix0,f_At[ix0]/(double)(GLB_X*GLB_Y*GLB_Z));	
   }

/*f_1 - NEED TO DO EACH color/dirac component separately, then combine at the end*/
   /*U' and P+ on prop at T-2 (actually P- since there is a g5 that needs to be commuted through)*/

   for(source=0;source<4*NF;source++)
   {
	start_sf_sendrecv(&prop[source]);
	if(COORD[0]==0)
	{

	_PIECE_FOR((&prop[source])->type,i)
	{
		_SITE_FOR((&prop[source])->type,i)
		{
	for(ix1=0;ix1<X;ix1++)
	  {
	    for(ix2=0;ix2<Y;ix2++)
	    {
	      for(ix3=0;ix3<Z;ix3++)
	      {
		if (ipt(2,ix1,ix2,ix3)==i)
		{
		      U1 = pu_gauge_f(idn(i,0),0);
		      S_y = _FIELD_AT(&prop[source],i);
		      S1 = (*S_y);
		      for(j=0;j<4;j++)
		      {
        		 _suNf_multiply(Stmp1P.c[j],(*U1),(*S_y).c[j]);
		      }
		      _vector_lc_f((*S_y).c[0],0.5,Stmp1P.c[0],+0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[1],0.5,Stmp1P.c[1],+0.5,Stmp1P.c[3]);
		      _vector_lc_f((*S_y).c[2],+0.5,Stmp1P.c[0],0.5,Stmp1P.c[2]);
		      _vector_lc_f((*S_y).c[3],+0.5,Stmp1P.c[1],0.5,Stmp1P.c[3]);
		      _spinor_mul_add_assign_f((*S_bottom[source]),1.0,(*S_y));
		 }
	       }
	     }
	   }
	      }
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(&prop[source]);
		}
	}
	}
    }

   for(source=0;source<4*NF;source++)
   {
	global_sum((double*)S_bottom[source],(int)(sizeof(suNf_spinor)/sizeof(double)));
   }
	if(PID==0)
	{
		f_1t=0;
		for(source=0;source<4*NF;source++)
		{
			_spinor_prod_re_f(temp,(*S_bottom[source]),(*S_bottom[source]));
			f_1t+=temp;
		}
	}
	lprintf("PC_wall_AC",0,"f1_neg = %.15f\n",f_1t/((double)(GLB_X*GLB_Y*GLB_Z*GLB_X*GLB_Y*GLB_Z)));
	lprintf("PC_wall_AC",0,"ZP_neg = %.15f\n",(sqrt(f_1t)/(f_Pt[(int)(GLB_X/2)])));

	lprintf("PC_wall_AC",0,"Z_P = %.15f\n",0.5*(sqrt(f_1t)/(f_Pt[(int)(GLB_X/2)]))+0.5*(sqrt(f_1)/(f_P[(int)(GLB_X/2)])));

	for (ix0=2;ix0<GLB_T-3;ix0++)
	{
	lprintf("PC_wall_AC",0,"PCACnegt%d = %f\n",ix0,(double)(f_At[(int)(ix0)+1] - f_At[(int)(ix0)-1])/(4*f_Pt[(int)(ix0)]));
	}

	free_spinor_field(prop_source);
	free_spinor_field(prop);

	return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]);
}
