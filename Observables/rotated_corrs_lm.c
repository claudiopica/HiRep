#ifdef ROTATED_SF
#include "global.h"
#include "communications.h"
#include "observables.h"
#include "logger.h"
void rotated_lXuum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *lauu_m=corr_mem->g1, *lvuu_m=corr_mem->g2, *ltuu_m=corr_mem->g3, *lttuu_m=corr_mem->g4, **M=corr_mem->M;

  complex *lauu_1m=corr_mem->l11, *lvuu_1m=corr_mem->l21, *ltuu_1m=corr_mem->l31, *lttuu_1m=corr_mem->l41;
  complex *lauu_2m=corr_mem->l12, *lvuu_2m=corr_mem->l22, *ltuu_2m=corr_mem->l32, *lttuu_2m=corr_mem->l42;
  complex *lauu_3m=corr_mem->l13, *lvuu_3m=corr_mem->l23, *ltuu_3m=corr_mem->l33, *lttuu_3m=corr_mem->l43;

  complex ***la1_ij=corr_mem->l11_ij, ***la2_ij=corr_mem->l12_ij, ***la3_ij=corr_mem->l13_ij;
  complex ***lv1_ij=corr_mem->l21_ij, ***lv2_ij=corr_mem->l22_ij, ***lv3_ij=corr_mem->l23_ij;
  complex ***lt1_ij=corr_mem->l31_ij, ***lt2_ij=corr_mem->l32_ij, ***lt3_ij=corr_mem->l33_ij;
  complex ***ltt1_ij=corr_mem->g1_ij,  ***ltt2_ij=corr_mem->g2_ij,  ***ltt3_ij=corr_mem->g3_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;

/*********************   lY correlations  

 Note: Since g1 and g3 appear allways twice in the correlation functions, we can simply ignore the imaginary
      components ans multiply the correlation by an overall (-1).                                             

***************************/


  /**************************************************
    set:   lauu_m, lpuu_m, lsuu_m, lvuu_m

lauu- = - 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   gk     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q+ csi  

    lauu_1- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g1     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q+ csi
    lauu_2- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g2     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q+ csi
    lauu_3- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g3     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q+ csi

lvuu- = 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5gk     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q+ csi  

    lvuu_1- = 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5g1     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q+ csi
    lvuu_2- = 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5g2     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q+ csi
    lvuu_3- = 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5g3     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q+ csi

ltuu- = - i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig0k     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q+ csi  

    ltuu_1- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig01     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q+ csi
    ltuu_2- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig02     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q+ csi
    ltuu_3- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig03     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q+ csi

lttuu- = - i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)    sig0k     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q+ csi  

    lttuu_1- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   sig01     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q+ csi
    lttuu_2- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   sig02     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q+ csi
    lttuu_3- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   sig03     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q+ csi

															(sig0k=ig0gk)
  **************************************************/

  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(lauu_m[ix0]);
      _complex_0(lauu_1m[ix0]);
      _complex_0(lauu_2m[ix0]);
      _complex_0(lauu_3m[ix0]);

      _complex_0(lvuu_m[ix0]);
      _complex_0(lvuu_1m[ix0]);
      _complex_0(lvuu_2m[ix0]);
      _complex_0(lvuu_3m[ix0]);

      _complex_0(ltuu_m[ix0]);
      _complex_0(ltuu_1m[ix0]);
      _complex_0(ltuu_2m[ix0]);
      _complex_0(ltuu_3m[ix0]);

      _complex_0(lttuu_m[ix0]);
      _complex_0(lttuu_1m[ix0]);
      _complex_0(lttuu_2m[ix0]);
      _complex_0(lttuu_3m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(la1_ij[ix0][s1][s2]);
	  _complex_0(la2_ij[ix0][s1][s2]);
	  _complex_0(la3_ij[ix0][s1][s2]);

	  _complex_0(lv1_ij[ix0][s1][s2]);
	  _complex_0(lv2_ij[ix0][s1][s2]);
	  _complex_0(lv3_ij[ix0][s1][s2]);

	  _complex_0(lt1_ij[ix0][s1][s2]);
	  _complex_0(lt2_ij[ix0][s1][s2]);
	  _complex_0(lt3_ij[ix0][s1][s2]);

	  _complex_0(ltt1_ij[ix0][s1][s2]);
	  _complex_0(ltt2_ij[ix0][s1][s2]);
	  _complex_0(ltt3_ij[ix0][s1][s2]);
	}
    }


/********            lY_1           *******************/
  /*  Gamma matrix structure:    \csi^dag_j  g5 g1 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g1/i Q+*/
        /*Q+*/
        stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g1/i*/
        stmp1[1].c[0]=stmp1[0].c[3];
        stmp1[1].c[1]=stmp1[0].c[2];
        stmp1[1].c[2]=stmp1[0].c[1];
        stmp1[1].c[3]=stmp1[0].c[0];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[1],-1.0,stmp1[1].c[1]);
        _vector_mul_f(stmp1[1].c[2],-1.0,stmp1[1].c[2]);
        _vector_mul_f(stmp1[1].c[3],-1.0,stmp1[1].c[3]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l1_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	      /*lauu_1m*/
	      /*gamma_1/i*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 


	      /*lvuu_1m*/
	      /*g5 gamma_1/i*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 

	      /*ltuu_1m*/
	      /*g5 sig01*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 

	      /*lttuu_1m*/
	      /*sig01*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     
	    }
        // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
	_complex_mulr(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }

  /***** Here we do the contractions of lY1_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lauu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvuu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltuu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttuu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}


/********            lY_2           *******************/
  /*  Gamma matrix structure:    \csi^dag_j  g5 g2 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g2 Q+*/
        /*Q+*/
        stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g2*/
        stmp1[1].c[0]=stmp1[0].c[3];
        stmp1[1].c[1]=stmp1[0].c[2];
        stmp1[1].c[2]=stmp1[0].c[1];
        stmp1[1].c[3]=stmp1[0].c[0];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[2],-1.0,stmp1[1].c[2]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l2_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	      /*lauu_2m*/
	      /*gamma_2*/

	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

	      /*lvuu_2m*/
	      /*g5 gamma_2*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 

	      /*ltuu_2m*/
	      /*g5 sig02*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 

	      /*lttuu_2m*/
	      /*sig02*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     


	    }
        // Normalization with -1.0/(6.0*GLB_VOL3)
	_complex_mulr(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }

  /***** Here we do the contractions of lY2_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lauu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvuu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltuu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttuu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}




/********            lY_3           *******************/
  /*  Gamma matrix structure:    \csi^dag_j  g5 g3 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g3/i Q+*/
        /*Q+*/
        stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g3/i*/

        stmp1[1].c[0]=stmp1[0].c[2];
        stmp1[1].c[1]=stmp1[0].c[3];
        stmp1[1].c[2]=stmp1[0].c[0];
        stmp1[1].c[3]=stmp1[0].c[1];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[2],-1.0,stmp1[1].c[2]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l3_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){

	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	      /*lauu_3m*/
	      /*gamma_3/i*/

	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     



	      /*lvuu_3m*/
	      /*g5 gamma_3/i*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 

	      /*ltuu_3m*/
	      /*g5 sig03*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp); 

	      /*lttuu_3m*/
	      /*sig03*/
  	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     



	    }
        // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
	_complex_mulr(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }

  /***** Here we do the contractions of lY3_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lauu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvuu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltuu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttuu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}


/*******   Sum over k ************/
   for(ix0=0;ix0<T;ix0++){
	_complex_add_assign(lauu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lauu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lauu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lauu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lauu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lauu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);

	_complex_add_assign(lvuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvuu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lvuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvuu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lvuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvuu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);

	_complex_add_assign(ltuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltuu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ltuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltuu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ltuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltuu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);

	_complex_add_assign(lttuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttuu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lttuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttuu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lttuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttuu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)lauu_m,2*GLB_T);
  global_sum((double*)lvuu_m,2*GLB_T);
  global_sum((double*)ltuu_m,2*GLB_T);
  global_sum((double*)lttuu_m,2*GLB_T);

  /************************************ END of set lY_uu_-*********************/




  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," la_uu_-[%d] = %.10e,%.10e\n",ix0,lauu_m[ix0].re,lauu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lv_uu_-[%d] = %.10e,%.10e\n",ix0,lvuu_m[ix0].re,lvuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lt_uu_-[%d] = %.10e,%.10e\n",ix0,ltuu_m[ix0].re,ltuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ltt_uu_-[%d] = %.10e,%.10e\n",ix0,lttuu_m[ix0].re,lttuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}





 


void rotated_lXddm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *ladd_m=corr_mem->g1, *lvdd_m=corr_mem->g2, *ltdd_m=corr_mem->g3, *lttdd_m=corr_mem->g4, **M=corr_mem->M;

  complex *ladd_1m=corr_mem->l11, *lvdd_1m=corr_mem->l21, *ltdd_1m=corr_mem->l31, *lttdd_1m=corr_mem->l41;
  complex *ladd_2m=corr_mem->l12, *lvdd_2m=corr_mem->l22, *ltdd_2m=corr_mem->l32, *lttdd_2m=corr_mem->l42;
  complex *ladd_3m=corr_mem->l13, *lvdd_3m=corr_mem->l23, *ltdd_3m=corr_mem->l33, *lttdd_3m=corr_mem->l43;

  complex ***la1_ij=corr_mem->l11_ij, ***la2_ij=corr_mem->l12_ij, ***la3_ij=corr_mem->l13_ij;
  complex ***lv1_ij=corr_mem->l21_ij, ***lv2_ij=corr_mem->l22_ij, ***lv3_ij=corr_mem->l23_ij;
  complex ***lt1_ij=corr_mem->l31_ij, ***lt2_ij=corr_mem->l32_ij, ***lt3_ij=corr_mem->l33_ij;
  complex ***ltt1_ij=corr_mem->g1_ij,  ***ltt2_ij=corr_mem->g2_ij,  ***ltt3_ij=corr_mem->g3_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;



   /**************************************************

      set:   ladd_m, lpdd_m, lsdd_m, lvd_m


  ladd- = - 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   gk     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q- csi

      ladd_1- = - 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g1     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q- csi
      ladd_2- = - 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g2     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q- csi
      ladd_3- = - 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g3     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q- csi



  lvdd- = 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5gk     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q- csi

      lvdd_1- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5g1     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q- csi
      lvdd_2- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5g2     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q- csi
      lvdd_3- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5g3     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q- csi

  ltdd- = i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig0k     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q- csi

      ltdd_1- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig01     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q- csi
      ltdd_2- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig02     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q- csi
      ltdd_3- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig03     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q- csi

  lttdd- = i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)    sig0k     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 gk Q- csi

      lttdd_1- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   sig01     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g1 Q- csi
      lttdd_2- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   sig02     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g2 Q- csi
      lttdd_3- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   sig03     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g5 g3 Q- csi

  															(sig0k=ig0gk)
    **************************************************/

    /****Initialization*****/
    for(ix0=0;ix0<GLB_T;ix0++)
      {
        _complex_0(ladd_m[ix0]);
        _complex_0(ladd_1m[ix0]);
        _complex_0(ladd_2m[ix0]);
        _complex_0(ladd_3m[ix0]);

        _complex_0(lvdd_m[ix0]);
        _complex_0(lvdd_1m[ix0]);
        _complex_0(lvdd_2m[ix0]);
        _complex_0(lvdd_3m[ix0]);

        _complex_0(ltdd_m[ix0]);
        _complex_0(ltdd_1m[ix0]);
        _complex_0(ltdd_2m[ix0]);
        _complex_0(ltdd_3m[ix0]);

        _complex_0(lttdd_m[ix0]);
        _complex_0(lttdd_1m[ix0]);
        _complex_0(lttdd_2m[ix0]);
        _complex_0(lttdd_3m[ix0]);
     
        for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
  	  _complex_0(la1_ij[ix0][s1][s2]);
  	  _complex_0(la2_ij[ix0][s1][s2]);
  	  _complex_0(la3_ij[ix0][s1][s2]);

  	  _complex_0(lv1_ij[ix0][s1][s2]);
  	  _complex_0(lv2_ij[ix0][s1][s2]);
  	  _complex_0(lv3_ij[ix0][s1][s2]);

  	  _complex_0(lt1_ij[ix0][s1][s2]);
  	  _complex_0(lt2_ij[ix0][s1][s2]);
  	  _complex_0(lt3_ij[ix0][s1][s2]);

  	  _complex_0(ltt1_ij[ix0][s1][s2]);
  	  _complex_0(ltt2_ij[ix0][s1][s2]);
  	  _complex_0(ltt3_ij[ix0][s1][s2]);
  	}
      }


  /********            lY_1           *******************/
    /*  Gamma matrix structure:    \csi^dag_j  g5 g1 Q- csi_i    */

    for(s1=0;s1<4*NF;s1++){
      for(s2=0;s2<4*NF;s2++){
  	/*g5 g1/i Q-*/
          /*Q-*/
          stmp1[0]=chi[s2];
        _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
        _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
        _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
        _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
        _spinor_mul_f(stmp1[0],.5,stmp1[0]);

  	/*g5g1/i*/
          stmp1[1].c[0]=stmp1[0].c[3];
          stmp1[1].c[1]=stmp1[0].c[2];
          stmp1[1].c[2]=stmp1[0].c[1];
          stmp1[1].c[3]=stmp1[0].c[0];
          _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
          _vector_mul_f(stmp1[1].c[1],-1.0,stmp1[1].c[1]);
          _vector_mul_f(stmp1[1].c[2],-1.0,stmp1[1].c[2]);
          _vector_mul_f(stmp1[1].c[3],-1.0,stmp1[1].c[3]);
        _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

      }
    }
    /**** construction of l1_ij with open indeces ij ******/

    for(s2=0;s2<4*NF;s2++){
      for(s1=0;s1<4*NF;s1++){

        for(ix0=0;ix0<T;ix0++){

  	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
  	      i=ipt(ix0,ix1,ix2,ix3);
  	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
  	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

  	      /*ladd_1m*/
  	      /*gamma_1/i*/
  	      stmp1[0].c[0]=sptr1[0]->c[3];
  	      stmp1[0].c[1]=sptr1[0]->c[2];
  	      stmp1[0].c[2]=sptr1[0]->c[1];
  	      stmp1[0].c[3]=sptr1[0]->c[0];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);


  	      /*lvdd_1m*/
  	      /*g5 gamma_1/i*/
  	      stmp1[0].c[0]=sptr1[0]->c[3];
  	      stmp1[0].c[1]=sptr1[0]->c[2];
  	      stmp1[0].c[2]=sptr1[0]->c[1];
  	      stmp1[0].c[3]=sptr1[0]->c[0];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
  	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
  	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*ltdd_1m*/
  	      /*g5 sig01*/
  	      stmp1[0].c[0]=sptr1[0]->c[1];
  	      stmp1[0].c[1]=sptr1[0]->c[0];
  	      stmp1[0].c[2]=sptr1[0]->c[3];
  	      stmp1[0].c[3]=sptr1[0]->c[2];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
  	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
  	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*lttdd_1m*/
  	      /*sig01*/
  	      stmp1[0].c[0]=sptr1[0]->c[1];
  	      stmp1[0].c[1]=sptr1[0]->c[0];
  	      stmp1[0].c[2]=sptr1[0]->c[3];
  	      stmp1[0].c[3]=sptr1[0]->c[2];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
  	    }
          // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
  	_complex_mulr(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
        }
      }
    }

    /***** Here we do the contractions of lY1_ij M_ij    ****/


    for(s2=0;s2<4*NF;s2++){
      for(s1=0;s1<4*NF;s1++){
        for(ix0=0;ix0<T;ix0++){
  	// Contract and accumulate
  	_complex_mul(temp_comp,la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(ladd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(lvdd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(ltdd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(lttdd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }


  /********            lY_2           *******************/
    /*  Gamma matrix structure:    \csi^dag_j  g5 g2 Q- csi_i    */

    for(s1=0;s1<4*NF;s1++){
      for(s2=0;s2<4*NF;s2++){
  	/*g5 g2 Q-*/
          /*Q-*/
          stmp1[0]=chi[s2];
        _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
        _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
        _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
        _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
        _spinor_mul_f(stmp1[0],.5,stmp1[0]);

  	/*g5g2*/
          stmp1[1].c[0]=stmp1[0].c[3];
          stmp1[1].c[1]=stmp1[0].c[2];
          stmp1[1].c[2]=stmp1[0].c[1];
          stmp1[1].c[3]=stmp1[0].c[0];
          _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
          _vector_mul_f(stmp1[1].c[2],-1.0,stmp1[1].c[2]);
        _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

      }
    }
    /**** construction of l2_ij with open indeces ij ******/

    for(s2=0;s2<4*NF;s2++){
      for(s1=0;s1<4*NF;s1++){

        for(ix0=0;ix0<T;ix0++){
  	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
  	      i=ipt(ix0,ix1,ix2,ix3);
  	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
  	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

  	      /*ladd_2m*/
  	      /*gamma_2*/

  	      stmp1[0].c[0]=sptr1[0]->c[3];
  	      stmp1[0].c[1]=sptr1[0]->c[2];
  	      stmp1[0].c[2]=sptr1[0]->c[1];
  	      stmp1[0].c[3]=sptr1[0]->c[0];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*lvdd_2m*/
  	      /*g5 gamma_2*/
  	      stmp1[0].c[0]=sptr1[0]->c[3];
  	      stmp1[0].c[1]=sptr1[0]->c[2];
  	      stmp1[0].c[2]=sptr1[0]->c[1];
  	      stmp1[0].c[3]=sptr1[0]->c[0];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*ltdd_2m*/
  	      /*g5 sig02*/
  	      stmp1[0].c[0]=sptr1[0]->c[1];
  	      stmp1[0].c[1]=sptr1[0]->c[0];
  	      stmp1[0].c[2]=sptr1[0]->c[3];
  	      stmp1[0].c[3]=sptr1[0]->c[2];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*lttdd_2m*/
  	      /*sig02*/
  	      stmp1[0].c[0]=sptr1[0]->c[1];
  	      stmp1[0].c[1]=sptr1[0]->c[0];
  	      stmp1[0].c[2]=sptr1[0]->c[3];
  	      stmp1[0].c[3]=sptr1[0]->c[2];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);


  	    }
          // Normalization with -1.0/(6.0*GLB_VOL3)
  	_complex_mulr(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
        }
      }
    }

    /***** Here we do the contractions of lY2_ij M_ij    ****/


    for(s2=0;s2<4*NF;s2++){
      for(s1=0;s1<4*NF;s1++){
        for(ix0=0;ix0<T;ix0++){
  	// Contract and accumulate
  	_complex_mul(temp_comp,la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(ladd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(lvdd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(ltdd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(lttdd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }




  /********            lY_3           *******************/
    /*  Gamma matrix structure:    \csi^dag_j  g5 g3 Q- csi_i    */

    for(s1=0;s1<4*NF;s1++){
      for(s2=0;s2<4*NF;s2++){
  	/*g5 g3/i Q-*/
          /*Q-*/
          stmp1[0]=chi[s2];
        _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
        _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
        _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
        _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
        _spinor_mul_f(stmp1[0],.5,stmp1[0]);

  	/*g5g3/i*/

          stmp1[1].c[0]=stmp1[0].c[2];
          stmp1[1].c[1]=stmp1[0].c[3];
          stmp1[1].c[2]=stmp1[0].c[0];
          stmp1[1].c[3]=stmp1[0].c[1];
          _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
          _vector_mul_f(stmp1[1].c[2],-1.0,stmp1[1].c[2]);
        _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

      }
    }
    /**** construction of l3_ij with open indeces ij ******/

    for(s2=0;s2<4*NF;s2++){
      for(s1=0;s1<4*NF;s1++){

        for(ix0=0;ix0<T;ix0++){
  	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
  	      i=ipt(ix0,ix1,ix2,ix3);
  	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
  	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

  	      /*ladd_3m*/
  	      /*gamma_3/i*/

  	      stmp1[0].c[0]=sptr1[0]->c[2];
  	      stmp1[0].c[1]=sptr1[0]->c[3];
  	      stmp1[0].c[2]=sptr1[0]->c[0];
  	      stmp1[0].c[3]=sptr1[0]->c[1];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);



  	      /*lvdd_3m*/
  	      /*g5 gamma_3/i*/
  	      stmp1[0].c[0]=sptr1[0]->c[2];
  	      stmp1[0].c[1]=sptr1[0]->c[3];
  	      stmp1[0].c[2]=sptr1[0]->c[0];
  	      stmp1[0].c[3]=sptr1[0]->c[1];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*ltdd_3m*/
  	      /*g5 sig03*/
  	      stmp1[0].c[0]=sptr1[0]->c[0];
  	      stmp1[0].c[1]=sptr1[0]->c[1];
  	      stmp1[0].c[2]=sptr1[0]->c[2];
  	      stmp1[0].c[3]=sptr1[0]->c[3];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

  	      /*lttdd_3m*/
  	      /*sig03*/
    	      stmp1[0].c[0]=sptr1[0]->c[0];
  	      stmp1[0].c[1]=sptr1[0]->c[1];
  	      stmp1[0].c[2]=sptr1[0]->c[2];
  	      stmp1[0].c[3]=sptr1[0]->c[3];
  	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
  	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

  	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
  	      //  acumulate for the spatial average
  	      _complex_add_assign(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);



  	    }
          // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
  	_complex_mulr(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
  	_complex_mulr(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
        }
      }
    }

    /***** Here we do the contractions of lY3_ij M_ij    ****/


    for(s2=0;s2<4*NF;s2++){
      for(s1=0;s1<4*NF;s1++){
        for(ix0=0;ix0<T;ix0++){
  	// Contract and accumulate
  	_complex_mul(temp_comp,la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(ladd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(lvdd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(ltdd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

  	_complex_mul(temp_comp,ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
  	_complex_add_assign(lttdd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }


  /*******   Sum over k ************/
     for(ix0=0;ix0<T;ix0++){
  	_complex_add_assign(ladd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ladd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(ladd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ladd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(ladd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ladd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);

  	_complex_add_assign(lvdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvdd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(lvdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvdd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(lvdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvdd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);

  	_complex_add_assign(ltdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltdd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(ltdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltdd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(ltdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltdd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);

  	_complex_add_assign(lttdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttdd_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(lttdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttdd_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  	_complex_add_assign(lttdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttdd_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
    }



    /*****   Average over processors in temporal direction *****/
    global_sum((double*)ladd_m,2*GLB_T);
    global_sum((double*)lvdd_m,2*GLB_T);
    global_sum((double*)ltdd_m,2*GLB_T);
    global_sum((double*)lttdd_m,2*GLB_T);

    /************************************ END of set lY_dd_-*********************/
  

  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," la_dd_-[%d] = %.10e,%.10e\n",ix0,ladd_m[ix0].re,ladd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lv_dd_-[%d] = %.10e,%.10e\n",ix0,lvdd_m[ix0].re,lvdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lt_dd_-[%d] = %.10e,%.10e\n",ix0,ltdd_m[ix0].re,ltdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ltt_dd_-[%d] = %.10e,%.10e\n",ix0,lttdd_m[ix0].re,lttdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}




void rotated_lXudm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *laud_m=corr_mem->g1, *lvud_m=corr_mem->g2, *ltud_m=corr_mem->g3, *lttud_m=corr_mem->g4, **M=corr_mem->M;

  complex *laud_1m=corr_mem->l11, *lvud_1m=corr_mem->l21, *ltud_1m=corr_mem->l31, *lttud_1m=corr_mem->l41;
  complex *laud_2m=corr_mem->l12, *lvud_2m=corr_mem->l22, *ltud_2m=corr_mem->l32, *lttud_2m=corr_mem->l42;
  complex *laud_3m=corr_mem->l13, *lvud_3m=corr_mem->l23, *ltud_3m=corr_mem->l33, *lttud_3m=corr_mem->l43;

  complex ***la1_ij=corr_mem->l11_ij, ***la2_ij=corr_mem->l12_ij, ***la3_ij=corr_mem->l13_ij;
  complex ***lv1_ij=corr_mem->l21_ij, ***lv2_ij=corr_mem->l22_ij, ***lv3_ij=corr_mem->l23_ij;
  complex ***lt1_ij=corr_mem->l31_ij, ***lt2_ij=corr_mem->l32_ij, ***lt3_ij=corr_mem->l33_ij;
  complex ***ltt1_ij=corr_mem->g1_ij,  ***ltt2_ij=corr_mem->g2_ij,  ***ltt3_ij=corr_mem->g3_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp,temp_comp2;





  /**************************************************

    set:   laud_m, lpud_m, lsud_m, lvud_m

laud- = - i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   gk     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q+ csi

    lauu_1- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g1     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q+ csi
    lauu_2- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g2     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q+ csi
    lauu_3- = - i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g3     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q+ csi

lvud- = i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5gk     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q+ csi

    lvuu_1- = i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5g1     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q+ csi
    lvuu_2- = i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5g2     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q+ csi
    lvuu_3- = i 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5g3     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q+ csi

ltud- = - 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig0k     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q+ csi

    ltuu_1- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig01     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q+ csi
    ltuu_2- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig02     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q+ csi
    ltuu_3- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   g5 sig03     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q+ csi


lttud- = - 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)    sig0k     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q+ csi

    lttuu_1- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   sig01     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q+ csi
    lttuu_2- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   sig02     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q+ csi
    lttuu_3- = - 1/6 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;x)   sig03     Hdd^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q+ csi

															(sig0k=ig0gk)
  **************************************************/

  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(laud_m[ix0]);
      _complex_0(laud_1m[ix0]);
      _complex_0(laud_2m[ix0]);
      _complex_0(laud_3m[ix0]);

      _complex_0(lvud_m[ix0]);
      _complex_0(lvud_1m[ix0]);
      _complex_0(lvud_2m[ix0]);
      _complex_0(lvud_3m[ix0]);

      _complex_0(ltud_m[ix0]);
      _complex_0(ltud_1m[ix0]);
      _complex_0(ltud_2m[ix0]);
      _complex_0(ltud_3m[ix0]);

      _complex_0(lttud_m[ix0]);
      _complex_0(lttud_1m[ix0]);
      _complex_0(lttud_2m[ix0]);
      _complex_0(lttud_3m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(la1_ij[ix0][s1][s2]);
	  _complex_0(la2_ij[ix0][s1][s2]);
	  _complex_0(la3_ij[ix0][s1][s2]);

	  _complex_0(lv1_ij[ix0][s1][s2]);
	  _complex_0(lv2_ij[ix0][s1][s2]);
	  _complex_0(lv3_ij[ix0][s1][s2]);

	  _complex_0(lt1_ij[ix0][s1][s2]);
	  _complex_0(lt2_ij[ix0][s1][s2]);
	  _complex_0(lt3_ij[ix0][s1][s2]);

	  _complex_0(ltt1_ij[ix0][s1][s2]);
	  _complex_0(ltt2_ij[ix0][s1][s2]);
	  _complex_0(ltt3_ij[ix0][s1][s2]);
	}
    }


/********            lY_1           *******************/
  /*  Gamma matrix structure:    \csi^dag_j g1 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g1/i Q+*/
        /*Q+*/
        stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g1/i*/
        stmp1[1].c[0]=stmp1[0].c[3];
        stmp1[1].c[1]=stmp1[0].c[2];
        stmp1[1].c[2]=stmp1[0].c[1];
        stmp1[1].c[3]=stmp1[0].c[0];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[1],-1.0,stmp1[1].c[1]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l1_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	      /*laud_1m*/
	      /*gamma_1/i*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);


	      /*lvuu_1m*/
	      /*g5 gamma_1/i*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*ltuu_1m*/
	      /*g5 sig01*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lttuu_1m*/
	      /*sig01*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
	    }
        // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
	_complex_mulr(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }

  /***** Here we do the contractions of lY1_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(laud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}


/********            lY_2           *******************/
  /*  Gamma matrix structure:    \csi^dag_j  g2 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g2 Q+*/
        /*Q+*/
        stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g2*/
        stmp1[1].c[0]=stmp1[0].c[3];
        stmp1[1].c[1]=stmp1[0].c[2];
        stmp1[1].c[2]=stmp1[0].c[1];
        stmp1[1].c[3]=stmp1[0].c[0];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[3],-1.0,stmp1[1].c[3]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l2_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	      /*laud_2m*/
	      /*gamma_2*/

	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lvud_2m*/
	      /*g5 gamma_2*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*ltud_2m*/
	      /*g5 sig02*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lttud_2m*/
	      /*sig02*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);


	    }
        // Normalization with -1.0/(6.0*GLB_VOL3)
	_complex_mulr(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }

  /***** Here we do the contractions of lY2_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(laud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}




/********            lY_3           *******************/
  /*  Gamma matrix structure:    \csi^dag_j g3 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g3/i Q+*/
        /*Q+*/
        stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g3/i*/

        stmp1[1].c[0]=stmp1[0].c[2];
        stmp1[1].c[1]=stmp1[0].c[3];
        stmp1[1].c[2]=stmp1[0].c[0];
        stmp1[1].c[3]=stmp1[0].c[1];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[3],-1.0,stmp1[1].c[3]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l3_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){

	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	      /*lauu_3m*/
	      /*gamma_3/i*/

	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);



	      /*lvud_3m*/
	      /*g5 gamma_3/i*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*ltud_3m*/
	      /*g5 sig03*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lttud_3m*/
	      /*sig03*/
  	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);



	    }
        // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
	_complex_mulr(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }

  /***** Here we do the contractions of lY3_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(laud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}


/*******   Sum over k and multiply by i************/

   _complex_i(temp_comp);

	
   for(ix0=0;ix0<T;ix0++){
	_complex_add_assign(laud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],laud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(laud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],laud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(laud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],laud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,laud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	laud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;

	_complex_add_assign(lvud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lvud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lvud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,lvud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	lvud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;

	_complex_add_assign(ltud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ltud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ltud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,ltud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	ltud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;

	_complex_add_assign(lttud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttud_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lttud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttud_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lttud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttud_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,lttud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	lttud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;
  }


    /*****   Average over processors in temporal direction *****/
  global_sum((double*)laud_m,2*GLB_T);
  global_sum((double*)lvud_m,2*GLB_T);
  global_sum((double*)ltud_m,2*GLB_T);
  global_sum((double*)lttud_m,2*GLB_T);



  /************************************ END of set lY_ud_-*********************/



  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," la_ud_-[%d] = %.10e,%.10e\n",ix0,laud_m[ix0].re,laud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lv_ud_-[%d] = %.10e,%.10e\n",ix0,lvud_m[ix0].re,lvud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lt_ud_-[%d] = %.10e,%.10e\n",ix0,ltud_m[ix0].re,ltud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ltt_ud_-[%d] = %.10e,%.10e\n",ix0,lttud_m[ix0].re,lttud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}




void rotated_lXdum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *ladu_m=corr_mem->g1, *lvdu_m=corr_mem->g2, *ltdu_m=corr_mem->g3, *lttdu_m=corr_mem->g4, **M=corr_mem->M;

  complex *ladu_1m=corr_mem->l11, *lvdu_1m=corr_mem->l21, *ltdu_1m=corr_mem->l31, *lttdu_1m=corr_mem->l41;
  complex *ladu_2m=corr_mem->l12, *lvdu_2m=corr_mem->l22, *ltdu_2m=corr_mem->l32, *lttdu_2m=corr_mem->l42;
  complex *ladu_3m=corr_mem->l13, *lvdu_3m=corr_mem->l23, *ltdu_3m=corr_mem->l33, *lttdu_3m=corr_mem->l43;

  complex ***la1_ij=corr_mem->l11_ij, ***la2_ij=corr_mem->l12_ij, ***la3_ij=corr_mem->l13_ij;
  complex ***lv1_ij=corr_mem->l21_ij, ***lv2_ij=corr_mem->l22_ij, ***lv3_ij=corr_mem->l23_ij;
  complex ***lt1_ij=corr_mem->l31_ij, ***lt2_ij=corr_mem->l32_ij, ***lt3_ij=corr_mem->l33_ij;
  complex ***ltt1_ij=corr_mem->g1_ij,  ***ltt2_ij=corr_mem->g2_ij,  ***ltt3_ij=corr_mem->g3_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp,temp_comp2;



  /**************************************************
ladu- = i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   gk     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q- csi

    ladu_1- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g1     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q- csi
    ladu_2- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g2     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q- csi
    ladu_3- = i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g3     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q- csi

lvdu- = -i 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5gk     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q- csi

    lvdu_1- = -i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5g1     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q- csi
    lvdu_2- = -i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5g2     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q- csi
    lvdu_3- = -i 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5g3     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q- csi

ltdu- = 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig0k     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q- csi

    ltdu_1- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig01     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q- csi
    ltdu_2- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig02     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q- csi
    ltdu_3- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   g5 sig03     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q- csi

lttdu- = 1/6 sum_k=1^3 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)    sig0k     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag gk Q- csi

    lttdu_1- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   sig01     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g1 Q- csi
    lttdu_2- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   sig02     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g2 Q- csi
    lttdu_3- = 1/6 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;x)   sig03     Huu^(-1)(x:2,y) U0(1,y)^dag psi       psi^dag g3 Q- csi

															(sig0k=ig0gk)
  **************************************************/

  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(ladu_m[ix0]);
      _complex_0(ladu_1m[ix0]);
      _complex_0(ladu_2m[ix0]);
      _complex_0(ladu_3m[ix0]);

      _complex_0(lvdu_m[ix0]);
      _complex_0(lvdu_1m[ix0]);
      _complex_0(lvdu_2m[ix0]);
      _complex_0(lvdu_3m[ix0]);

      _complex_0(ltdu_m[ix0]);
      _complex_0(ltdu_1m[ix0]);
      _complex_0(ltdu_2m[ix0]);
      _complex_0(ltdu_3m[ix0]);

      _complex_0(lttdu_m[ix0]);
      _complex_0(lttdu_1m[ix0]);
      _complex_0(lttdu_2m[ix0]);
      _complex_0(lttdu_3m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(la1_ij[ix0][s1][s2]);
	  _complex_0(la2_ij[ix0][s1][s2]);
	  _complex_0(la3_ij[ix0][s1][s2]);

	  _complex_0(lv1_ij[ix0][s1][s2]);
	  _complex_0(lv2_ij[ix0][s1][s2]);
	  _complex_0(lv3_ij[ix0][s1][s2]);

	  _complex_0(lt1_ij[ix0][s1][s2]);
	  _complex_0(lt2_ij[ix0][s1][s2]);
	  _complex_0(lt3_ij[ix0][s1][s2]);

	  _complex_0(ltt1_ij[ix0][s1][s2]);
	  _complex_0(ltt2_ij[ix0][s1][s2]);
	  _complex_0(ltt3_ij[ix0][s1][s2]);
	}
    }


/********            lY_1           *******************/
  /*  Gamma matrix structure:    \csi^dag_j g1 Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g1/i Q-*/
        /*Q-*/
        stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g1/i*/
        stmp1[1].c[0]=stmp1[0].c[3];
        stmp1[1].c[1]=stmp1[0].c[2];
        stmp1[1].c[2]=stmp1[0].c[1];
        stmp1[1].c[3]=stmp1[0].c[0];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[1],-1.0,stmp1[1].c[1]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l1_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

	      /*ladu_1m*/
	      /*gamma_1/i*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);


	      /*lvdu_1m*/
	      /*g5 gamma_1/i*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*ltdu_1m*/
	      /*g5 sig01*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lttdu_1m*/
	      /*sig01*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[1],-1.0,stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
	    }
        // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
	_complex_mulr(la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }

  /***** Here we do the contractions of lY1_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ladu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvdu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltdu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttdu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}


/********            lY_2           *******************/
  /*  Gamma matrix structure:    \csi^dag_j  g2 Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g2 Q-*/
        /*Q-*/
        stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g2*/
        stmp1[1].c[0]=stmp1[0].c[3];
        stmp1[1].c[1]=stmp1[0].c[2];
        stmp1[1].c[2]=stmp1[0].c[1];
        stmp1[1].c[3]=stmp1[0].c[0];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[3],-1.0,stmp1[1].c[3]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l2_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

	      /*ladu_2m*/
	      /*gamma_2*/

	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lvdu_2m*/
	      /*g5 gamma_2*/
	      stmp1[0].c[0]=sptr1[0]->c[3];
	      stmp1[0].c[1]=sptr1[0]->c[2];
	      stmp1[0].c[2]=sptr1[0]->c[1];
	      stmp1[0].c[3]=sptr1[0]->c[0];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*ltdu_2m*/
	      /*g5 sig02*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lttdu_2m*/
	      /*sig02*/
	      stmp1[0].c[0]=sptr1[0]->c[1];
	      stmp1[0].c[1]=sptr1[0]->c[0];
	      stmp1[0].c[2]=sptr1[0]->c[3];
	      stmp1[0].c[3]=sptr1[0]->c[2];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);


	    }
        // Normalization with -1.0/(6.0*GLB_VOL3)
	_complex_mulr(la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }

  /***** Here we do the contractions of lY2_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ladu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvdu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltdu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttdu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}




/********            lY_3           *******************/
  /*  Gamma matrix structure:    \csi^dag_j g3 Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
	/*g5 g3/i Q-*/
        /*Q-*/
        stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

	/*g5g3/i*/

        stmp1[1].c[0]=stmp1[0].c[2];
        stmp1[1].c[1]=stmp1[0].c[3];
        stmp1[1].c[2]=stmp1[0].c[0];
        stmp1[1].c[3]=stmp1[0].c[1];
        _vector_mul_f(stmp1[1].c[0],-1.0,stmp1[1].c[0]);
        _vector_mul_f(stmp1[1].c[3],-1.0,stmp1[1].c[3]);
      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of l3_ij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){

	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

	      /*ladu_3m*/
	      /*gamma_3/i*/

	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);



	      /*lvdu_3m*/
	      /*g5 gamma_3/i*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*ltdu_3m*/
	      /*g5 sig03*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[2],-1.0,stmp1[0].c[2]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*lttdu_3m*/
	      /*sig03*/
  	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_f(stmp1[0].c[0],-1.0,stmp1[0].c[0]);
	      _vector_mul_f(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);



	    }
        // Normalization with 1.0/(6.0*GLB_VOL3), (Remember that the minus signs pf lY1 and lY3 cancels out from the imaginary part of the gamma matrix
	_complex_mulr(la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],1.0/(6.0*GLB_VOL3),la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-1.0/(6.0*GLB_VOL3),ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }

  /***** Here we do the contractions of lY3_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,la3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ladu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lv3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lvdu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,lt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(ltdu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,ltt3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(lttdu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
    }
  }
}


/*******   Sum over k and multiply by i************/

   _complex_i(temp_comp);

	
   for(ix0=0;ix0<T;ix0++){
	_complex_add_assign(ladu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ladu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ladu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ladu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ladu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ladu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,ladu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	ladu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;

	_complex_add_assign(lvdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvdu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lvdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvdu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lvdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lvdu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,lvdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	lvdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;

	_complex_add_assign(ltdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltdu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ltdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltdu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(ltdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],ltdu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,ltdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	ltdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;

	_complex_add_assign(lttdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttdu_1m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lttdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttdu_2m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(lttdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],lttdu_3m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_mul(temp_comp2,temp_comp,lttdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	lttdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]=temp_comp2;
  }


    /*****   Average over processors in temporal direction *****/
  global_sum((double*)ladu_m,2*GLB_T);
  global_sum((double*)lvdu_m,2*GLB_T);
  global_sum((double*)ltdu_m,2*GLB_T);
  global_sum((double*)lttdu_m,2*GLB_T);



  /************************************ END of set lY_du_-*********************/




  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," la_du_-[%d] = %.10e,%.10e\n",ix0,ladu_m[ix0].re,ladu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lv_du_-[%d] = %.10e,%.10e\n",ix0,lvdu_m[ix0].re,lvdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," lt_du_-[%d] = %.10e,%.10e\n",ix0,ltdu_m[ix0].re,ltdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ltt_du_-[%d] = %.10e,%.10e\n",ix0,lttdu_m[ix0].re,lttdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}

#endif
