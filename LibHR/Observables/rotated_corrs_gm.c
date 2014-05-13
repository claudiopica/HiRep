#ifdef ROTATED_SF

#include "global.h"
#include "communications.h"
#include "observables.h"
#include "logger.h"




/***************************************************************

	Begin the computation of the g_X_f1f2 correlation functions with the inverted boundary projectors



*******************************************************************************************************/

void rotated_gXuum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *gsuu_m=corr_mem->g1, *gpuu_m=corr_mem->g2, *gvuu_m=corr_mem->g3, *gauu_m=corr_mem->g4, **M=corr_mem->M;

  complex ***gs_ij=corr_mem->g1_ij,  ***gp_ij=corr_mem->g2_ij,  ***gv_ij=corr_mem->g3_ij,  ***ga_ij=corr_mem->g4_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;

  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;


  /**************************************************
    set:   gauu_p, gpuu_p, gsuu_p, gvuu_p

gpuu- = -  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       II      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       

gauu- =  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       

gsuu- = -  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       

gvuu- = -  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gauu_m[ix0]);
      _complex_0(gpuu_m[ix0]);
      _complex_0(gsuu_m[ix0]);
      _complex_0(gvuu_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(ga_ij[ix0][s1][s2]);
	  _complex_0(gp_ij[ix0][s1][s2]);
	  _complex_0(gv_ij[ix0][s1][s2]);
	  _complex_0(gs_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  g0 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q+*/
      stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      /*gamma_0*/
      stmp1[1].c[0]=stmp1[0].c[2];
      stmp1[1].c[1]=stmp1[0].c[3];
      stmp1[1].c[2]=stmp1[0].c[0];
      stmp1[1].c[3]=stmp1[0].c[1];
      _spinor_mul_f(stmp1[1],-1.0,stmp1[1]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
 



  /**** construction of gij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

	      /*gsuu*/						
	      /*gamma_5*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_g(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_g(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*gvuu*/
	      /*gamma_5*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_minus_g(stmp1[0].c[0],stmp1[0].c[0]);
	      _vector_minus_g(stmp1[0].c[1],stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],.5/GLB_VOL3,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gauu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gpuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gsuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gauu_m,2*GLB_T);
  global_sum((double*)gpuu_m,2*GLB_T);
  global_sum((double*)gsuu_m,2*GLB_T);
  global_sum((double*)gvuu_m,2*GLB_T);

  /************************************ END of set gX_uu_-*********************/
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_uu_-[%d] = %.10e,%.10e\n",ix0,gauu_m[ix0].re,gauu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_uu_-[%d] = %.10e,%.10e\n",ix0,gvuu_m[ix0].re,gvuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_uu_-[%d] = %.10e,%.10e\n",ix0,gpuu_m[ix0].re,gpuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_uu_-[%d] = %.10e,%.10e\n",ix0,gsuu_m[ix0].re,gsuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
}


void rotated_gXddm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *gsdd_m=corr_mem->g1, *gpdd_m=corr_mem->g2, *gvdd_m=corr_mem->g3, *gadd_m=corr_mem->g4, **M=corr_mem->M;

  complex ***gs_ij=corr_mem->g1_ij,  ***gp_ij=corr_mem->g2_ij,  ***gv_ij=corr_mem->g3_ij,  ***ga_ij=corr_mem->g4_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;


  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;




  /**************************************************
    set:   gadd_m, gpdd_m, gsdd_m, gvdd_m

gpdd- = -  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       II      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       

gadd- =  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       

gsdd- = -  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       

gvdd- = -  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gadd_m[ix0]);
      _complex_0(gpdd_m[ix0]);
      _complex_0(gsdd_m[ix0]);
      _complex_0(gvdd_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(ga_ij[ix0][s1][s2]);
	  _complex_0(gp_ij[ix0][s1][s2]);
	  _complex_0(gv_ij[ix0][s1][s2]);
	  _complex_0(gs_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  g0 Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q-*/
      stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      /*gamma_0*/
      stmp1[1].c[0]=stmp1[0].c[2];
      stmp1[1].c[1]=stmp1[0].c[3];
      stmp1[1].c[2]=stmp1[0].c[0];
      stmp1[1].c[3]=stmp1[0].c[1];
      _spinor_mul_f(stmp1[1],-1.0,stmp1[1]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
 


  /**** construction of gij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

	      /*gsuu*/						
	      /*gamma_5*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_g(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_g(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*gvuu*/
	      /*gamma_5*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_minus_g(stmp1[0].c[0],stmp1[0].c[0]);
	      _vector_minus_g(stmp1[0].c[1],stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],.5/GLB_VOL3,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gadd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gpdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gsdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gadd_m,2*GLB_T);
  global_sum((double*)gpdd_m,2*GLB_T);
  global_sum((double*)gsdd_m,2*GLB_T);
  global_sum((double*)gvdd_m,2*GLB_T);

  /************************************ END of set gX_dd_-*********************/
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_dd_-[%d] = %.10e,%.10e\n",ix0,gadd_m[ix0].re,gadd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_dd_-[%d] = %.10e,%.10e\n",ix0,gvdd_m[ix0].re,gvdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_dd_-[%d] = %.10e,%.10e\n",ix0,gpdd_m[ix0].re,gpdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_dd_-[%d] = %.10e,%.10e\n",ix0,gsdd_m[ix0].re,gsdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
}


void rotated_gXudm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *gsud_m=corr_mem->g1, *gpud_m=corr_mem->g2, *gvud_m=corr_mem->g3, *gaud_m=corr_mem->g4, **M=corr_mem->M;

  complex ***gs_ij=corr_mem->g1_ij,  ***gp_ij=corr_mem->g2_ij,  ***gv_ij=corr_mem->g3_ij,  ***ga_ij=corr_mem->g4_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;


  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;




  /**************************************************
    set:   gaud_m, gpud_m, gsud_m, gvud_m

gpud- =  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       II      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       

gaud- = - 1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       

gsud- =   1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       

gvud- =   1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {

      _complex_0(gsud_m[ix0]);
      _complex_0(gvud_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(ga_ij[ix0][s1][s2]);
	  _complex_0(gp_ij[ix0][s1][s2]);
	  _complex_0(gv_ij[ix0][s1][s2]);
	  _complex_0(gs_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q+*/
      stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[0]);

    }
  }
 


  /**** construction of gij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

	      /*gsuu*/						
	      /*gamma_5*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_g(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_g(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*gvuu*/
	      /*gamma_5*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_minus_g(stmp1[0].c[0],stmp1[0].c[0]);
	      _vector_minus_g(stmp1[0].c[1],stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gaud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gpud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gsud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gaud_m,2*GLB_T);
  global_sum((double*)gpud_m,2*GLB_T);
  global_sum((double*)gsud_m,2*GLB_T);
  global_sum((double*)gvud_m,2*GLB_T);

  /************************************ END of set gX_ud_-*********************/
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_ud_-[%d] = %.10e,%.10e\n",ix0,gaud_m[ix0].re,gaud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_ud_-[%d] = %.10e,%.10e\n",ix0,gvud_m[ix0].re,gvud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_ud_-[%d] = %.10e,%.10e\n",ix0,gpud_m[ix0].re,gpud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_ud_-[%d] = %.10e,%.10e\n",ix0,gsud_m[ix0].re,gsud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
}



void rotated_gXdum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {

  complex *gsdu_m=corr_mem->g1, *gpdu_m=corr_mem->g2, *gvdu_m=corr_mem->g3, *gadu_m=corr_mem->g4, **M=corr_mem->M;

  complex ***gs_ij=corr_mem->g1_ij,  ***gp_ij=corr_mem->g2_ij,  ***gv_ij=corr_mem->g3_ij,  ***ga_ij=corr_mem->g4_ij;

  int i,ix0,ix1,ix2,ix3,s1,s2;


  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;




  /**************************************************
    set:   gadu_m, gpdu_m, gsdu_m, gvdu_m

gpud- =  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       II      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       

gaud- = - 1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       

gsud- =   1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       

gvud- =   1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gadu_m[ix0]);
      _complex_0(gpdu_m[ix0]);
      _complex_0(gsdu_m[ix0]);
      _complex_0(gvdu_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(ga_ij[ix0][s1][s2]);
	  _complex_0(gp_ij[ix0][s1][s2]);
	  _complex_0(gv_ij[ix0][s1][s2]);
	  _complex_0(gs_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q-*/
      stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[0]);

    }
  }
 


  /**** construction of gij with open indeces ij ******/

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

	      /*gsuu*/						
	      /*gamma_5*/
	      stmp1[0].c[0]=sptr1[0]->c[0];
	      stmp1[0].c[1]=sptr1[0]->c[1];
	      stmp1[0].c[2]=sptr1[0]->c[2];
	      stmp1[0].c[3]=sptr1[0]->c[3];
	      _vector_mul_g(stmp1[0].c[2],-1.0,stmp1[0].c[2]);
	      _vector_mul_g(stmp1[0].c[3],-1.0,stmp1[0].c[3]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

	      /*gvuu*/
	      /*gamma_5*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _vector_minus_g(stmp1[0].c[0],stmp1[0].c[0]);
	      _vector_minus_g(stmp1[0].c[1],stmp1[0].c[1]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,ga_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gadu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gp_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gpdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gs_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gsdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gv_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gadu_m,2*GLB_T);
  global_sum((double*)gpdu_m,2*GLB_T);
  global_sum((double*)gsdu_m,2*GLB_T);
  global_sum((double*)gvdu_m,2*GLB_T);

  /************************************ END of set gX_du_-*********************/

  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_du_-[%d] = %.10e,%.10e\n",ix0,gadu_m[ix0].re,gadu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_du_-[%d] = %.10e,%.10e\n",ix0,gvdu_m[ix0].re,gvdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_du_-[%d] = %.10e,%.10e\n",ix0,gpdu_m[ix0].re,gpdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_du_-[%d] = %.10e,%.10e\n",ix0,gsdu_m[ix0].re,gsdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
}





void rotated_gvtuum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex *gvtuu_m=corr_mem->g1, *gvt1uu_m=corr_mem->g2, *gvt2uu_m=corr_mem->g3, **M=corr_mem->M;

  complex ***gvt1_ij=corr_mem->g1_ij,  ***gvt2_ij=corr_mem->g2_ij;

  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;


  
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;


/****************    Computation of gvt_uu_-   

gvtuu- = - gvt1uu- + gvt2uu-


    gvt1uu- = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_uu        g5 P- U(x)    H(x+a0;a,y)_uu U^dag(0,y)psi      psi^dag g0 Q+ csi)

    gvt2uu- = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_uu     g5 P+ U^dag(x)     H(x;a,y)_uu U^dag(0,y)psi    psi^dag g0 Q+ csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gvt1uu_m[ix0]);
      _complex_0(gvt2uu_m[ix0]);
      _complex_0(gvtuu_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(gvt1_ij[ix0][s1][s2]);
	  _complex_0(gvt2_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  g0 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q+*/
      stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      /*gamma_0*/
      stmp1[1].c[0]=stmp1[0].c[2];
      stmp1[1].c[1]=stmp1[0].c[3];
      stmp1[1].c[2]=stmp1[0].c[0];
      stmp1[1].c[3]=stmp1[0].c[1];
      _spinor_mul_f(stmp1[1],-1.0,stmp1[1]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of gij with open indeces ij ******/

  /****  gvt1uu_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i2);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	     /****  U0(x) prop_uu  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
                _suNf_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]);
             } 

	    stmp1[1]=stmp1[0];
 	    /****  P-  *****/
	    _vector_add_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_add_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_add_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_add_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2uu_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i2);

	     /****  U0(x) prop_uu  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
               _suNf_inverse_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]); 
             }

	    stmp1[1]=stmp1[0];   
	     /****  P+  *****/
	    _vector_sub_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_sub_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_sub_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_sub_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T-1;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt1uu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt2uu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

/*******   Sum over k ************/
   for(ix0=0;ix0<T-1;ix0++){
	_complex_sub_assign(gvtuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt1uu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(gvtuu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt2uu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
   }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gvtuu_m,2*GLB_T);
  /************************************ END of set gvt_uu_-*********************/
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_uu_-[%d] = %.10e,%.10e\n",ix0,gvtuu_m[ix0].re,gvtuu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}



void rotated_gvtddm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex *gvtdd_m=corr_mem->g1, *gvt1dd_m=corr_mem->g2, *gvt2dd_m=corr_mem->g3, **M=corr_mem->M;

  complex ***gvt1_ij=corr_mem->g1_ij,  ***gvt2_ij=corr_mem->g2_ij;

  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;



/****************    Computation of gvt_dd_- 

gvtdd- = - gvt1dd- + gvt2dd-


    gvt1dd- = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_dd        g5 P- U(x)    H(x+a0;a,y)_dd U^dag(0,y)psi      psi^dag g0 Q- csi)

    gvt2dd- = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_dd     g5 P+ U^dag(x)     H(x;a,y)_dd U^dag(0,y)psi    psi^dag g0 Q- csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gvt1dd_m[ix0]);
      _complex_0(gvt2dd_m[ix0]);
      _complex_0(gvtdd_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(gvt1_ij[ix0][s1][s2]);
	  _complex_0(gvt2_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  g0 Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q-*/
      stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      /*gamma_0*/
      stmp1[1].c[0]=stmp1[0].c[2];
      stmp1[1].c[1]=stmp1[0].c[3];
      stmp1[1].c[2]=stmp1[0].c[0];
      stmp1[1].c[3]=stmp1[0].c[1];
      _spinor_mul_f(stmp1[1],-1.0,stmp1[1]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }
  /**** construction of gij with open indeces ij ******/

  /****  gvt1dd_- *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i2);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

	     /****  U0(x) prop_dd  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
                _suNf_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]);
             } 

	    stmp1[1]=stmp1[0];
 	    /****  P-  *****/
	    _vector_add_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_add_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_add_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_add_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2dd_- *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i2);

	     /****  U0(x) prop_dd  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
               _suNf_inverse_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]); 
             }

	    stmp1[1]=stmp1[0];   
	     /****  P+  *****/
	    _vector_sub_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_sub_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_sub_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_sub_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T-1;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt1dd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt2dd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

/*******   Sum over k ************/
   for(ix0=0;ix0<T-1;ix0++){
	_complex_sub_assign(gvtdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt1dd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(gvtdd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt2dd_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
   }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gvtdd_m,2*GLB_T);
  /************************************ END of set gvt_dd_-*********************/


  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_dd_-[%d] = %.10e,%.10e\n",ix0,gvtdd_m[ix0].re,gvtdd_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}



void rotated_gvtudm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex *gvtud_m=corr_mem->g1, *gvt1ud_m=corr_mem->g2, *gvt2ud_m=corr_mem->g3, **M=corr_mem->M;

  complex ***gvt1_ij=corr_mem->g1_ij,  ***gvt2_ij=corr_mem->g2_ij;

  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;





/****************    Computation of gvt_ud_-   

gvtud- = - gvt1ud- + gvt2ud-


    gvt1ud- = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_uu        g5 P- U(x)    H(x+a0;a,y)_dd U^dag(0,y)psi      psi^dag Q+ csi)

    gvt2ud- = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_uu     g5 P+ U^dag(x)     H(x;a,y)_dd U^dag(0,y)psi    psi^dag Q+ csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gvt1ud_m[ix0]);
      _complex_0(gvt2ud_m[ix0]);
      _complex_0(gvtud_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(gvt1_ij[ix0][s1][s2]);
	  _complex_0(gvt2_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q+*/
      stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[0]);

    }
  }
  /**** construction of gij with open indeces ij ******/

  /****  gvt1ud_- *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i2);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);

	     /****  U0(x) prop_ud  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
                _suNf_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]);
             } 

	    stmp1[1]=stmp1[0];
 	    /****  P-  *****/
	    _vector_add_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_add_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_add_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_add_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2ud_- *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i2);

	     /****  U0(x) prop_ud  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
               _suNf_inverse_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]); 
             }

	    stmp1[1]=stmp1[0];   
	     /****  P+  *****/
	    _vector_sub_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_sub_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_sub_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_sub_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T-1;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt1ud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt2ud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

/*******   Sum over k ************/
   for(ix0=0;ix0<T-1;ix0++){
	_complex_sub_assign(gvtud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt1ud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(gvtud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt2ud_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
   }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gvtud_m,2*GLB_T);
  /************************************ END of set gvt_ud_-*********************/




 for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_ud_-[%d] = %.10e,%.10e\n",ix0,gvtud_m[ix0].re,gvtud_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}



void rotated_gvtdum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex *gvtdu_m=corr_mem->g1, *gvt1du_m=corr_mem->g2, *gvt2du_m=corr_mem->g3, **M=corr_mem->M;

  complex ***gvt1_ij=corr_mem->g1_ij,  ***gvt2_ij=corr_mem->g2_ij;

  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;



/****************    Computation of gvt_du_-  

gvtdu- = - gvt1du- + gvt2du-


    gvt1du- = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_dd       g5 P- U(x)    H(x+a0;a,y)_uu U^dag(0,y)psi      psi^dag Q- csi)

    gvt2du- = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_dd     g5 P+ U^dag(x)     H(x;a,y)_uu U^dag(0,y)psi    psi^dag Q- csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(gvt1du_m[ix0]);
      _complex_0(gvt2du_m[ix0]);
      _complex_0(gvtdu_m[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(gvt1_ij[ix0][s1][s2]);
	  _complex_0(gvt2_ij[ix0][s1][s2]);
	}
    }

  /*  Gamma matrix structure:    \csi^dag_j  Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
        
      /*Q-*/
      stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[0]);

    }
  }
  /**** construction of gij with open indeces ij ******/

  /****  gvt1du_- *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i2);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);

	     /****  U0(x) prop_du  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
                _suNf_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]);
             } 

	    stmp1[1]=stmp1[0];
 	    /****  P-  *****/
	    _vector_add_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_add_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_add_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_add_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2du_- *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T-1;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);

	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i2);

	     /****  U0(x) prop_du  ****/
             uptr=pu_gauge_f(i,0);    
             for(j=0;j<4;j++) {
               _suNf_inverse_multiply(stmp1[0].c[j],*uptr,sptr1[0]->c[j]); 
             }

	    stmp1[1]=stmp1[0];   
	     /****  P+  *****/
	    _vector_sub_assign_f(stmp1[1].c[0],stmp1[0].c[2]);
	    _vector_sub_assign_f(stmp1[1].c[1],stmp1[0].c[3]);
	    _vector_sub_assign_f(stmp1[1].c[2],stmp1[0].c[0]);
	    _vector_sub_assign_f(stmp1[1].c[3],stmp1[0].c[1]);
	    _spinor_mul_f(stmp1[1],.5,stmp1[1]);

	      /*gamma_5*/
	      stmp1[2].c[0]=stmp1[1].c[0];
	      stmp1[2].c[1]=stmp1[1].c[1];
	      stmp1[2].c[2]=stmp1[1].c[2];
	      stmp1[2].c[3]=stmp1[1].c[3];
	      _vector_minus_g(stmp1[2].c[2],stmp1[2].c[2]);
	      _vector_minus_g(stmp1[2].c[3],stmp1[2].c[3]);


	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[2]);
	      _complex_add_assign(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T-1;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,gvt1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt1du_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,gvt2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(gvt2du_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

/*******   Sum over k ************/
   for(ix0=0;ix0<T-1;ix0++){
	_complex_sub_assign(gvtdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt1du_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
	_complex_add_assign(gvtdu_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],gvt2du_m[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
   }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)gvtdu_m,2*GLB_T);
  /************************************ END of set gvt_du_-*********************/

 for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_du_-[%d] = %.10e,%.10e\n",ix0,gvtdu_m[ix0].re,gvtdu_m[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}



void rotated_g1uum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex g1uu_m;

  complex **g1_ij=corr_mem->g1_ij[0], **M=corr_mem->M;

  int i,j,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  


  /*************************  g1uu_m **********************************/

  /*   g1uu- = -1/2 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;T-2,z')     U0(T-2,z') g0 Q- U0(T-2,y')^dag    
       Huu^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi     */


  /****Initialization*****/
  _complex_0(g1uu_m);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[s1][s2]);
      _complex_0(M[s1][s2]);
    }
  
  
  /*  Gamma matrix structure:    \csi^dag_j  g0 Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
      
      /*Q+*/
      stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      /*gamma_0*/
      stmp1[1].c[0]=stmp1[0].c[2];
      stmp1[1].c[1]=stmp1[0].c[3];
      stmp1[1].c[2]=stmp1[0].c[0];
      stmp1[1].c[3]=stmp1[0].c[1];
      _spinor_mul_f(stmp1[1],-1.0,stmp1[1]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }


  /*  U0^dagg prop_uu[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr1[0] = _FIELD_AT(&prop_uu[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi1[4*NF].c[j],*uptr,sptr1[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi1[s],chi1[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);}
  }
  /*  U0^dagg prop2[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr2[0] = _FIELD_AT(&prop_dd[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi2[4*NF].c[j],*uptr,sptr2[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi2[s],chi2[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);}
  }


  global_sum((double*)chi1,4*NF*4*NF*2);    
  global_sum((double*)chi2,4*NF*4*NF*2);
        
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      stmp1[0]=chi1[s1];
      /*Q-*/
      _vector_i_sub_assign_f(stmp1[0].c[0],chi1[s1].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi1[s1].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi1[s1].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi1[s1].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);
      /*g0*/
      chi1[4*NF].c[0]=stmp1[0].c[2];
      chi1[4*NF].c[1]=stmp1[0].c[3];
      chi1[4*NF].c[2]=stmp1[0].c[0];
      chi1[4*NF].c[3]=stmp1[0].c[1];
      _spinor_mul_f(chi1[4*NF],-1.0,chi1[4*NF]);


      _spinor_prod_f(g1_ij[s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[s2][s1],M[s1][s2]);
      _complex_add_assign(g1uu_m,temp_comp);

    }
  }
  _complex_mulr(g1uu_m,(-0.5/GLB_VOL3)/GLB_VOL3,g1uu_m);


  /*****   Average over processors in temporal direction *****/

  /************************************ END of set gvt_du_-*********************/

  lprintf("PC_twisted_AC",10," g1_uu_- = %.10e,%.10e\n",g1uu_m.re,g1uu_m.im);	
  lprintf("PC_twisted_AC",10,"\n");
  
}



void rotated_g1ddm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex g1dd_m;

  complex  **g1_ij=corr_mem->g1_ij[0], **M=corr_mem->M;

  int i,j,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  


  /*************************  g1dd_m **********************************/

  /*   g1dd- = -1/2 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;T-2,z')     U0(T-2,z') g0 Q+ U0(T-2,y')^dag    
       Hdd^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi     */


  /****Initialization*****/
  _complex_0(g1dd_m);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[s1][s2]);
      _complex_0(M[s1][s2]);
    }
  
  
  /*  Gamma matrix structure:    \csi^dag_j  g0 Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
      
      /*Q-*/
      stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      /*gamma_0*/
      stmp1[1].c[0]=stmp1[0].c[2];
      stmp1[1].c[1]=stmp1[0].c[3];
      stmp1[1].c[2]=stmp1[0].c[0];
      stmp1[1].c[3]=stmp1[0].c[1];
      _spinor_mul_f(stmp1[1],-1.0,stmp1[1]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[1]);

    }
  }


  /*  U0^dagg prop_dd[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr1[0] = _FIELD_AT(&prop_dd[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi1[4*NF].c[j],*uptr,sptr1[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi1[s],chi1[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);}
  }
  /*  U0^dagg prop2[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr2[0] = _FIELD_AT(&prop_uu[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi2[4*NF].c[j],*uptr,sptr2[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi2[s],chi2[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);}
  }


  global_sum((double*)chi1,4*NF*4*NF*2);    
  global_sum((double*)chi2,4*NF*4*NF*2);
        
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      stmp1[0]=chi1[s1];
      /*Q+*/
      _vector_i_add_assign_f(stmp1[0].c[0],chi1[s1].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi1[s1].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi1[s1].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi1[s1].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);
      /*g0*/
      chi1[4*NF].c[0]=stmp1[0].c[2];
      chi1[4*NF].c[1]=stmp1[0].c[3];
      chi1[4*NF].c[2]=stmp1[0].c[0];
      chi1[4*NF].c[3]=stmp1[0].c[1];
      _spinor_mul_f(chi1[4*NF],-1.0,chi1[4*NF]);


      _spinor_prod_f(g1_ij[s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[s2][s1],M[s1][s2]);
      _complex_add_assign(g1dd_m,temp_comp);

    }
  }
  _complex_mulr(g1dd_m,(-0.5/GLB_VOL3)/GLB_VOL3,g1dd_m);

  lprintf("PC_twisted_AC",10," g1_dd_- = %.10e,%.10e\n",g1dd_m.re,g1dd_m.im);	
  lprintf("PC_twisted_AC",10,"\n");
}


void rotated_g1udm(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex g1ud_m;

  complex  **g1_ij=corr_mem->g1_ij[0],  **M=corr_mem->M;

  int i,j,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  

  /*************************  g1ud_m **********************************/

  /*   g1ud- = -1/2 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;T-2,z')     U0(T-2,z') Q- U0(T-2,y')^dag    
       Huu^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dagQ- csi     */


  /****Initialization*****/
  _complex_0(g1ud_m);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[s1][s2]);
      _complex_0(M[s1][s2]);
    }
  
  
  /*  Gamma matrix structure:    \csi^dag_j Q- csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
      
      /*Q-*/
      stmp1[0]=chi[s2];
      _vector_i_sub_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[0]);

    }
  }


  /*  U0^dagg prop_uu[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr1[0] = _FIELD_AT(&prop_uu[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi1[4*NF].c[j],*uptr,sptr1[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi1[s],chi1[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);}
  }
  /*  U0^dagg prop2[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr2[0] = _FIELD_AT(&prop_uu[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi2[4*NF].c[j],*uptr,sptr2[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi2[s],chi2[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);}
  }


  global_sum((double*)chi1,4*NF*4*NF*2);    
  global_sum((double*)chi2,4*NF*4*NF*2);
        
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      stmp1[0]=chi1[s1];
      /*Q-*/
      _vector_i_sub_assign_f(stmp1[0].c[0],chi1[s1].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi1[s1].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi1[s1].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi1[s1].c[1]);
      _spinor_mul_f(chi1[4*NF],.5,stmp1[0]);


      _spinor_prod_f(g1_ij[s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[s2][s1],M[s1][s2]);
      _complex_add_assign(g1ud_m,temp_comp);

    }
  }
  _complex_mulr(g1ud_m,(0.5/GLB_VOL3)/GLB_VOL3,g1ud_m);

  lprintf("PC_twisted_AC",10," g1_ud_- = %.10e,%.10e\n",g1ud_m.re,g1ud_m.im);	
  lprintf("PC_twisted_AC",10,"\n");
}


void rotated_g1dum(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {


  complex g1du_m;

  complex  **g1_ij=corr_mem->g1_ij[0], **M=corr_mem->M;

  int i,j,ix1,ix2,ix3,s1,s2;

  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;
  suNf *uptr;
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  

  /*************************  g1du_m **********************************/

  /*   g1du- = -1/2 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;T-2,z')     U0(T-2,z') Q+ U0(T-2,y')^dag    
       Hdd^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dagQ+ csi     */


  /****Initialization*****/
  _complex_0(g1du_m);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[s1][s2]);
      _complex_0(M[s1][s2]);
    }
  
  
  /*  Gamma matrix structure:    \csi^dag_j Q+ csi_i    */

  for(s1=0;s1<4*NF;s1++){
    for(s2=0;s2<4*NF;s2++){
      
      /*Q+*/
      stmp1[0]=chi[s2];
      _vector_i_add_assign_f(stmp1[0].c[0],chi[s2].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi[s2].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi[s2].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi[s2].c[1]);
      _spinor_mul_f(stmp1[0],.5,stmp1[0]);

      _spinor_prod_f(M[s1][s2],chi[s1],stmp1[0]);

    }
  }

  /*  U0^dagg prop_dd[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr1[0] = _FIELD_AT(&prop_dd[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi1[4*NF].c[j],*uptr,sptr1[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi1[s],chi1[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi1[s]);}
  }
  /*  U0^dagg prop2[s]*/
  if(COORD[0]==NP_T-1){
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);
      for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	    i=ipt((T-2),ix1,ix2,ix3);
	    sptr2[0] = _FIELD_AT(&prop_dd[s],i);
	    uptr=pu_gauge_f(i,0);
	    for(j=0;j<4;j++) {
	      _suNf_inverse_multiply(chi2[4*NF].c[j],*uptr,sptr2[0]->c[j]);
	    }
        
	    _spinor_add_assign_f(chi2[s],chi2[4*NF]);
        
	  }
    }
  } else {
    for(int s=0;s<4*NF;s++){
      _spinor_zero_f(chi2[s]);}
  }


  global_sum((double*)chi1,4*NF*4*NF*2);    
  global_sum((double*)chi2,4*NF*4*NF*2);
        
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      stmp1[0]=chi1[s1];
      /*Q+*/
      _vector_i_add_assign_f(stmp1[0].c[0],chi1[s1].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi1[s1].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi1[s1].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi1[s1].c[1]);
      _spinor_mul_f(chi1[4*NF],.5,stmp1[0]);


      _spinor_prod_f(g1_ij[s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[s2][s1],M[s1][s2]);
      _complex_add_assign(g1du_m,temp_comp);

    }
  }
  _complex_mulr(g1du_m,(0.5/GLB_VOL3)/GLB_VOL3,g1du_m);

  lprintf("PC_twisted_AC",10," g1_du_- = %.10e,%.10e\n",g1du_m.re,g1du_m.im);	
  lprintf("PC_twisted_AC",10,"\n");
}




#endif
