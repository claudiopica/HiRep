#ifdef ROTATED_SF
#include "global.h"
#include "communications.h"
#include "observables.h"
#include "logger.h"
#include <stdlib.h>

chisf_mem * init_rotated_corr_mem(){


   int i,j;
   chisf_mem * corrmem=malloc(sizeof(chisf_mem));
   error(corrmem==NULL,1,"init_rotated_corr_mem","Could not allocate memory space for mem");   
   
   corrmem->g1=malloc(16*GLB_T*sizeof(complex));
   error(corrmem->g1==NULL,1,"init_rotated_corr_mem","Could not allocate memory space for mem");   
   corrmem->g4 =corrmem->g1+GLB_T;
   corrmem->g2 =corrmem->g4+GLB_T;
   corrmem->g3 =corrmem->g2+GLB_T;
   corrmem->l41=corrmem->g3+GLB_T;
   corrmem->l21=corrmem->l41+GLB_T;
   corrmem->l31=corrmem->l21+GLB_T;
   corrmem->l11=corrmem->l31+GLB_T;
   corrmem->l42=corrmem->l11+GLB_T;
   corrmem->l22=corrmem->l42+GLB_T;
   corrmem->l32=corrmem->l22+GLB_T;
   corrmem->l12=corrmem->l32+GLB_T;
   corrmem->l43=corrmem->l12+GLB_T;
   corrmem->l23=corrmem->l43+GLB_T;
   corrmem->l33=corrmem->l23+GLB_T;
   corrmem->l13=corrmem->l33+GLB_T;

   
   corrmem->g1_ij =malloc(13*GLB_T*sizeof(complex**));
   error(corrmem->g1_ij==NULL,1,"init_rotated_corr_mem","Could not allocate memory space for mem");   
   corrmem->g2_ij =corrmem->g1_ij +GLB_T;
   corrmem->g3_ij =corrmem->g2_ij +GLB_T;
   corrmem->g4_ij =corrmem->g3_ij +GLB_T;
   corrmem->l21_ij=corrmem->g4_ij +GLB_T;
   corrmem->l31_ij=corrmem->l21_ij+GLB_T;
   corrmem->l11_ij=corrmem->l31_ij+GLB_T;
   corrmem->l22_ij=corrmem->l11_ij+GLB_T;
   corrmem->l32_ij=corrmem->l22_ij+GLB_T;
   corrmem->l12_ij=corrmem->l32_ij+GLB_T;
   corrmem->l23_ij=corrmem->l12_ij+GLB_T;
   corrmem->l33_ij=corrmem->l23_ij+GLB_T;
   corrmem->l13_ij=corrmem->l33_ij+GLB_T;

   
   complex** tmp=malloc(GLB_T*13*4*NF*sizeof(complex*));
   complex* tmp1=malloc(4*NF*GLB_T*13*4*NF*sizeof(complex));
   error(tmp==NULL,1,"init_rotated_corr_mem","Could not allocate memory space for mem");
   error(tmp1==NULL,1,"init_rotated_corr_mem","Could not allocate memory space for mem");   
   for(i=0;i<GLB_T;i++){
     corrmem->g1_ij[i] =tmp+i*4*NF;
     corrmem->g4_ij[i] =tmp+i*4*NF+4*NF*GLB_T;
     corrmem->g2_ij[i] =tmp+i*4*NF+2*4*NF*GLB_T;
     corrmem->g3_ij[i] =tmp+i*4*NF+3*4*NF*GLB_T;
     corrmem->l21_ij[i]=tmp+i*4*NF+4*4*NF*GLB_T;
     corrmem->l31_ij[i]=tmp+i*4*NF+5*4*NF*GLB_T;
     corrmem->l11_ij[i]=tmp+i*4*NF+6*4*NF*GLB_T;
     corrmem->l22_ij[i]=tmp+i*4*NF+7*4*NF*GLB_T;
     corrmem->l32_ij[i]=tmp+i*4*NF+8*4*NF*GLB_T;
     corrmem->l12_ij[i]=tmp+i*4*NF+9*4*NF*GLB_T;
     corrmem->l23_ij[i]=tmp+i*4*NF+10*4*NF*GLB_T;
     corrmem->l33_ij[i]=tmp+i*4*NF+11*4*NF*GLB_T;
     corrmem->l13_ij[i]=tmp+i*4*NF+12*4*NF*GLB_T;

     for(j=0;j<4*NF;j++){
       corrmem->g1_ij[i][j] =tmp1+j*4*NF+i*4*NF*4*NF;		 
       corrmem->g2_ij[i][j] =tmp1+j*4*NF+i*4*NF*4*NF+4*NF*4*NF*GLB_T;	 
       corrmem->g4_ij[i][j] =tmp1+j*4*NF+i*4*NF*4*NF+4*NF*2*4*NF*GLB_T; 
       corrmem->g3_ij[i][j] =tmp1+j*4*NF+i*4*NF*4*NF+4*NF*3*4*NF*GLB_T; 
       corrmem->l21_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*4*4*NF*GLB_T; 
       corrmem->l31_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*5*4*NF*GLB_T; 
       corrmem->l11_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*6*4*NF*GLB_T; 
       corrmem->l22_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*7*4*NF*GLB_T; 
       corrmem->l32_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*8*4*NF*GLB_T; 
       corrmem->l12_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*9*4*NF*GLB_T; 
       corrmem->l23_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*10*4*NF*GLB_T;
       corrmem->l33_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*11*4*NF*GLB_T;
       corrmem->l13_ij[i][j]=tmp1+j*4*NF+i*4*NF*4*NF+4*NF*12*4*NF*GLB_T;
       
     }
   }
   (corrmem->l13_ij)[GLB_T-1][4*NF-1][4*NF-1].re=1.0;
   corrmem->M=malloc(4*NF*sizeof(complex*));
   for(j=0;j<4*NF;j++) corrmem->M[j]=malloc(4*NF*sizeof(complex));

   return corrmem; 
}

 void rotated_gXuup(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  
  /**************************************************
    set:   gauu_p, gpuu_p, gsuu_p, gvuu_p

gpuu+ = -  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       II      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       

gauu+ =  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       

gsuu+ = -  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       

gvuu+ = -  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi       
  **************************************************/
   complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij, ***g3_ij=corr_mem->g3_ij, ***g4_ij=corr_mem->g4_ij, *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, *g4=corr_mem->g4, **M=corr_mem->M;
  int i,ix0,ix1,ix2,ix3,s1,s2;
  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;

  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g4[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g1[ix0]);
      _complex_0(g3[ix0]);
    
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g4_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
	  _complex_0(g3_ij[ix0][s1][s2]);
	  _complex_0(g1_ij[ix0][s1][s2]);
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
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

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
	      _complex_add_assign(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],.5/GLB_VOL3,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g4[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g4,2*GLB_T);
  global_sum((double*)g2,2*GLB_T);
  global_sum((double*)g1,2*GLB_T);
  global_sum((double*)g3,2*GLB_T);

  /************************************ END of set gX_uu_+*********************/


  /* gA_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_uu_+[%d] = %.10e,%.10e\n",ix0,g4[ix0].re,g4[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  /* gV_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_uu_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
 
  /* gP_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_uu_+[%d] = %.10e,%.10e\n",ix0,g2[ix0].re,g2[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

  /* gS_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_uu_+[%d] = %.10e,%.10e\n",ix0,g1[ix0].re,g1[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}






 void rotated_gXddp(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij, ***g3_ij=corr_mem->g3_ij, ***g4_ij=corr_mem->g4_ij, *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, *g4=corr_mem->g4, **M=corr_mem->M;
  int i,ix0,ix1,ix2,ix3,s1,s2;
  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;

  /**************************************************
    set:   gadd_p, gpdd_p, gsdd_p, gvdd_p

gpdd+ = -  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       II      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       

gadd+ =  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       

gsdd+ = -  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       

gvdd+ = -  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g4[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g1[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g4_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
	  _complex_0(g3_ij[ix0][s1][s2]);
	  _complex_0(g1_ij[ix0][s1][s2]);
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
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

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
	      _complex_add_assign(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],.5/GLB_VOL3,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g4[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g4,2*GLB_T);
  global_sum((double*)g2,2*GLB_T);
  global_sum((double*)g1,2*GLB_T);
  global_sum((double*)g3,2*GLB_T);

  /************************************ END of set gX_dd_+*********************/


  /* gA_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_dd_+[%d] = %.10e,%.10e\n",ix0,g4[ix0].re,g4[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  /* gV_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_dd_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
 
  /* gP_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_dd_+[%d] = %.10e,%.10e\n",ix0,g2[ix0].re,g2[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

  /* gS_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_dd_+[%d] = %.10e,%.10e\n",ix0,g1[ix0].re,g1[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}

 void rotated_gXudp(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij, ***g3_ij=corr_mem->g3_ij, ***g4_ij=corr_mem->g4_ij, *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, *g4=corr_mem->g4, **M=corr_mem->M;
  int i,ix0,ix1,ix2,ix3,s1,s2;
  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;

  /**************************************************
    set:   gaud_p, gpud_p, gsud_p, gvud_p

gpud+ =  1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       II      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       

gaud+ = - 1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       

gsud+ =   1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       

gvud+ =   1/2 sum  csi^dag U0(1,z) Huu^(-1)(2,z;x)       g5g0      Hdd^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q- csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g4[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g1[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g4_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
	  _complex_0(g3_ij[ix0][s1][s2]);
	  _complex_0(g1_ij[ix0][s1][s2]);
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
	      sptr1[0] = _FIELD_AT(&prop_dd[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_dd[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

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
	      _complex_add_assign(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g4[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g4,2*GLB_T);
  global_sum((double*)g2,2*GLB_T);
  global_sum((double*)g1,2*GLB_T);
  global_sum((double*)g3,2*GLB_T);

  /************************************ END of set gX_ud_+*********************/


  /* gA_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_ud_+[%d] = %.10e,%.10e\n",ix0,g4[ix0].re,g4[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  /* gV_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_ud_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
 
  /* gP_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_ud_+[%d] = %.10e,%.10e\n",ix0,g2[ix0].re,g2[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

  /* gS_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_ud_+[%d] = %.10e,%.10e\n",ix0,g1[ix0].re,g1[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}

 void rotated_gXdup(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij, ***g3_ij=corr_mem->g3_ij, ***g4_ij=corr_mem->g4_ij, *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, *g4=corr_mem->g4, **M=corr_mem->M;
  int i,ix0,ix1,ix2,ix3,s1,s2;
  
  suNf_spinor stmp1[2];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  complex temp_comp;


  /**************************************************
    set:   gadu_p, gpdu_p, gsdu_p, gvdu_p

gpud+ =  1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       II      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       

gaud+ = - 1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       

gsud+ =   1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       

gvud+ =   1/2 sum  csi^dag U0(1,z) Hdd^(-1)(2,z;x)       g5g0      Huu^(-1)(x:2,y) U0(1,y)^dag psi      psi^dag Q+ csi       
  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g4[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g1[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g4_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
	  _complex_0(g3_ij[ix0][s1][s2]);
	  _complex_0(g1_ij[ix0][s1][s2]);
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
	      sptr1[0] = _FIELD_AT(&prop_uu[s1],i);
	      sptr2[0] = _FIELD_AT(&prop_uu[s2],i);
 

	      /*gpuu*/
	      _spinor_prod_f(temp_comp,*sptr2[0],*sptr1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  

	      /*gauu*/
	      /*gamma_0*/
	      stmp1[0].c[0]=sptr1[0]->c[2];
	      stmp1[0].c[1]=sptr1[0]->c[3];
	      stmp1[0].c[2]=sptr1[0]->c[0];
	      stmp1[0].c[3]=sptr1[0]->c[1];
	      _spinor_mul_f(stmp1[0],-1.0,stmp1[0]);

	      _spinor_prod_f(temp_comp,*sptr2[0],stmp1[0]);
	      //  acumulate for the spatial average
	      _complex_add_assign(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);     

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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);

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
	      _complex_add_assign(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
       
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
	_complex_mulr(g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g4_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g4[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g3_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }



  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g4,2*GLB_T);
  global_sum((double*)g2,2*GLB_T);
  global_sum((double*)g1,2*GLB_T);
  global_sum((double*)g3,2*GLB_T);

  /************************************ END of set gX_du_+*********************/


  /* gA_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," ga_du_+[%d] = %.10e,%.10e\n",ix0,g4[ix0].re,g4[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
  
  /* gV_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gv_du_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");
 
  /* gP_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gp_du_+[%d] = %.10e,%.10e\n",ix0,g2[ix0].re,g2[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

  /* gS_f1f2_+ */
  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gs_du_+[%d] = %.10e,%.10e\n",ix0,g1[ix0].re,g1[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}

 void rotated_gvtuup(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij, *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, **M=corr_mem->M;
  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;


  /****************    Computation of gvt_uu_+

gvtuu+ = - gvt1uu+ + gvt2uu+


    gvt1uu+ = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_uu        g5 P- U(x)    H(x+a0;a,y)_uu U^dag(0,y)psi      psi^dag g0 Q- csi)

    gvt2uu+ = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_uu     g5 P+ U^dag(x)     H(x;a,y)_uu U^dag(0,y)psi    psi^dag g0 Q- csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g1[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g1_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
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
 
  /****  gvt1uu_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
	for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	      i=ipt(ix0,ix1,ix2,ix3);
              i2=iup(i,0);
	      //ipt(ix0+1,ix1,ix2,ix3);

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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }


  /****  gvt2uu_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      }
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

  /*******   Sum over k ************/
  for(ix0=0;ix0<T;ix0++){
    _complex_sub_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
    _complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g3,2*GLB_T);
  /************************************ END of set gvt_uu_+*********************/

  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_uu_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");



}

 void rotated_gvtddp(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij, *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, **M=corr_mem->M;
  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;

  /****************    Computation of gvt_dd_+   

gvtdd+ = - gvt1dd+ + gvt2dd+


    gvt1dd+ = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_dd        g5 P- U(x)    H(x+a0;a,y)_dd U^dag(0,y)psi      psi^dag g0 Q+ csi)

    gvt2dd+ = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_dd     g5 P+ U^dag(x)     H(x;a,y)_dd U^dag(0,y)psi    psi^dag g0 Q+ csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g1[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g1_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
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

  /****  gvt1dd_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2dd_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],-0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

  /*******   Sum over k ************/
  for(ix0=0;ix0<T;ix0++){
    _complex_sub_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
    _complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g3,2*GLB_T);
  /************************************ END of set gvt_dd_+*********************/

  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_dd_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");



}






 void rotated_gvtudp(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij,  *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3, **M=corr_mem->M;
  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;




  /****************    Computation of gvt_ud_+   

gvtud+ = - gvt1ud+ + gvt2ud+


    gvt1ud+ = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_uu        g5 P- U(x)    H(x+a0;a,y)_dd U^dag(0,y)psi      psi^dag Q- csi)

    gvt2ud+ = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_uu     g5 P+ U^dag(x)     H(x;a,y)_dd U^dag(0,y)psi    psi^dag Q- csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g1[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g1_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
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

  /****  gvt1ud_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2ud_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

  /*******   Sum over k ************/
  for(ix0=0;ix0<T;ix0++){
    _complex_sub_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
    _complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g3,2*GLB_T);
  /************************************ END of set gvt_ud_+*********************/



  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_ud_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");



}


 void rotated_gvtdup(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, ***g2_ij=corr_mem->g2_ij,  *g1=corr_mem->g1, *g2=corr_mem->g2, *g3=corr_mem->g3,  **M=corr_mem->M;
  int i,j,i2,ix0,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;



  /****************    Computation of gvt_du_+   

gvtdu+ = - gvt1du+ + gvt2du+


    gvt1du+ = -1/2 Sum Tr (  csi^dag U(0,z) H(a,z;x)_dd       g5 P- U(x)    H(x+a0;a,y)_uu U^dag(0,y)psi      psi^dag Q+ csi)

    gvt2du+ = -1/2 Sum Tr (csi^dag U(0,z) H(a,z;x+a0)_dd     g5 P+ U^dag(x)     H(x;a,y)_uu U^dag(0,y)psi    psi^dag Q+ csi)
              

  **************************************************/


  /****Initialization*****/
  for(ix0=0;ix0<GLB_T;ix0++)
    {
      _complex_0(g1[ix0]);
      _complex_0(g2[ix0]);
      _complex_0(g3[ix0]);
     
      for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
	  _complex_0(g1_ij[ix0][s1][s2]);
	  _complex_0(g2_ij[ix0][s1][s2]);
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

  /****  gvt1du_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /****  gvt2du_+ *****/
  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){

      for(ix0=0;ix0<T;ix0++){
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
	      _complex_add_assign(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],temp_comp);  
	    }
        // Normalization with 0.5/GLB_VOL3
	_complex_mulr(g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],0.5/GLB_VOL3,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1]);
      } 
    }
  }


  /***** Here we do the contractions of gX_ij M_ij    ****/


  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      for(ix0=0;ix0<T;ix0++){
	// Contract and accumulate
	_complex_mul(temp_comp,g1_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);

	_complex_mul(temp_comp,g2_ij[(zerocoord[0]+ix0-1+GLB_T)%GLB_T][s2][s1],M[s1][s2]);
	_complex_add_assign(g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],temp_comp);
      }
    }
  }

  /*******   Sum over k ************/
  for(ix0=0;ix0<T;ix0++){
    _complex_sub_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g1[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
    _complex_add_assign(g3[(zerocoord[0]+ix0-1+GLB_T)%GLB_T],g2[(zerocoord[0]+ix0-1+GLB_T)%GLB_T]);
  }


  /*****   Average over processors in temporal direction *****/
  global_sum((double*)g3,2*GLB_T);
  /************************************ END of set gvt_du_+*********************/



  for(ix0=0;ix0<GLB_T-1;ix0++)
    lprintf("PC_twisted_AC",10," gvt_du_+[%d] = %.10e,%.10e\n",ix0,g3[ix0].re,g3[ix0].im);	
  lprintf("PC_twisted_AC",10,"\n");



}




 void rotated_g1uup(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij,  *g1=corr_mem->g1,  **M=corr_mem->M;
  int i,j,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;

  /*************************  g1uu_P **********************************/

  /*   g1uu+ = -1/2 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;T-2,z')     U0(T-2,z') g0 Q+ U0(T-2,y')^dag    
       Huu^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dag g0 Q- csi     */


  /****Initialization*****/
  _complex_0(g1[0]);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[0][s1][s2]);
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


      _spinor_prod_f(g1_ij[0][s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[0][s2][s1],M[s1][s2]);
      _complex_add_assign(g1[0],temp_comp);

    }
  }
  _complex_mulr(g1[0],(-0.5/GLB_VOL3)/GLB_VOL3,g1[0]);

  lprintf("PC_twisted_AC",10," g1_uu_+ = %.10e,%.10e\n",g1[0].re,g1[0].im);	
  lprintf("PC_twisted_AC",10,"\n");


}





 void rotated_g1ddp(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, *g1=corr_mem->g1, **M=corr_mem->M;
  int i,j,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;

  /*************************  g1dd_P **********************************/

  /*   g1dd+ = -1/2 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;T-2,z')     U0(T-2,z') g0 Q- U0(T-2,y')^dag
       Hdd^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dag g0 Q+ csi     */


  /****Initialization*****/
  _complex_0(g1[0]);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[0][s1][s2]);
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


      _spinor_prod_f(g1_ij[0][s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[0][s2][s1],M[s1][s2]);
      _complex_add_assign(g1[0],temp_comp);

    }
  }
  _complex_mulr(g1[0],(-0.5/GLB_VOL3)/GLB_VOL3,g1[0]);

  lprintf("PC_twisted_AC",10," g1_dd_+ = %.10e,%.10e\n",g1[0].re,g1[0].im);	
  lprintf("PC_twisted_AC",10,"\n");

}



 void rotated_g1udp(chisf_mem* corr_mem,suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij,  *g1=corr_mem->g1, **M=corr_mem->M;
  int i,j,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;

  /*************************  g1ud_P **********************************/

  /*   g1ud+ = -1/2 sum tr csi^dag U0(1,z) Hdd^(-1)(2,z;T-2,z')     U0(T-2,z') Q+ U0(T-2,y')^dag
       Huu^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dagQ+ csi     */


  /****Initialization*****/
  _complex_0(g1[0]);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[0][s1][s2]);
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
      /*Q+*/
      _vector_i_add_assign_f(stmp1[0].c[0],chi1[s1].c[2]);
      _vector_i_add_assign_f(stmp1[0].c[1],chi1[s1].c[3]);
      _vector_i_sub_assign_f(stmp1[0].c[2],chi1[s1].c[0]);
      _vector_i_sub_assign_f(stmp1[0].c[3],chi1[s1].c[1]);
      _spinor_mul_f(chi1[4*NF],.5,stmp1[0]);


      _spinor_prod_f(g1_ij[0][s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[0][s2][s1],M[s1][s2]);
      _complex_add_assign(g1[0],temp_comp);

    }
  }
  _complex_mulr(g1[0],(0.5/GLB_VOL3)/GLB_VOL3,g1[0]);

  lprintf("PC_twisted_AC",10," g1_ud_+ = %.10e,%.10e\n",g1[0].re,g1[0].im);	
  lprintf("PC_twisted_AC",10,"\n");



}

 void rotated_g1dup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd) {
  complex ***g1_ij=corr_mem->g1_ij, *g1=corr_mem->g1,**M=corr_mem->M;
  int i,j,ix1,ix2,ix3,s1,s2;
  suNf_spinor stmp1[3];
  suNf_spinor chi1[4*NF+1];  
  suNf_spinor chi2[4*NF+1];  
  suNf_spinor *sptr2[2];
  suNf_spinor *sptr1[2];
  suNf *uptr;
  complex temp_comp;

  /*************************  g1du_P **********************************/

  /*   g1du+ = -1/2 sum tr csi^dag U0(1,z) Huu^(-1)(2,z;T-2,z')     U0(T-2,z') Q- U0(T-2,y')^dag
       Hdd^(-1)(T-2,y',2,y) U0(1,y)^dag psi      psi^dagQ- csi     */


  /****Initialization*****/
  _complex_0(g1[0]);
  for(s1=0;s1<4*NF;s1++)for(s2=0;s2<4*NF;s2++){
      _complex_0(g1_ij[0][s1][s2]);
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
      /*Q-*/
      _vector_i_sub_assign_f(stmp1[0].c[0],chi1[s1].c[2]);
      _vector_i_sub_assign_f(stmp1[0].c[1],chi1[s1].c[3]);
      _vector_i_add_assign_f(stmp1[0].c[2],chi1[s1].c[0]);
      _vector_i_add_assign_f(stmp1[0].c[3],chi1[s1].c[1]);
      _spinor_mul_f(chi1[4*NF],.5,stmp1[0]);


      _spinor_prod_f(g1_ij[0][s2][s1],chi2[s2],chi1[4*NF]);
      //  acumulate for the spatial average
    }
  }

  for(s2=0;s2<4*NF;s2++){
    for(s1=0;s1<4*NF;s1++){
      // Contract and accumulate
      _complex_mul(temp_comp,g1_ij[0][s2][s1],M[s1][s2]);
      _complex_add_assign(g1[0],temp_comp);

    }
  }
  _complex_mulr(g1[0],(0.5/GLB_VOL3)/GLB_VOL3,g1[0]);

  lprintf("PC_twisted_AC",10," g1_du_+ = %.10e,%.10e\n",g1[0].re,g1[0].im);	
  lprintf("PC_twisted_AC",10,"\n");



}


#endif
