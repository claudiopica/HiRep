/******************************************************************************
 *
 * Methods to compute disconnected loops
 * Copyright (c) 2014, R. Arthur, V. Drach, A. Hietanen 
 * All rights reserved.
 * 
 *******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "disconnected.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#define PI 3.141592653589793238462643383279502884197


#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif


void measure_bilinear_loops_4spinorfield(spinor_field* prop,spinor_field* source,int k,int nm,int tau,int col,int eo)
{
				complex* corr[16];
				double* corr_re[16];
				double* corr_im[16];
				suNf_spin_matrix sma,smb, sm1;
				int size=nm*GLB_T;
				int i,ix,t,x,y,z,tc,beta;
				int tau_min=0;
				int tau_max=GLB_T;
			  complex tr;
				int NGamma=16;

				int offset=0 ; /* set to 0 just to copy paste some code. I don't know what it  is  originally used for */
				struct timeval start, end, etime;

				gettimeofday(&start,0);

				

				if (nm != 1) error(nm != 1, 1,"[measure_biliniear_loops_4spinorfield]", "Multimass not implemented !");

				for (i=0;i<NGamma;++i){
								corr[i]=(complex*) malloc(sizeof(complex)*size);
								corr_re[i]=(double*) malloc(sizeof(double)*size);
								corr_im[i]=(double*) malloc(sizeof(double)*size);
				}

				for (i=0;i<NGamma;++i) for (t=0;t<nm*GLB_T;t++)
				{ 
								corr[i][t].re=0.0;
								corr[i][t].im=0.0;
								corr_re[i][t]=0.0;
								corr_im[i][t]=0.0;
				}


				if (tau!=-1)
				{
								tau_min = tau;
								tau_max = tau+1;
				}

				/* loop on nm if the MMS is used */
				for(i=0; i<nm; i++) 
				{
								for (t=0; t<T; t++)
								{
												/* offset is set to zero here */
												tc = (zerocoord[0]+t+GLB_T)%GLB_T+i*GLB_T+offset;
												/* loop on spatial volume */
												for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
																ix=ipt(t,x,y,z);					

																for (beta=0;beta<4;beta++){
																				_spinmatrix_assign_row(sma, *_FIELD_AT(&source[beta*nm+i],ix), beta);
																				_spinmatrix_assign_row(smb, *_FIELD_AT(&prop[beta*nm+i],ix), beta);
																}

																_g5_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[0][tc], tr);

																_g1_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[1][tc], tr);
																_g2_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[2][tc], tr);
																_g3_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[3][tc], tr);

																_g5g0_spinmatrix(sm1, smb); // minus sign is missing
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[4][tc], tr);

																_g0g1_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[5][tc], tr);
																_g0g2_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[6][tc], tr);
																_g0g3_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[7][tc], tr);

																_spinmatrix_mul_trace(tr, sma, smb);
																_complex_add_assign(corr[8][tc], tr);

																_g5g1_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[9][tc], tr);
																_g5g2_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[10][tc], tr);
																_g5g3_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[11][tc], tr);


																_g0_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[12][tc], tr);
																_g5g0g1_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[13][tc], tr);
																_g5g0g2_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[14][tc], tr);
																_g5g0g3_spinmatrix(sm1, smb);
																_spinmatrix_mul_trace(tr, sma, sm1);
																_complex_add_assign(corr[15][tc], tr);
												}

								}/* end loop t */
				} /* end loop i (nm) */

				int iGamma;
				complex tmp;				
				for (t=0;t<T;t++){ 

								tc = (zerocoord[0]+t+GLB_T)%GLB_T+0*GLB_T+offset;
								for (iGamma=0;iGamma<NGamma;iGamma++) {
												if (iGamma!= 0  && iGamma!=1 &&iGamma!=2 && iGamma!=3 && iGamma!=8 && iGamma!=12)	
												{
																/* multiply by -i to match old convention [and to have hermitean operators ] */
																tmp.re = corr[iGamma][tc].im;
																tmp.im = -corr[iGamma][tc].re;
																corr[iGamma][tc].re=tmp.re;
																corr[iGamma][tc].im=tmp.im;

												}

												corr_re[iGamma][tc]=corr[iGamma][tc].re;
												corr_im[iGamma][tc]=corr[iGamma][tc].im;


								}
				}


				/* print into output file */
				for(iGamma=0; iGamma<NGamma; iGamma++) {
								global_sum(corr_re[iGamma],GLB_T*nm);
								global_sum(corr_im[iGamma],GLB_T*nm);
				}

				if (k==0  && tau==-1 &&  col==-1 && eo ==-1) lprintf("CORR",0," Output format [t iGamma iSource Re Im] \n");   // no dilution 
				if (k==0 && tau==0   && col==-1 && eo ==-1) lprintf("CORR",0,"Output format  [t iGamma iSource Re Im] \n"); // time dilution 
				if (k==0 && tau==0   && col==0 && eo ==-1) lprintf("CORR",0,"Output format  [t iGamma iSource col Re Im] \n"); // time + col dilution 
				if (k==0 && tau==0   && col==0 && eo ==0) lprintf("CORR",0,"Output format  [t iGamma iSource col eo Re Im] \n"); // time + col + eo dilution 
				if (k==0 && tau==-1   && col==0 && eo ==-1) lprintf("CORR",0,"Output format  [t iGamma iSource col Re Im] \n"); //   col dilution 
				if (k==0 && tau==-1  && col==0 && eo ==0) lprintf("CORR",0,"Output format  [t iGamma iSource col eo Re Im] \n"); // col + eo dilution 

				if (col ==-1&& eo ==-1)  for(t=tau_min;t<tau_max;++t) for(iGamma=0;iGamma<NGamma;iGamma++ ) lprintf("CORR",0,"%i %i %i %3.10e %3.10e \n",t,iGamma,k,corr_re[iGamma][t],corr_im[iGamma][t]);
				if (col !=-1&& eo ==-1)  for(t=tau_min;t<tau_max;++t) for(iGamma=0;iGamma<NGamma;iGamma++ ) lprintf("CORR",0,"%i %i %i %i %3.10e %3.10e \n",t,iGamma,k,col,corr_re[iGamma][t],corr_im[iGamma][t]);
				if (col ==-1&& eo !=-1)  for(t=tau_min;t<tau_max;++t) for(iGamma=0;iGamma<NGamma;iGamma++ ) lprintf("CORR",0,"%i %i %i %i %3.10e %3.10e \n",t,iGamma,k,eo,corr_re[iGamma][t],corr_im[iGamma][t]);
				if (col !=-1&& eo !=-1)  for(t=tau_min;t<tau_max;++t) for(iGamma=0;iGamma<NGamma;iGamma++ ) lprintf("CORR",0,"%i %i %i %i %i %3.10e %3.10e \n",t,iGamma,k,col,eo,corr_re[iGamma][t],corr_im[iGamma][t]);


				fflush(stdout); 
				gettimeofday(&end,0);
				timeval_subtract(&etime,&end,&start);
				if( tau != -1 ) lprintf("TIMING",0,"Contractions for source %i done and time %i [%ld sec %ld usec]\n",k,tau,etime.tv_sec,etime.tv_usec);
				if( tau == -1 ) lprintf("TIMING",0,"Contractions for source %i done and all timeslices [%ld sec %ld usec]\n",k,etime.tv_sec,etime.tv_usec);

}




void measure_loops(int nm, double* m, int nhits,int conf_num, double precision,int source_type,int n_mom)
{

				int k,l;
				int n_spinor;
				int eo,tau,col ;
				struct timeval start, end, etime;

				if (source_type ==0 )lprintf("CORR",0,"Pure volume source  will be used  \n");
				if (source_type ==1 )lprintf("CORR",0,"Gauge fixed source  with time and spin dilution will be used \n");
				if (source_type ==2 )lprintf("CORR",0,"Time and spin dilution  will be used \n");
				if (source_type ==3 )lprintf("CORR",0,"Time, spin and color dilution  will be used \n");
				if (source_type ==4 )lprintf("CORR",0,"Time, spin , color and eo dilution  will be used \n");
				if (source_type ==5 )lprintf("CORR",0,"Spin , color and eo dilution  will be used \n");

				gettimeofday(&start,0);
				init_propagator_eo(nm, m, precision);
				if (nm != 1 ) 
				{
								lprintf("ERR",0,"nm != 1  not tested !\n");
								exit(1);
				}



				for (k=0;k<nhits;k++){
								if (source_type == 0)  	/* generation of a volume source with Z2xZ2 noise */
								{
												spinor_field* source = alloc_spinor_field_f(1,&glattice);
												spinor_field* prop =  alloc_spinor_field_f(nm,&glattice);
												z2_spinor_field(source);
												start_sf_sendrecv(source);
												complete_sf_sendrecv(source);

												calc_propagator(prop,source,1);// No dilution 
												start_sf_sendrecv(prop);
												complete_sf_sendrecv(prop);

												lprintf("CORR",0,"Start to perform the contractions ... \n");
												measure_bilinear_loops_spinorfield(prop,source,k,nm,n_mom);
												lprintf("CORR",0,"Contraction done\n");
												free_spinor_field_f(source);
												free_spinor_field_f(prop);

								}

								if (source_type==1)   // experimental Gauge Fixed Wall sources
								{
												spinor_field* source = alloc_spinor_field_f(4,&glattice);
												spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
												suNg_field* u_gauge_old=alloc_gfield(&glattice);

												tau = 0;

												suNg_field_copy(u_gauge_old,u_gauge);
												//Fix the Gauge
												double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
																				1.8,    //overrelax
																				10000,  //maxit
																				1e-12, //tolerance
																				u_gauge //gauge
																				);
												lprintf("GFWALL",0,"Gauge fixed action  %1.6f\n",act);
												double p2 = calc_plaq(u_gauge);
												lprintf("TEST",0,"fixed_gauge plaq %1.6f\n",p2);
												full_plaquette();
												represent_gauge_field();
												init_propagator_eo(nm, m, precision);

												for (tau=0;tau<GLB_T;++tau){
																for (l=0;l<NF;++l){
																				create_gauge_fixed_wall_source(source, tau, l);
																				calc_propagator(prop,source,4);//4 for spin dilution
																				create_point_source(source,tau,l); //to get the contraction right
																				//		measure_discon(prop,source,nm,tau);

																}}
												//This gets the norm of the 2pt wrong by a factor GLB_VOL3 but the norm of the disconnected right
												//			print_mesons(GLB_T,conf_num,nm,m,"DISCON_GFWALL"); 

												suNg_field_copy(u_gauge,u_gauge_old);
												represent_gauge_field();

												//	free_propagator_eo(); 
												free_spinor_field_f(source);
												free_spinor_field_f(prop);
												//	free_gfield(u_gauge_old);

												//	measure_spectrum_discon_gfwall(nm,m,conf_num,precision); 
								}  /* gfwall */

								if (source_type==2)   
								{

												spinor_field* source = alloc_spinor_field_f(4,&glattice);
												spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);

												tau = 0;
												init_propagator_eo(nm, m, precision);

												for (tau=0;tau<GLB_T;++tau){	
																create_diluted_source_equal_atau(source,tau);	  

																calc_propagator(prop,source,4);//4 for spin dilution
																measure_bilinear_loops_4spinorfield(prop,source,k,nm,tau,-1,-1);

												} 
												free_spinor_field_f(source);
												free_spinor_field_f(prop);






								}  /* time + spin dilution  */
								if (source_type==3)   
								{

												spinor_field* source = alloc_spinor_field_f(4,&glattice);
												spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);

												tau = 0;
												init_propagator_eo(nm, m, precision);

												for (tau=0;tau<GLB_T;++tau){	 
																for (col=0;col<NF;++col){
																				create_diluted_source_equal_atau_col(source,tau,col);	  

																				calc_propagator(prop,source,4);//4 for spin dilution
																				measure_bilinear_loops_4spinorfield(prop,source,k,nm,tau,col,-1);

																}	
												} 

												free_spinor_field_f(source);
												free_spinor_field_f(prop);





								}  /* time + spin + color dilution  */
								if (source_type==4)   
								{

												spinor_field* source = alloc_spinor_field_f(4,&glattice);
												spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);


												tau = 0;
												n_spinor = 4;
												init_propagator_eo(nm, m, precision);

												for (tau=0;tau<GLB_T;++tau){	 
																for (col=0;col<NF;++col){
																				for (eo=0;eo<2;++eo)
																				{
																								create_diluted_source_equal_atau_col(source,tau,col);
																								zero_even_or_odd_site_spinorfield(source,n_spinor,eo);
																								calc_propagator(prop,source,4);//4 for spin dilution
																								measure_bilinear_loops_4spinorfield(prop,source,k,nm,tau,col,eo);

																				}

																}


												} 
												free_spinor_field_f(source);
												free_spinor_field_f(prop);



								}  /* time + spin + color +eo  dilution  */

								if (source_type==5)   
								{
												spinor_field* source = alloc_spinor_field_f(4,&glattice);
												spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);


												n_spinor = 4;
												init_propagator_eo(nm, m, precision);

												for (col=0;col<NF;++col){
																for (eo=0;eo<2;++eo)
																{

																				create_noise_source_equal_col_dil(source,col);
																				zero_even_or_odd_site_spinorfield(source,n_spinor,eo); //set even or odd site to zero
																				calc_propagator(prop,source,4);//4 for spin dilution
																				measure_bilinear_loops_4spinorfield(prop,source,k,nm,-1,col,eo);
																}
												}
												free_spinor_field_f(source);
												free_spinor_field_f(prop);


								}  /* volume source + spin + color + eo  dilution  */


				}

				free_propagator_eo(); 


				gettimeofday(&end,0);
				timeval_subtract(&etime,&end,&start);
				lprintf("TIMING",0,"Sources generation, invert and contract for %i sources done [%ld sec %ld usec]\n",nhits,etime.tv_sec,etime.tv_usec);


}


// set to zero even or odd site of the source.
void zero_even_or_odd_site_spinorfield(spinor_field *source,int nspinor,int eo)
{
				int c[4];
				int i;

				for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
								if(((zerocoord[0]+c[0]+zerocoord[1]+c[1]+zerocoord[2]+c[2]+zerocoord[3]+c[3])&1)==eo){
												for (i=0;i<nspinor;++i)
												{
																_vector_zero_g((_FIELD_AT(&source[0],ipt(c[0],c[1],c[2],c[3])))->c[i]);   
												}
								}
				} 
}






void measure_bilinear_loops_spinorfield(spinor_field* prop,spinor_field* source,int k,int nm, int n_mom)
{
  int px,py,pz,ip;
  int n_mom_tot = n_mom*n_mom*n_mom;
  double pdotx;
  complex phase;
  complex* corr[n_mom_tot][16];
  double* corr_re[n_mom_tot][16];
  double* corr_im[n_mom_tot][16];

  int pt[4];
  pt[0]=pt[1]=pt[2]=pt[3]=0;       

  int size=nm*GLB_T;
  int i,j,ix,t,x,y,z,tc;
  
  int NGamma=16;

  int offset=0 ; 
  suNf_spinor tmp_spinor;
  complex tmp;
  struct timeval start, end, etime;

  
  gettimeofday(&start,0);
  
  if (nm != 1) error(nm != 1, 1,"[measure_biliniear_loops_spinorfield]", "Multimass not implemented !");

 for (j=0;j<n_mom_tot;++j) for (i=0;i<NGamma;++i){
    corr[j][i]=(complex*) malloc(sizeof(complex)*size);
    corr_re[j][i]=(double*) malloc(sizeof(double)*size);
    corr_im[j][i]=(double*) malloc(sizeof(double)*size);
  }
  
  for (j=0;j<n_mom_tot;++j) for (i=0;i<NGamma;++i) for (t=0;t<nm*GLB_T;t++)
			   { 
			     corr[j][i][t].re=0.0;
			     corr[j][i][t].im=0.0;
			     corr_re[j][i][t]=0.0;
			     corr_im[j][i][t]=0.0;
			   }
  
				/* loop on nm if the MMS is used */
				for(i=0; i<nm; i++) 
				{
				  for (t=0; t<T; t++)
				  {        
				  /* offset set to zero here */
					  tc = (zerocoord[0]+t+GLB_T)%GLB_T+i*GLB_T+offset;
					  /* loop on spatial volume */
					  for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 

						  ix=ipt(t,x,y,z);                                        
						  ip = 0;
						  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz) { 
							  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1]-pt[1])/GLB_X + ((double) py)*(y+zerocoord[2]-pt[2])/GLB_Y + ((double) pz)*(z+zerocoord[3]-pt[3])/GLB_Z);
							  phase.re = cos(pdotx);
							  phase.im = sin(pdotx);

							  _spinor_g5_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][0][tc],phase,tmp);

							  _spinor_g1_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][1][tc],phase,tmp);

							  _spinor_g2_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][2][tc],phase,tmp);

							  _spinor_g3_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][3][tc],phase,tmp);

							  _spinor_g0g5_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][4][tc],phase,tmp);

							  _spinor_g0g1_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][5][tc],phase,tmp);

							  _spinor_g0g2_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][6][tc],phase,tmp);

							  _spinor_g0g3_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][7][tc],phase,tmp);

							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),*_FIELD_AT(prop,ix)); 
							  _complex_mul_assign(corr[ip][8][tc],phase,tmp);

							  _spinor_g5g1_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][9][tc],phase,tmp);

							  _spinor_g5g2_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][10][tc],phase,tmp);

							  _spinor_g5g3_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][11][tc],phase,tmp);

							  _spinor_g0_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor);
							  _complex_mul_assign(corr[ip][12][tc],phase,tmp);

							  _spinor_g5g0g1_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][13][tc],phase,tmp);

							  _spinor_g5g0g2_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][14][tc],phase,tmp);

							  _spinor_g5g0g3_f(tmp_spinor,*_FIELD_AT(prop,ix));
							  _spinor_prod_f(tmp,*_FIELD_AT(source,ix),tmp_spinor); 
							  _complex_mul_assign(corr[ip][15][tc],phase,tmp);


							  ip=ip+1;
						  }


					  }/* end loop t */
				  } /* end loop px,py,pz */
				} /* end loop i (nm) */


				// the following have to be modified as well !
				int iGamma;
				for (t=0;t<T;t++){ 

					tc = (zerocoord[0]+t+GLB_T)%GLB_T+0*GLB_T+offset;
					for (j=0;j<n_mom_tot;++j) 	for (iGamma=0;iGamma<NGamma;iGamma++) {

						if (iGamma!= 0  && iGamma!=1 &&iGamma!=2 && iGamma!=3 && iGamma!=8 && iGamma!=12)       
						{
							/* multiply by -i to match old convention [and to have hermitean operators ] */
							tmp.re = corr[j][iGamma][tc].im;
							tmp.im = -corr[j][iGamma][tc].re;
							corr[j][iGamma][tc].re=tmp.re;
							corr[j][iGamma][tc].im=tmp.im;

						}

						corr_re[j][iGamma][tc]=corr[j][iGamma][tc].re;
						corr_im[j][iGamma][tc]=corr[j][iGamma][tc].im;


					}

				}


				/* print into output file */
				for (j=0;j<n_mom_tot;++j)  for(iGamma=0; iGamma<NGamma; iGamma++) {
					global_sum(corr_re[j][iGamma],GLB_T*nm);
					global_sum(corr_im[j][iGamma],GLB_T*nm);
				}

				if (k==0 ) lprintf("CORR",0,"loops for one noise vector, all local bilinear [t iGamma iSource Re Im] \n"); 
					
				for(t=0;t<GLB_T;++t) for(iGamma=0;iGamma<NGamma;iGamma++ ){
					ip=0;				
					for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz) 
					{
						lprintf("CORR",0,"%i %i %i %i %i %i %3.10e %3.10e \n",t,iGamma,k,px,py,pz,corr_re[ip][iGamma][t],corr_im[ip][iGamma][t]);
						ip = ip+1;
					}
				}
				fflush(stdout); 
				gettimeofday(&end,0);
				timeval_subtract(&etime,&end,&start);
				lprintf("TIMING",0,"Contractions for source %i done [%ld sec %ld usec]\n",k,etime.tv_sec,etime.tv_usec);

}
