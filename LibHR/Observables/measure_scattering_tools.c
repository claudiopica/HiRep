#include <string.h>
#include "utils.h"
#include "suN.h"
#include "memory.h"
#include "io.h"
#include "logger.h"
#include "global.h"
#include "meson_observables.h"
#include "observables.h"
#include "scattering.h"
#include "random.h"
#include "spinor_field.h"
#include <communications.h>
#include "linear_algebra.h"
#include "hr_complex.h"
/**
 *
 * @file scatter_functions.h
 * Contains structures and functions used in the scattering code
 *
 * @author Tadeusz Janowski
 */



/**
 * @brief Converts a string of "(px,py,pz)(px2,py2,pz2)..." into a 2D array of integers.
 * @param momstring input string
 * @param N number of momenta in the string (used as an output)
 */
int** getmomlist(char* momstring, int* N){
    char* tmp = momstring;
    *N = 0;
    while(tmp != NULL){
        tmp = strchr(tmp+1,'(');
        (*N)++;
    }
    int** plist = (int**) malloc(*N*sizeof(int*));
    int i=0;
    tmp = momstring;
    lprintf("getmomlist",0,"%d %s\n",*N,tmp);
    while(tmp != NULL){
        plist[i] = (int*) malloc(3*sizeof(int));
        sscanf(tmp,"(%d,%d,%d)", plist[i], plist[i]+1, plist[i]+2);
        lprintf("getmomlist",0,"(%d,%d,%d)\n",*(plist[i]),*(plist[i]+1),*(plist[i]+2));
        tmp = strchr(tmp+1,'(');
        i++;
    }
    return plist;
}

/** 
 * @brief Frees the 2D array of momenta allocated by getmomlist.
 * @param p array of momenta
 * @param N number of momenta
 * @see getmomlist
 */
void freep(int **p, int N){
    for(int i=0; i<N;i++){
        free(p[i]);
    }
    free(p);
}
//#endif


/** 
 * @brief Function for initiating meson observable (used to store the correlation function)
 * @param mo meson_observable to initialise.
 * @param name name of the channel (e.g. "pi").
 * @param size Size of the output, typically GLB_T*(number of output momenta)
 */
void init_mo(meson_observable* mo, char* name, int size)
{
  int i;
  //ind1 and ind2 don't do anything for the moment
  mo->ind1 = _g5;
  mo->ind2 = _g5;
  strcpy(mo->channel_name,name);
  strcpy(mo->channel_type, "Pi Pi scattering");
  mo->sign=1.0;
  mo->corr_size=size;
  mo->corr_re = (double * ) malloc(size * sizeof(double));
  mo->corr_im = (double * ) malloc(size * sizeof(double));
  
  if (mo->corr_re == NULL || mo->corr_im == NULL)
  {
    fprintf(stderr, "malloc failed in init_mo \n");
    return;
  }
  mo->next=NULL;
  for (i=0; i<size; ++i)
  {
    mo->corr_re[i]=0.0;
    mo->corr_im[i]=0.0;
  }
}

/**
 * @brief Resets the meson_observable object (sets all entries to 0).
 * @param mo meson_observable to reset
 */
void reset_mo(meson_observable* mo)
{
  int i;
  for (i=0; i< mo->corr_size; ++i)
  {
    mo->corr_re[i]=0.0;
    mo->corr_im[i]=0.0;
  }
}

/**
 * @brief Sums the entries of meson_observable over all MPI processes
 * @param mo meson_observable to sum over
 * @param norm number to multiply by after the sum
 */
static void do_global_sum(meson_observable* mo, double norm){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
      global_sum(motmp->corr_re,motmp->corr_size);
      global_sum(motmp->corr_im,motmp->corr_size);

      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] *= norm;
	motmp->corr_im[i] *= norm;

      }
    motmp=motmp->next;
 }
}

/**
 * @brief Flips the boundary conditions in the T direction between periodic and anti-periodic.
 * @param tau Time slice corresponding to the boundary.
 */
static void flip_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  tau-=1;
  if (tau<0) tau+= GLB_T;
  lprintf("meson_measurements",15,"Flipping the boundary at global time slice %d\n",tau);
  fflush(stdout);
  if((zerocoord[0]-1<=tau && zerocoord[0]+T>tau) || (zerocoord[0]==0 && tau==GLB_T-1)) { 
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
	  if ((tau==zerocoord[0]-1) || (zerocoord[0]==0 && tau==GLB_T-1)){
	    index=ipt_ext(0,ix,iy,iz);
	  }
	  else{
	    index=ipt_ext(T_BORDER+tau-zerocoord[0],ix,iy,iz);
	  }
	  if(index!=-1) {
	    u=pu_gauge_f(index,0);
	    _suNf_minus(*u,*u);
	  }
	}
  }
  lprintf("meson_measurements",50,"Flipping DONE!\n");
}

/**
 * @brief Frees the memory allocated by the meson observable.
 */
void free_mo(meson_observable* mo)
{
  free(mo->corr_re);
  free(mo->corr_im);
  free(mo);
}

/// \cond
#define BASENAME(filename) (strrchr((filename),'/') ? strrchr((filename),'/')+1 : filename )

#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*GLB_T*(nm)+(py)*(n_mom)*GLB_T*(nm)+(pz)*GLB_T*(nm)+ ((cm)*GLB_T) +(tc))
/// \endcond

/// \cond
#define INDEX(px,py,pz,n_mom,tc) ((px + n_mom)*(2*n_mom+1)*(2*n_mom+1)*(GLB_T)+(py + n_mom)*(2*n_mom+1)*(GLB_T)+(pz + n_mom)*(GLB_T)+ (tc))
/// \endcond

// Sources

/**
 * @brief Allocates memory and creates point sources with momentum 0 (used for testing).
 * @param src output
 * @param tau time slice at which the sources are generated
 */
void init_src_common_point(struct src_common* src, int tau){
	// Malloc and initialise
	src->src_0 = alloc_spinor_field_f(4, &glattice);
	src->src_0_eta = alloc_spinor_field_f(4, &glattice);
	src->src_0_0 = alloc_spinor_field_f(4, &glattice);
	for(int i=0; i< 4; i++){
		spinor_field_zero_f(&(src->src_0[i]));
		spinor_field_zero_f(&(src->src_0_eta[i]));
		spinor_field_zero_f(&(src->src_0_0[i]));
	}
	
	// Create point source
  	create_point_source(src->src_0, 0, 0);
  	create_point_source(src->src_0_eta, 0, 0);
}

/**
 * @brief Allocates memory and creates noise sources with momentum 0.
 * @param src output
 * @param tau time slice at which the sources are generated
 */
void init_src_common(struct src_common* src, int tau){
	// Malloc and initialise
	src->src_0 = alloc_spinor_field_f(4, &glattice);
	src->src_0_eta = alloc_spinor_field_f(4, &glattice);
	src->src_0_0 = alloc_spinor_field_f(4, &glattice);
	for(int i=0; i< 4;i++){
		spinor_field_zero_f(&(src->src_0[i]));
		spinor_field_zero_f(&(src->src_0_eta[i]));
		spinor_field_zero_f(&(src->src_0_0[i]));
	}

	// Create point source
  	create_diluted_source_equal_atau(src->src_0, tau);
  	create_diluted_source_equal_atau(src->src_0_eta, tau);
}

/**
 * @brief Allocates memory and creates noise sources with momentum p.
 * @param srcp output
 * @param src0 source at zero momentum. This needs to initialised by init_src_common prior to executing this function.
 * @param tau time slice at which the sources are generated
 * @param px,py,pz momentum of the generated propagators
 * @see init_src_common
 */
void init_src_p(struct src_p* srcp, struct src_common* src0, int px, int py, int pz){
	srcp->p[0]=px;
	srcp->p[1]=py;
	srcp->p[2]=pz;

	srcp->src_p = alloc_spinor_field_f(4, &glattice);
	srcp->src_mp = alloc_spinor_field_f(4, &glattice);
	srcp->src_0_p = alloc_spinor_field_f(4, &glattice);
	srcp->src_p_0 = alloc_spinor_field_f(4, &glattice);
	srcp->src_0_mp = alloc_spinor_field_f(4, &glattice);
	srcp->src_mp_0 = alloc_spinor_field_f(4, &glattice);

	for(int i=0; i< 4; i++){
		spinor_field_zero_f(&(srcp->src_p[i]));
		spinor_field_zero_f(&(srcp->src_mp[i]));
		spinor_field_zero_f(&(srcp->src_p_0[i]));
		spinor_field_zero_f(&(srcp->src_0_p[i]));
		spinor_field_zero_f(&(srcp->src_mp_0[i]));
		spinor_field_zero_f(&(srcp->src_0_mp[i]));
	}

	add_momentum(srcp->src_p, src0->src_0, px, py, pz);
	add_momentum(srcp->src_mp, src0->src_0, -px, -py, -pz);
    lprintf("init_src_p",0,"Source p generated\n");
}

/**
 * @brief Frees memory allocated by to src_common
 */
void free_src_common(struct src_common* src){
	free_spinor_field_f(src->src_0);
	free_spinor_field_f(src->src_0_eta);
	free_spinor_field_f(src->src_0_0);
    lprintf("free_src_common",0,"Freed memory\n");
}
/**
 * @brief Frees memory allocated by to src_p
 */
void free_src_p(struct src_p* src){
	free_spinor_field_f(src->src_p);
	free_spinor_field_f(src->src_mp);
	free_spinor_field_f(src->src_p_0);
	free_spinor_field_f(src->src_mp_0);
	free_spinor_field_f(src->src_0_p);
	free_spinor_field_f(src->src_0_mp);
    lprintf("free_src_p",0,"Freed memory\n");
}

// Propagators

/**
 * @brief Creates a propagator with periodic boundary conditions
 * @param prop output
 * @param src source of the propagator
 * @param ndilute number of dilution vectors, e.g. 4 for spin dilution
 * @param tau origin time slice
 */
void make_propagator_P(spinor_field* prop, spinor_field* src, int ndilute, int tau){
	calc_propagator(prop,src,ndilute);
}

/**
 * @brief Creates a propagator with P+A boundary conditions
 * @param prop output
 * @param src source of the propagator
 * @param ndilute number of dilution vectors, e.g. 4 for spin dilution
 * @param tau origin time slice
 */
void make_propagator_PA(spinor_field* prop, spinor_field* src, int ndilute, int tau){
	spinor_field* proptmp = alloc_spinor_field_f(4,&glattice);

	calc_propagator(prop,src,ndilute);
	flip_T_bc(tau);
	calc_propagator(proptmp,src,ndilute);
	flip_T_bc(tau);
	for(int l=0;l<4;++l){
		spinor_field_add_assign_f(&prop[l],&proptmp[l]);
		spinor_field_mul_f(&prop[l],0.5,&prop[l]);
	}

    free_spinor_field_f(proptmp);
}

/**
 * @brief Creates a propagator bundle with zero momentum
 * @param prop output
 * @param src0 source of the propagators
 * @param ndilute number of dilution vectors, e.g. 4 for spin dilution
 * @param tau origin time slice
 * @param bc if set to "PA" generates P+A bcs in time direction, defaults to periodic otherwise
 */
void make_prop_common(struct prop_common* prop, struct src_common* src0, int ndilute, int tau, char* bc){
    void (*fun) (spinor_field*, spinor_field*,int,int);
    //if(bc == "PA"){
    if(strcmp(bc,"PA") == 0){
        lprintf("make_prop_common",0,"Inverting propagator with P+A boundary conditions");
        fun = &make_propagator_PA;
    } else {
        lprintf("make_prop_common",0,"Inverting propagator with P boundary conditions");
        fun = &make_propagator_P;
    }
	prop->Q_0 = alloc_spinor_field_f(4,&glattice);
    fun(prop->Q_0, src0->src_0, ndilute, tau);
	prop->Q_0_eta = alloc_spinor_field_f(4,&glattice);
    fun(prop->Q_0_eta, src0->src_0_eta, ndilute, tau);
    prop->W_0_0 = (spinor_field**) malloc(GLB_T*sizeof(spinor_field*));
    for(int t=0; t<GLB_T; t++){
        prop->W_0_0[t] = alloc_spinor_field_f(4,&glattice);
        spinor_field* tmp = alloc_spinor_field_f(4,&glattice);
        create_sequential_source_stoch(tmp,(t+tau)%GLB_T,prop->Q_0);
        fun(prop->W_0_0[t], tmp, 4, tau);
        free_spinor_field_f(tmp);
    }
    create_sequential_source_stoch(src0->src_0_0,tau,prop->Q_0);
    lprintf("MAIN",0,"Propagator with momentum 0 inverted\n");
}

/**
 * @brief Creates a propagator bundle with momentum p
 * @param prop output
 * @param srcp source of the propagators with momentum p
 * @param src0 source of the propagators with momentum 0
 * @param ndilute number of dilution vectors, e.g. 4 for spin dilution
 * @param tau origin time slice
 * @param bc if set to "PA" generates P+A bcs in time direction, defaults to periodic otherwise
 */
void make_prop_p(struct prop_p* prop, struct src_p* srcp, struct src_common* src0, int ndilute, int tau, char* bc){
    void (*fun) (spinor_field*, spinor_field*,int,int);
    if(strcmp(bc,"PA") == 0){
        lprintf("make_prop_p",0,"Inverting propagator with P+A boundary conditions");
        fun = &make_propagator_PA;
    } else {
        lprintf("make_prop_p",0,"Inverting propagator with P boundary conditions");
        fun = &make_propagator_P;
    }
    /// \cond
#define PROPLIST X(Q_p) X(Q_mp) X(W_0_p) X(W_0_mp) X(W_p_0) X(W_mp_0)
#define X(NAME) \
	prop->NAME = alloc_spinor_field_f(4, &glattice);
	PROPLIST
#undef X
#undef PROPLIST
    fun(prop->Q_p, srcp->src_p, ndilute, tau);
    fun(prop->Q_mp, srcp->src_mp, ndilute, tau);

	create_sequential_source_stoch(srcp->src_p_0,tau,prop->Q_p);
	create_sequential_source_stoch(srcp->src_mp_0,tau,prop->Q_mp);
	add_momentum(srcp->src_0_p, src0->src_0_0, srcp->p[0], srcp->p[1], srcp->p[2]);
	add_momentum(srcp->src_0_mp, src0->src_0_0,-(srcp->p[0]), -(srcp->p[1]), -(srcp->p[2]));
#define PROP X(0_p) X(0_mp) X(p_0) X(mp_0)
#define X(NAME) \
    fun(prop -> W_##NAME, srcp->src_##NAME, ndilute, tau);
	PROP
#undef X
#undef PROP
        /// \endcond
    lprintf("MAIN",0,"Propagator with momentum p inverted\n");

}

/**
 * @brief Frees memory associated with zero-momentum propagator bundle
 */
void free_prop_common(struct prop_common* prop){
	free_spinor_field_f( prop -> Q_0);
	free_spinor_field_f( prop -> Q_0_eta);
    for(int t=0; t<GLB_T; t++){
        free_spinor_field_f( prop -> W_0_0[t]);
    }
    free(prop -> W_0_0);
    lprintf("free_prop_common",0,"Freed memory\n");
}

/**
 * @brief Frees memory associated with momentum-p propagator bundle
 */
void free_prop_p(struct prop_p* prop){
	free_spinor_field_f( prop -> Q_p);
	free_spinor_field_f( prop -> Q_mp);
	free_spinor_field_f( prop -> W_0_p);
	free_spinor_field_f( prop -> W_0_mp);
	free_spinor_field_f( prop -> W_p_0);
	free_spinor_field_f( prop -> W_mp_0);
    lprintf("free_prop_p",0,"Freed memory\n");
}

// Meson observable stuff

/**
 * @brief Initialises bundle of zero-momentum meson observables.
 */
void init_mo_0(struct mo_0* mo){
	mo->pi = (meson_observable*) malloc(sizeof(meson_observable));
    lprintf("init_mo_0",0,"pi initiated\n");
	init_mo(mo->pi,"pi",27*GLB_T);
	for(int i=0; i<3; i++) for (int j=0;j<3;j++){
		mo->rho[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
		init_mo(mo->rho[i][j],"rho",27*GLB_T);
      	mo->rho[i][j]->ind1 = i+3;
		mo->rho[i][j]->ind2 = j+3;
        lprintf("init_mo_0",0,"rho %d,%d initiated\n",i,j);
	}
    lprintf("MAIN",0,"Meson observable 0 initiated!\n");
}

/**
 * @brief Initialises a bundle of meson_observables with momentum p
 * @param px,py,pz momentum components
 */
void init_mo_p(struct mo_p* mo, int px, int py, int pz){
	mo->p[0] = px;
	mo->p[1] = py;
	mo->p[2] = pz;
/// \cond
#define R X(d) X(pi) X(r1) X(r2) X(r3) X(r4)
#define X(NAME)\
       	mo->NAME = (meson_observable*) malloc(sizeof(meson_observable));\
	init_mo(mo->NAME,#NAME,27*GLB_T);
	R
#undef X
#undef R
/// \endcond
	for(int i=0; i<3; i++){
		mo->t1[i] = (meson_observable*) malloc(sizeof(meson_observable));
		init_mo(mo->t1[i],"t1",27*GLB_T);
		mo->t1[i]->ind2 = i+3;

		mo->t2[i] = (meson_observable*) malloc(sizeof(meson_observable));
		init_mo(mo->t2[i],"t2",27*GLB_T);
		mo->t2[i]->ind2 = i+3;
        lprintf("init_mo_p",0,"Allocated t gamma=%d\n",i);
        
		for(int j=0; j<3; j++){
			mo->rho[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
			init_mo(mo->rho[i][j],"rho",27*GLB_T);
			mo->rho[i][j]->ind1 = i+3;
			mo->rho[i][j]->ind2 = j+3;
            lprintf("init_mo_p",0,"Allocated rho gamma=%d,%d\n",i,j);
		}
	}
    lprintf("MAIN",0,"Meson observable p initiated!\n");
}

/**
 * @brief Generates zero-momentum contractions
 * @param mo output
 * @param p0 zero-momentum propagator bundle
 * @param s0 zero-momentum source
 * @param tau origin time slice
 */
void gen_mo_0(struct mo_0* mo, struct prop_common* p0, struct src_common* s0, int tau){
	measure_mesons_core(p0->Q_0, p0->Q_0, s0->src_0, mo->pi, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo->pi,1.0);
	for(int i=0;i<3; i++) for(int j=0; j<3; j++){
		measure_mesons_core(p0->Q_0, p0->Q_0, s0->src_0, mo->rho[i][j], 1, tau, 2, 0, GLB_T);
		do_global_sum(mo->rho[i][j],1.0);
	}
    lprintf("MAIN",0,"Generated mo 0\n");
}

/**
 * @brief Generates contractions with source momentum p
 * @param mo output
 * @param p0 zero-momentum propagator bundle
 * @param pp p-momentum propagator bundle
 * @param s0 zero-momentum source
 * @param tau origin time slice
 */
void gen_mo_p(struct mo_p* mo, struct prop_common* p0, struct prop_p* pp, struct src_common* s0, int tau){
	measure_mesons_core(p0->Q_0, pp->Q_p, s0->src_0, mo->pi, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo->pi,1.0);
	measure_scattering_AD_core(mo->d, pp->Q_p, p0->Q_0, p0->Q_0_eta, p0->Q_0_eta, tau, 0, 1, mo->p[0], mo->p[1], mo->p[2] );
	do_global_sum(mo->d,1.0);
	for(int i=0;i<3;i++){
		measure_mesons_core(pp->W_0_mp, p0->Q_0, s0->src_0, mo->t1[i],1,tau, 2, 0, GLB_T);
		do_global_sum(mo->t1[i],1.0);

		measure_mesons_core(p0->Q_0, pp->W_0_p, s0->src_0, mo->t2[i],1,tau, 2, 0, GLB_T);
		do_global_sum(mo->t2[i],1.0);

		for(int j=0; j<3; j++){
			measure_mesons_core(p0->Q_0, pp->Q_p, s0->src_0, mo->rho[i][j],1,tau, 2, 0, GLB_T);
			do_global_sum(mo->rho[i][j],1.0);
		}
	}

    /// \cond
#define R X(r1) X(r2) X(r3) X(r4)
#define X(NAME) meson_observable* mo_tmp##NAME;
    R
#undef X
	for(int t=0; t<GLB_T; t++){
#define X(NAME)\
        mo_tmp##NAME = (meson_observable*) malloc(sizeof(meson_observable));\
		init_mo(mo_tmp##NAME,"tmp",27*GLB_T);
        R
#undef X
		measure_mesons_core(pp->W_mp_0, p0->W_0_0[t], s0->src_0, mo_tmpr1, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(p0->W_0_0[t], pp->W_p_0, s0->src_0, mo_tmpr2, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(p0->W_0_0[t], pp->W_0_p, s0->src_0, mo_tmpr3, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(pp->W_0_mp, p0->W_0_0[t], s0->src_0, mo_tmpr4, 1, tau, 2, 0, GLB_T);
		for(int px_=0;px_<2;++px_)for(int py_=0;py_<2;++py_)for(int pz_=0;pz_<2;++pz_){
#define X(NAME)\
				mo->NAME->corr_re[corr_ind(px_,py_,pz_,2,t,1,0)] = mo_tmp##NAME->corr_re[corr_ind(px_,py_,pz_,2,t,1,0)];\
				mo->NAME->corr_im[corr_ind(px_,py_,pz_,2,t,1,0)] = mo_tmp##NAME->corr_im[corr_ind(px_,py_,pz_,2,t,1,0)];
            R
#undef X
		}

#define X(NAME)\
		free_mo(mo_tmp##NAME);
        R
#undef X
	}
#define X(NAME)\
	do_global_sum(mo->NAME, 1.0); 
    R
#undef X
    lprintf("MAIN",0,"Generated mo p\n");
#undef R
    /// \endcond
}	

/**
 * @brief Resets the mo_p structure, setting all correlation functions to 0.
 */
void reset_mo_p(struct mo_p* mo){
/// \cond
#define R X(d) X(pi) X(r1) X(r2) X(r3) X(r4)
#define X(NAME)\
	reset_mo(mo->NAME);
	R
#undef X
#undef R
/// \endcond
	for(int i=0; i<3; i++){
		reset_mo(mo->t1[i]);
		reset_mo(mo->t2[i]);
		for(int j=0; j<3; j++){
			reset_mo(mo->rho[i][j]);
		}
	}
}

/**
 * @brief copy meson_observable structure.
 */
void copy_mo(meson_observable* mo_dest,meson_observable* mo_src ){
	
	for(int t=0; t<GLB_T; t++){
		mo_dest->corr_re[corr_ind(0,0,0,2,t,1,0)] = mo_src->corr_re[corr_ind(0,0,0,2,t,1,0)];
		mo_dest->corr_im[corr_ind(0,0,0,2,t,1,0)] = mo_src->corr_im[corr_ind(0,0,0,2,t,1,0)];
	}
}

/**
 * @brief Frees the memory allocated by mo_0 structure.
 */
void free_mo_0(struct mo_0* mo){
	free_mo(mo->pi);
	for(int i=0; i<3; i++) for (int j=0;j<3;j++){
		free_mo(mo->rho[i][j]);
	}
	free(mo);
}

/**
 * @brief Frees the memory allocated by mo_p structure.
 */
void free_mo_p(struct mo_p* mo){
/// \cond
#define R X(d) X(pi) X(r1) X(r2) X(r3) X(r4)
#define X(NAME)\
	free_mo(mo->NAME);
	R
#undef X
#undef R
/// \endcond
	for(int i=0; i<3; i++){
		free_mo(mo->t1[i]);
		free_mo(mo->t2[i]);
		for(int j=0; j<3; j++){
			free_mo(mo->rho[i][j]);
		}
	}
	free(mo);
}


/**
 * @brief Prints the 2-point function to a file. Not used for JSON output.
 * @param mo meson_observable to print
 * @param pmax maximum momentum at the sink
 * @param sourceno noise source number
 * @param path directory to which the output file will be saved
 * @param name name of the output file
 * @param cnfg_filename name of the configuration
 */
void io2pt(meson_observable* mo, int pmax, int sourceno, char* path, char* name,char * cnfg_filename)
{
	FILE* file;
	char outfile[256] = {};
	int px,py,pz,t;
	if(PID==0){
		sprintf(outfile,"%s/%s_src_%d_%s", path, name, sourceno, BASENAME(cnfg_filename) );
		file=fopen(outfile,"w+");
		//Factor of 2 to correct for the noise source normalisation
		for(px=0;px<pmax;++px) for(py=0;py<pmax;++py) for(pz=0;pz<pmax;++pz) for(t=0;t<GLB_T;++t) fprintf(file,"%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,2*(mo->corr_re[corr_ind(px,py,pz,pmax,t,1,0)]), 2*(mo->corr_im[corr_ind(px,py,pz,pmax,t,1,0)]));
		fclose(file);
	}
	return;
}

/**
 * @brief Prints the 4-point function (d) to a file. Not used for JSON output.
 * @param mo meson_observable to print
 * @param pmax maximum momentum at the sink
 * @param sourceno noise source number
 * @param path directory to which the output file will be saved
 * @param name name of the output file
 * @param cnfg_filename name of the configuration
 */
void io4pt(meson_observable* mo, int pmax, int sourceno, char* path, char* name,char * cnfg_filename)
{
	FILE* file;
	char outfile[256] = {};
	int px,py,pz,t;
	if(PID==0){
		sprintf(outfile,"%s/%s_src_%d_%s", path, name, sourceno, BASENAME(cnfg_filename) );
		file=fopen(outfile,"w+");
		//Factor of 4 to correct for the noise source normalisation
		for(px=-pmax;px<=pmax;++px) for(py=-pmax;py<=pmax;++py) for(pz=-pmax;pz<=pmax;++pz) for(t=0;t<GLB_T;++t) fprintf(file, "%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,4*(mo->corr_re[INDEX(px,py,pz, pmax,t)]), 4*(mo->corr_im[INDEX(px,py,pz,pmax,t)]));
		fclose(file);
	}
	return;
}

/**
 * @brief "Old style" IO where each correlation function is saved to separate file. Prints zero-momentum only.
 * @param molist an array of mo_0 objects, where each index corresponds to a different noise source
 * @param numsources number of noise sources
 * @param path path to which the file should be saved 
 * @param cnfg_filename name of the configuration
 * @see IO_json_0
 */
void IOold_0(struct mo_0* molist[], int numsources, char* path, char* cnfg_filename){
	for(int src=0; src<numsources; src++){
        lprintf("IOold_0",0,"Printing pi for source %d\n",src);
		io2pt(molist[src]->pi, 2, src, path, "pi",cnfg_filename);
		for(int i=0;i<3; i++){
			char tmp[100];
            lprintf("IOold_0",0,"Printing rho %d for source %d\n",i+1,src);
			sprintf(tmp, "rho_p0_g%d", i+1);
			io2pt(molist[src]->rho[i][i], 2, src, path, tmp,cnfg_filename);
		}
	}
}

/**
 * @brief "Old style" IO where each correlation function is saved to separate file. Prints momentum p contractions.
 * @param molist an array of mo_p objects, where each index corresponds to a different noise source
 * @param numsources number of noise sources
 * @param path path to which the file should be saved 
 * @param cnfg_filename name of the configuration
 *
 * @see IO_json_p
 */
void IOold_p(struct mo_p* molist[], int numsources, char* path, char* cnfg_filename	){
	char tmp[100];
	for(int src=0; src<numsources; src++){
        int px = molist[src]->p[0];
        int py = molist[src]->p[1];
        int pz = molist[src]->p[2];
        lprintf("IOold_p",0,"Printing pi for source %d momentum (%d,%d,%d)\n",src,px,py,pz);
		sprintf(tmp, "pi_p(%d,%d,%d)",px,py,pz);
		io2pt(molist[src]->pi, 2, src, path, tmp,cnfg_filename);
        lprintf("IOold_p",0,"Printing d for source %d momentum (%d,%d,%d)\n",src,px,py,pz);
		sprintf(tmp, "d_p(%d,%d,%d)",px,py,pz);
		io4pt(molist[src]->d, 1, src, path, tmp,cnfg_filename);
        lprintf("IOold_p",0,"Printing r1 for source %d momentum (%d,%d,%d)\n",src,px,py,pz);
		sprintf(tmp, "r1_p(%d,%d,%d)",px,py,pz);
		io2pt(molist[src]->r1,2,src,path, tmp,cnfg_filename);
        lprintf("IOold_p",0,"Printing r2 for source %d momentum (%d,%d,%d)\n",src,px,py,pz);
		sprintf(tmp, "r2_p(%d,%d,%d)",px,py,pz);
		io2pt(molist[src]->r2,2,src,path, tmp,cnfg_filename);
        lprintf("IOold_p",0,"Printing r3 for source %d momentum (%d,%d,%d)\n",src,px,py,pz);
		sprintf(tmp, "r3_p(%d,%d,%d)",px,py,pz);
		io2pt(molist[src]->r3,2,src,path, tmp,cnfg_filename);
        lprintf("IOold_p",0,"Printing r4 for source %d momentum (%d,%d,%d)\n",src,px,py,pz);
		sprintf(tmp, "r4_p(%d,%d,%d)",px,py,pz);
		io2pt(molist[src]->r4,2,src,path, tmp,cnfg_filename);

		for(int i=0;i<3; i++){
            lprintf("IOold_p",0,"Printing for source %d momentum (%d,%d,%d) gamma %d\n",src,px,py,pz,i+1);
			sprintf(tmp, "t1_p(%d,%d,%d)_g%d", px, py, pz, i+1);
			io2pt(molist[src]->t1[i], 2, src, path, tmp,cnfg_filename);
			sprintf(tmp, "t2_p(%d,%d,%d)_g%d", px, py, pz, i+1);
			io2pt(molist[src]->t2[i], 2, src, path, tmp,cnfg_filename);
		       	for(int j=0;j<3; j++){
				sprintf(tmp, "rho_p(%d,%d,%d)_g%d%d", px, py, pz, i+1, j+1);
                lprintf("IOold_p",0,"Printing for source %d momentum (%d,%d,%d) gamma %d\n",src,px,py,pz,i+1, j+1);
				io2pt(molist[src]->rho[i][j], 2, src, path, tmp,cnfg_filename);
			}
		}
	}
}


/// \cond
#define JSON(STRUCT, NAME) \
        fprintf(f,"\t\"%s\":{\n\t\t",NAME);\
		for(px=0;px<pmax;++px) for(py=0;py<pmax;++py) for(pz=0;pz<pmax;++pz){\
            fprintf(f,"\"(%d,%d,%d)\":[[\n",px,py,pz);\
            for(int src=0; src<numsources; src++){\
                for(int t=0;t<GLB_T;++t){\
                    fprintf(f,"\t\t\t\t[%e,%e]", 2*(molist[src]->STRUCT->corr_re[corr_ind(px,py,pz,pmax,t,1,0)]), 2*(molist[src]->STRUCT->corr_im[corr_ind(px,py,pz,pmax,t,1,0)]));\
                    if(t<GLB_T-1){ \
                        fprintf(f,",\n");\
                    } else {\
                        fprintf(f,"]");\
                    }\
                }\
                if(src<numsources-1){\
                    fprintf(f,",[\n");\
                } else {\
                    fprintf(f,"]");\
                }\
            }\
            if(px==pmax-1 && py==pmax-1 && pz==pmax-1){\
                fprintf(f,"\n");\
            } else {\
                fprintf(f,",\n\t\t");\
            }\
        }\
        fprintf(f,"\n\t}");

#define JSON_4PT(STRUCT, NAME) \
        fprintf(f,"\t\"%s\":{\n\t\t",NAME);\
		for(px=0;px<pmax;++px) for(py=0;py<pmax;++py) for(pz=0;pz<pmax;++pz){\
            fprintf(f,"\"(%d,%d,%d)\":[[\n",px,py,pz);\
            for(int src=0; src<numsources; src++){\
                for(int t=0;t<GLB_T;++t){\
                    fprintf(f,"\t\t\t\t[%e,%e]", 4*(molist[src]->STRUCT->corr_re[INDEX(px,py,pz,pmax,t)]), 4*(molist[src]->STRUCT->corr_im[INDEX(px,py,pz,pmax,t)]));\
                    if(t<GLB_T-1){ \
                        fprintf(f,",\n");\
                    } else {\
                        fprintf(f,"]");\
                    }\
                }\
                if(src<numsources-1){\
                    fprintf(f,",[\n");\
                } else {\
                    fprintf(f,"]");\
                }\
            }\
            if(px==pmax-1 && py==pmax-1 && pz==pmax-1){\
                fprintf(f,"\n");\
            } else {\
                fprintf(f,",\n\t\t");\
            }\
        }\
        fprintf(f,"\n\t}");
/// \endcond

/**
 * @brief Prints a bundle of meson_observables with momentum 0 to a json file.
 * @param molist an array of mo_0 objects, where each index corresponds to a different noise source
 * @param numsources number of noise sources
 * @param path path to which the file should be saved 
 * @param cnfg_filename name of the configuration 
 */
void IO_json_0(struct mo_0* molist[], int numsources, char* path,char * cnfg_filename){
    FILE* f;
	char outfile[256] = {};
	int px,py,pz;
    int pmax = 2;
	if(PID==0){
		sprintf(outfile,"%s/p_(0,0,0)_n%s.json", path, strrchr(cnfg_filename,'n') + 1 );
        lprintf("IO_json_0",0,"%s\n",outfile);
		f=fopen(outfile,"w+");
        fprintf(f,"{\n");
        JSON(pi, "pi"); fprintf(f,",\n");
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            char c[256];
            sprintf(c,"rho_g%d_g%d",i+1,j+1);
            JSON(rho[i][j], c); 
            if(i!=2 || j!=2) fprintf(f,",");
            fprintf(f,"\n");
        }

        fprintf(f,"}");
		fclose(f);
	}
	return;
}

/**
 * @brief Prints a bundle of meson_observables with momentum p to a json file.
 * @param molist an array of mo_p objects, where each index corresponds to a different noise source
 * @param numsources number of noise sources
 * @param path path to which the file should be saved 
 * @param cnfg_filename name of the configuration 
 */
void IO_json_p(struct mo_p* molist[], int numsources, char* path, char* cnfg_filename){
    FILE* f;
	char outfile[256] = {};
	int px,py,pz;
    int pmax = 2;
	if(PID==0){
		sprintf(outfile,"%s/p_(%d,%d,%d)_n%s.json", path, molist[0]->p[0], molist[0]->p[1], molist[0]->p[2],strrchr(cnfg_filename,'n') + 1 );
        lprintf("IO_json_0",0,"%s\n",outfile);
		f=fopen(outfile,"w+");
        fprintf(f,"{\n");
        JSON(pi, "pi"); fprintf(f,",\n");
        JSON_4PT(d, "d"); fprintf(f,",\n");
        JSON(r1, "r1"); fprintf(f,",\n");
        JSON(r2, "r2"); fprintf(f,",\n");
        JSON(r3, "r3"); fprintf(f,",\n");
        JSON(r4, "r4"); fprintf(f,",\n");
        for(int i=0;i<3;i++){
            char c[256];
            sprintf(c,"t1_g%d",i+1);
            JSON(t1[i], c); fprintf(f,",\n");
            sprintf(c,"t2_g%d",i+1);
            JSON(t2[i], c); fprintf(f,",\n");
            for(int j=0;j<3;j++){
                sprintf(c,"rho_g%d_g%d",i+1,j+1);
                JSON(rho[i][j], c); 
                if(i!=2 || j!=2) fprintf(f,",");
                fprintf(f,"\n");
            }
        }

        fprintf(f,"}");
		fclose(f);
	}
	return;
}
//#endif


/* Random timeslice not previously chosen */
static int random_tau(){
  static int* slices=NULL;
  if (slices == NULL) slices = (int*) malloc(GLB_T*sizeof(int));
  static int counter = 0;
  int itmp,tau,i;
  double ran;

  if (counter == 0){
    for (i=0;i<GLB_T;++i){
      slices[i]=i;
    }
    counter=GLB_T;
  }
  do{
    ranlxd(&ran,1);
    itmp=(int)(ran*counter);
  } while(itmp==counter);
  counter--;
  tau = slices[itmp];
  slices[itmp]=slices[counter];
  slices[counter]=tau;
  bcast_int(&tau,1);
  return tau;
}

static int random_tau2(){
double ran;
int tau;
ranlxd(&ran,1);
tau = (int)(ran*GLB_T);

bcast_int(&tau,1);

return tau;	
}

static int gi(int i){
	if(i==1) return _g1;
	if(i==2) return _g2;
	if(i==3) return _g3;

	return -1;
}

void  measure_pion_scattering_I2(double* m, int numsources, double precision,char* path,char* cnfg_filename, meson_observable** mo_arr){
	int ts;
	meson_observable *rho1[3][3],*rho2[3][3];
	meson_observable *pi1,*pi2;
	meson_observable *AD;
	meson_observable *BC;
	spinor_field* source_ts1 = alloc_spinor_field_f(4,&glattice);
	spinor_field* source_ts2 = alloc_spinor_field_f(4,&glattice);
	char auxname[256];

	spinor_field* prop_ts1 =  alloc_spinor_field_f(4 ,&glattice);
	spinor_field* prop_ts2=  alloc_spinor_field_f(4 ,&glattice);

	pi1 = (meson_observable*) malloc(sizeof(meson_observable));
	pi2 = (meson_observable*) malloc(sizeof(meson_observable));
	AD = (meson_observable*) malloc(sizeof(meson_observable));
	BC = (meson_observable*) malloc(sizeof(meson_observable));

	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			rho1[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
			rho2[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
		}
	}
	init_mo(pi1,"Pi1",GLB_T);
	init_mo(pi2,"Pi2",GLB_T);
	init_mo(AD,"AD",GLB_T);
	init_mo(BC,"BC",GLB_T);

	pi1->ind1 = _g5;
	pi1->ind2 = _g5;
	pi2->ind1 = _g5;
	pi2->ind2 = _g5;
		
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
				init_mo(rho1[i][j],"rho1",GLB_T);
				init_mo(rho2[i][j],"rho2",GLB_T);
				rho1[i][j]->ind1 = gi(i+1);
				rho1[i][j]->ind2 = gi(j+1);
				rho2[i][j]->ind1 = gi(i+1);
				rho2[i][j]->ind2 = gi(j+1);
		}
	}
	init_propagator_eo(1, m, precision);

	for (int src=0;src<numsources;++src)
   	{
		reset_mo(pi1);
		reset_mo(pi2);
		reset_mo(AD);
		reset_mo(BC);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				reset_mo(rho1[i][j]);
				reset_mo(rho2[i][j]);
		}
	}
	
	
	ts=random_tau();
	lprintf("MAIN",0,"ts = %d \n",ts);
	spinor_field_zero_f(source_ts1);
	spinor_field_zero_f(prop_ts1);
	create_diluted_source_equal_atau(source_ts1, ts);
	calc_propagator(prop_ts1,source_ts1,4);
	spinor_field_zero_f(source_ts2);
	create_diluted_source_equal_atau(source_ts2, ts);
	calc_propagator(prop_ts2,source_ts2,4);
	lprintf("MAIN",0,"Start to perform the contractions ...\n");
	
	// "standard" two points : pi and rho 
	measure_mesons_core(prop_ts1,prop_ts1,source_ts1,pi1,1,ts,1,0,GLB_T);
	measure_mesons_core(prop_ts2,prop_ts2,source_ts2,pi2,1,ts,1,0,GLB_T);
	do_global_sum(pi1,1.0);
	do_global_sum(pi2,1.0);

	for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				measure_mesons_core(prop_ts1,prop_ts1,source_ts1,rho1[i][j],1,ts,1,0,GLB_T);
				measure_mesons_core(prop_ts2,prop_ts2,source_ts2,rho2[i][j],1,ts,1,0,GLB_T);
				do_global_sum(rho1[i][j],1.0);
				do_global_sum(rho2[i][j],1.0);

			}
	}

	// contraction 4 particles two-points.
	measure_scattering_AD_core(AD, prop_ts1,prop_ts1,prop_ts2,prop_ts2, ts, 0,0,0,0,0); 
	measure_scattering_BC_core(BC, prop_ts1,prop_ts1,prop_ts2,prop_ts2, ts, 0,0,0,0,0);
	
	lprintf("MAIN",0,"Contraction done\n");
	// Printing.
	if (path!=NULL)
	{
		io2pt(pi1,1,src,path,"pi1",cnfg_filename);
		io2pt(pi2,1,src,path,"pi2",cnfg_filename);
		for(int i=0; i<3; i++){	
			for(int j=0; j<3; j++){
				sprintf(auxname, "rho1_%d%d",i+1,j+1);
				io2pt(rho1[i][j],1,src,path,auxname,cnfg_filename);
				sprintf(auxname, "rho2_%d%d",i+1,j+1);
				io2pt(rho2[i][j],1,src,path,auxname,cnfg_filename);
			}
		}
	
		io4pt(AD,0,src,path,"AD",cnfg_filename);
		io4pt(BC,0,src,path,"BC",cnfg_filename);
	}
	
	
	}
	if (mo_arr != NULL){
		copy_mo(mo_arr[0],AD);
		copy_mo(mo_arr[1],BC);
		copy_mo(mo_arr[2],pi1);
		copy_mo(mo_arr[3],pi2);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++)	{	
			copy_mo(mo_arr[4+j+3*i],rho1[i][j]);
			copy_mo(mo_arr[13+j+3*i],rho2[i][j]);
			}
		}
	}

	//free memory
	free_propagator_eo();
  	free_spinor_field_f(source_ts1);
  	free_spinor_field_f(source_ts2);
	free_spinor_field_f(prop_ts1);
  	free_spinor_field_f(prop_ts2);
	free_mo(pi1);
	free_mo(pi2);
	free_mo(AD);
	free_mo(BC);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			free_mo(rho1[i][j]);
			free_mo(rho2[i][j]);
		}
	}

}


double complex Ralt_contract(spinor_field* prop1,spinor_field* prop2,int tref ){
	int ix0,ix1,ix2,ix3,tc,i;
	double complex z=0.0;
	suNf_spinor* sptr1; 
	suNf_spinor* sptr2;
	
	for (ix0=0;ix0<T;ix0++){
		tc = (zerocoord[0] + ix0 + GLB_T ) % GLB_T;
		if(tc==tref){
			for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
	           i=ipt(ix0,ix1,ix2,ix3);
	           sptr1 = _FIELD_AT(prop1,i);
			   sptr2 = _FIELD_AT(prop2,i);
			   _spinor_prod_assign_f(z,*sptr1,*sptr2);
			}
		}
	}

	return z;
}

void measure_pion_scattering_I0(double* m, int numsources, double precision,char* path,char* cnfg_filename, int seq_prop, meson_observable** mo_arr){
	int ts=0;
	double complex A,B;
	int ts_vec[numsources];
	double memsize=0.;
	char auxname[256];
	meson_observable *tmp_mo;
	meson_observable *rho1[3][3];
	meson_observable *pi1;
    meson_observable *D;
    meson_observable *C;
 	meson_observable *R;
	meson_observable *Ralt;
    meson_observable *V;
	meson_observable *disc;
	spinor_field** source_ts1; 
	spinor_field* source_ts2= alloc_spinor_field_f(4,&glattice);
	spinor_field*** prop_ts1;
	spinor_field* prop_ts2= alloc_spinor_field_f(4,&glattice);
	spinor_field* seq_0=NULL;
	spinor_field* seq_t=NULL;
	spinor_field* seq_source=NULL;	
	

	// a lot of memory can be saved when the sequential sources are used. 
 	long int spinor_field_size = 4*sizeof(suNf_spinor)*GLB_VOLUME; // 4 for the spin dilution.

	if (seq_prop == 0) memsize = ((GLB_T + GLB_T*numsources + 2 ) *spinor_field_size); // more memory could be save in that case !!!
	if (seq_prop == 1) memsize = ((GLB_T + GLB_T*numsources + 5 ) *spinor_field_size);  // invert sequencial sources
	
	lprintf("MAIN",0,"measure_pion_scattering_I0 can be memory intensive! \n");
	lprintf("MAIN",0,"Total memory used (by spinor_fields) is  %1.6g GB \n",(double) memsize/1.e9);
	lprintf("MAIN",0,"Memory used (by spinor_fields) per task is %1.6g GB\n",(double) memsize/(1.e9*MPI_WORLD_SIZE));

	source_ts1= (spinor_field**)malloc(sizeof(spinor_field*)*GLB_T );
	for(int t=0; t<GLB_T; t++)  source_ts1[t] = alloc_spinor_field_f(4,&glattice);

	if (seq_prop){ // if seq_prop==1 
		seq_0= alloc_spinor_field_f(4,&glattice);
		seq_t= alloc_spinor_field_f(4,&glattice);
		seq_source =  alloc_spinor_field_f(4,&glattice);
	}

	prop_ts1= (spinor_field***)malloc(sizeof(spinor_field**)*numsources);
	for(int src=0; src<numsources; src++) prop_ts1[src] = (spinor_field**)malloc(sizeof(spinor_field*)*GLB_T);
	for(int src=0; src<numsources; src++) for(int t=0; t<GLB_T; t++)  prop_ts1[src][t] = alloc_spinor_field_f(4,&glattice);


	pi1 = (meson_observable*) malloc(sizeof(meson_observable));
	D = (meson_observable*) malloc(sizeof(meson_observable));
	C = (meson_observable*) malloc(sizeof(meson_observable));
	R = (meson_observable*) malloc(sizeof(meson_observable));
	Ralt= (meson_observable*) malloc(sizeof(meson_observable));
	tmp_mo = (meson_observable*) malloc(sizeof(meson_observable));
	V = (meson_observable*) malloc(sizeof(meson_observable));
	disc = (meson_observable*) malloc(sizeof(meson_observable)); // just for the g5 loop

	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			rho1[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
		}
	}
	init_mo(pi1,"Pi1",GLB_T);
	init_mo(D,"D",GLB_T);
	init_mo(C,"C",GLB_T);
	init_mo(R,"R",GLB_T);
	init_mo(Ralt,"Ralt",GLB_T);
	init_mo(tmp_mo,"tmp_mo",GLB_T);
	init_mo(V,"V",GLB_T);
	init_mo(disc,"disc",GLB_T);
	disc->ind1 = _g5;
	disc->ind2 = _NOGAMMA;
	pi1->ind1 = _g5;
	pi1->ind2 = _g5;
			
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
				init_mo(rho1[i][j],"rho1",GLB_T);
				rho1[i][j]->ind1 = gi(i+1);
				rho1[i][j]->ind2 = gi(j+1);
		}
	}
	init_propagator_eo(1, m, precision);
	for (int src=0;src<numsources;++src)
   	{
		reset_mo(pi1);
		reset_mo(D);
		reset_mo(C);
		reset_mo(R);
		reset_mo(Ralt);
		reset_mo(V);
		reset_mo(disc);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				reset_mo(rho1[i][j]);
			}
		}
		
		
		//ts=random_tau();
		ts = random_tau2(); // randomly chosen timeslice 
		ts_vec[src] = ts;

		lprintf("MAIN",0,"Random timeslice chosen : ts=%d\n",ts);	
		for (int t=0;t<GLB_T;++t){
			create_diluted_source_equal_atau(source_ts1[t], (t+ts)%GLB_T) ;
			calc_propagator(prop_ts1[src][t],source_ts1[t],4);
		}
		
		create_diluted_source_equal_atau(source_ts2, ts);
		calc_propagator(prop_ts2,source_ts2,4);
		if (seq_prop==1){ // if seq_prop= true
			lprintf("MAIN",0,"Creating sequential source...\n");	
			create_sequential_source_stoch(seq_source,ts,prop_ts1[src][0]);
			calc_propagator(seq_0,seq_source,4);
			lprintf("MAIN",0,"Sequential source inverted.\n");

		
 			lprintf("MAIN",0,"Contraction R...\n");
			for (int t=0;t<GLB_T;++t){
				create_sequential_source_stoch(seq_source,(ts+t)%GLB_T,prop_ts1[src][0]);
				calc_propagator(seq_t,seq_source,4);
				
				reset_mo(tmp_mo);
				measure_mesons_core(seq_0,seq_t,source_ts1[0],tmp_mo,1,ts,1,0,GLB_T);
				do_global_sum(tmp_mo,1.0); //

				R->corr_re[corr_ind(0,0,0,1,t,1,0)] = tmp_mo -> corr_re[corr_ind(0,0,0,1,t,1,0)];
				R->corr_im[corr_ind(0,0,0,1,t,1,0)] = tmp_mo -> corr_im[corr_ind(0,0,0,1,t,1,0)];
			}
		}
		//  V
		lprintf("MAIN",0,"Contraction V...\n");
		for (int t=0;t<GLB_T;++t){
			reset_mo(tmp_mo);
			measure_mesons_core(prop_ts1[src][t],prop_ts1[src][t],source_ts1[t],tmp_mo,1,(ts+t)%GLB_T,1,0,GLB_T); 
			do_global_sum(tmp_mo,1.0); //
			V->corr_re[corr_ind(0,0,0,1,(ts+t)%GLB_T,1,0)] = tmp_mo -> corr_re[corr_ind(0,0,0,1,0,1,0)];
			V->corr_im[corr_ind(0,0,0,1,(ts+t)%GLB_T,1,0)] = tmp_mo -> corr_im[corr_ind(0,0,0,1,0,1,0)];
		}
		// use the opportunity to compute disconnected loops. add all the other loops ? add test in check_scattering_length_I0
		for (int t=0;t<GLB_T;++t){
				measure_mesons_core(prop_ts1[src][t],prop_ts1[src][t],source_ts1[t],disc,1,0,1,0,GLB_T);	 
		}
		do_global_sum(disc,1.0); //

		lprintf("MAIN",0,"Contraction pi,rho ...\n");
		// "standard" two points : pi and rho 
		for (int t=0;t<GLB_T;++t){
			measure_mesons_core(prop_ts1[src][t],prop_ts1[src][t],source_ts1[t],pi1,1,(ts+t)%GLB_T,1,0,GLB_T); // this is summing over t
		}
		measure_mesons_core(prop_ts2,prop_ts2,source_ts2,pi1,1,ts,1,0,GLB_T);
		do_global_sum(pi1,1.0/(double)(GLB_T+1)); // this is averaging over GLB_T+1 timeslices !
		
		
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				for (int t=0;t<GLB_T;++t){
					measure_mesons_core(prop_ts1[src][t],prop_ts1[src][t],source_ts1[t],rho1[i][j],1,(ts+t)%GLB_T,1,0,GLB_T);
				}
				measure_mesons_core(prop_ts2,prop_ts2,source_ts2,rho1[i][j],1,ts,1,0,GLB_T);
				do_global_sum(rho1[i][j],1.0/(double)(GLB_T+1));
			}
		}

		lprintf("MAIN",0,"Contraction D,C...\n");
		// contraction 4 particles two-points.
		measure_scattering_AD_core(D, prop_ts1[src][0],prop_ts1[src][0],prop_ts2,prop_ts2, ts, 0,0,0,0,0); 
		measure_scattering_BC_core(C, prop_ts1[src][0],prop_ts1[src][0],prop_ts2,prop_ts2, ts, 0,0,0,0,0);
		

		lprintf("MAIN",0,"Contraction done\n");
		// Printing.
		if (path!=NULL) {
			io2pt(pi1,1,src,path,"pi1",cnfg_filename);
			io2pt(V,1,src,path,"V",cnfg_filename);
			io2pt(disc,1,src,path,"disc",cnfg_filename);

			for(int i=0; i<3; i++){	
				for(int j=0; j<3; j++){
					sprintf(auxname, "rho1_%d%d",i+1,j+1);
					io2pt(rho1[i][j],1,src,path,auxname,cnfg_filename);
				}
			}
			io4pt(D,0,src,path,"D",cnfg_filename);
			io4pt(C,0,src,path,"C",cnfg_filename);
		 	if (seq_prop==1) io4pt(R,0,src,path,"R",cnfg_filename); // if seq_prop= true
		}
	}
	
    if (mo_arr != NULL){
		copy_mo(mo_arr[0],D);
		copy_mo(mo_arr[1],C);
		copy_mo(mo_arr[2],R);
		copy_mo(mo_arr[3],V);
		copy_mo(mo_arr[4],pi1);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++)	{	
			copy_mo(mo_arr[5+j+3*i],rho1[i][j]);
			}
		}
		copy_mo(mo_arr[15],disc);
	}
	if (seq_prop != 2) //  if seq_prop==2,  neither R or Ralt are computed. Ralt is expensive for large number of sources.
	{
		lprintf("MAIN",0,"Contraction R alternative...\n");
	

		reset_mo(Ralt);
		int ts1,ts2;
		for (int src1=0;src1<numsources;++src1) for (int src2=0;src2<numsources;++src2){
			ts1 = ts_vec[src1];
			ts2 = ts_vec[src2];
			for (int t=0;t<GLB_T;++t) for (int dt=0;dt<GLB_T;++dt){
				for (int alpha=0;alpha<4;alpha++) for (int beta=0;beta<4;beta++) { // dilution indices 				
					A = Ralt_contract(prop_ts1[src1][t]+alpha,prop_ts1[src2][(t+dt)%GLB_T]+beta,(t+ts1)%GLB_T );
					B = Ralt_contract(prop_ts1[src2][(t+dt)%GLB_T]+beta,prop_ts1[src1][t]+alpha,(t+dt+ts2)%GLB_T );
					
					global_sum((double *)(&A),2);
					global_sum((double *)(&B),2);

					Ralt->corr_re[corr_ind(0,0,0,1,(dt+ ts2-ts1 + GLB_T)%GLB_T,1,0)] += creal(A*B)/(GLB_T*numsources*numsources);
					Ralt->corr_im[corr_ind(0,0,0,1,(dt+ ts2-ts1 + GLB_T)%GLB_T,1,0)] += cimag(A*B)/(GLB_T*numsources*numsources);

				}
			}
				
		}
	
	
	if (mo_arr != NULL) copy_mo(mo_arr[14],Ralt);
		
	if (path!=NULL) io4pt(Ralt,0,0,path,"Ralt",cnfg_filename);
	}

	//free memory
	free_propagator_eo();
	for (int t=0;t<GLB_T;++t){
	  	free_spinor_field_f(source_ts1[t]);
		for (int src=0;src<numsources;++src) free_spinor_field_f(prop_ts1[src][t]);
	}
  	free_spinor_field_f(source_ts2);
  	free_spinor_field_f(prop_ts2);
	if (seq_prop){ 
		free_spinor_field_f(seq_0);
		free_spinor_field_f(seq_t);
		free_spinor_field_f(seq_source);
	}

	free_mo(pi1);
	free_mo(D);
	free_mo(disc);
	free_mo(C);
	free_mo(R);
	free_mo(Ralt);
	free_mo(V);
	free_mo(tmp_mo);

	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			free_mo(rho1[i][j]);
		}
	}

}


void measure_pion_scattering_I0_TS(double* m, int numsources, double precision,char* path,char* cnfg_filename, int seq_prop, meson_observable** mo_arr){
        int ts=0;
        int ts_vec[numsources];
        meson_observable *sigmaconn, *sigmadisc,*Tr;


        spinor_field* source_disc= alloc_spinor_field_f(4,&glattice);
        spinor_field* source_tri= alloc_spinor_field_f(4,&glattice);
        spinor_field* source_seq= alloc_spinor_field_f(4,&glattice);

        spinor_field* prop_disc= alloc_spinor_field_f(4,&glattice);
        spinor_field* prop_tri= alloc_spinor_field_f(4,&glattice);
        spinor_field* seq_0 = alloc_spinor_field_f(4,&glattice);


        sigmaconn = (meson_observable*) malloc(sizeof(meson_observable));
        sigmadisc = (meson_observable*) malloc(sizeof(meson_observable));
        Tr  = (meson_observable*) malloc(sizeof(meson_observable));


        init_mo(sigmaconn,"sigmaconn",GLB_T);
        init_mo(sigmadisc,"sigmadisc",GLB_T);
        init_mo(Tr,"T",GLB_T);

        sigmaconn->ind1 = _id;
        sigmadisc->ind1 = _id;
        sigmaconn->ind2 = _id;
        sigmadisc->ind2 = _id;


        Tr->ind2 = _id;
        Tr->ind1 = _g5;



        for (int src=0;src<numsources;++src)
        {
                reset_mo(sigmaconn);
                reset_mo(sigmadisc);
                reset_mo(Tr);
                
                init_propagator_eo(1, m, precision);
                ts = random_tau2(); // randomly chosen timeslice 
                ts_vec[src] = ts;

                lprintf("MAIN",0,"Random timeslice chosen : ts=%d\n",ts);       

                // connected sigma and triangle
                create_diluted_source_equal_atau(source_tri, ts) ;
                calc_propagator(prop_tri,source_tri,4);
                create_sequential_source_stoch(source_seq,ts,prop_tri);
                calc_propagator(seq_0,source_seq,4);

                // disconected sigma
                create_noise_source_equal_eo(source_disc);
                calc_propagator(prop_disc,source_disc,4);


                // contractions
                measure_mesons_core(prop_disc,prop_disc,source_disc,sigmadisc,1,0,1,0,GLB_T); //disc
                measure_mesons_core(prop_tri,prop_tri,source_tri,sigmaconn,1,ts,1,0,GLB_T); //conn
                measure_mesons_core(prop_tri, seq_0,source_tri,Tr,1,ts,1,0,GLB_T); //triangle

                
                lprintf("MAIN",0,"Contractions...\n");
                do_global_sum(sigmadisc,1.0);
                do_global_sum(sigmaconn,0.5); //this is necessary for the factor 2 in io2pt
                do_global_sum(Tr,1.0);


                if (path!=NULL){
                io2pt( sigmaconn ,1,src,path,"sigmaconn",cnfg_filename);
                io2pt( sigmadisc ,1,src,path,"sigmadisc",cnfg_filename);
                io2pt( Tr ,1,src,path,"T",cnfg_filename);
                }

                
        }

		if (mo_arr != NULL){
			copy_mo(mo_arr[0],sigmaconn);
			copy_mo(mo_arr[1],sigmadisc);
			copy_mo(mo_arr[2],Tr);
		}
       
        //free memory
        free_spinor_field_f(source_disc);
        free_spinor_field_f(source_tri);
        free_spinor_field_f(source_seq);
        free_spinor_field_f(prop_disc);
        free_spinor_field_f(prop_tri);
        free_spinor_field_f(seq_0);

        free_mo(sigmaconn);
        free_mo(sigmadisc);
        free_mo(Tr);

}
