/***************************************************************************\
* Copyright (c) 2013, Ari Hietanen & Ulrik Ishoej Soendergaard              *   
* All rights reserved.                                                      * 
\***************************************************************************/
#include <stdlib.h>

#include "geometry.h" 
#include "global.h" 
#include "error.h"
#include "logger.h"
#include "inverters.h" 
#include "linear_algebra.h"
#include "memory.h"
#include "update.h"

static spinor_field *P=NULL;

void DDalphaAMG_setup(mshift_par *par, spinor_operator M, int N, int nu, int n_inv){ 

	spinor_field *v,*tmp;				    // v is to contain the locally orthonormalized near kernel vectors
												// P is to contain the locally orthonormalized near kernel vectors
if (P==NULL){
	P=alloc_spinor_field_f(N,&glattice);		// &glattice might be changed to &even	
}
	v=alloc_spinor_field_f(N+1,&glattice);
	tmp=v+N;
	
	for (int i=0;i<N;i++){
		gaussian_spinor_field(&v[i]);		
			for (int eta=1;eta<4;eta++){		// Could be optimized
				SAP_prec(eta,&cg_mshift,par, M, &v[i], tmp);
				spinor_field_copy_f(&v[i],tmp);
			}
	}	// Line 6
	
	for (int eta=1;eta<=n_inv;eta++){
	// Line 8: Reconstruct P from v-vectors
	
	
	}
	
	
}


// 1 V-cycle (C operation)
void DDalphaAMG(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){ 
	SAP_prec(3,&cg_mshift,par, M, in, out);
}

void DDalphaAMG_finalize(){
	free_spinor_field_f(P);
	P=NULL;
}