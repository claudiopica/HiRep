/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File plaquette.c
*
* Routines for the average plaquette
*
*******************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"




void cplaq_sblk(double complex *ret,int ix,int mu,int nu)
{
  int iy,iz;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;
  double tmpre=0.;
  double tmpim=0.;
  
  iy=iup_sblk(ix,mu);
  iz=iup_sblk(ix,nu);

  v1=pu_gauge_sblk(ix,mu);
  v2=pu_gauge_sblk(iy,nu);
  v3=pu_gauge_sblk(iz,mu);
  v4=pu_gauge_sblk(ix,nu);

  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,(*v4),(*v3));
  _suNg_times_suNg_dagger(w3,w1,w2);

  _suNg_trace_re(tmpre,w3);

#ifndef GAUGE_SON
  _suNg_trace_im(tmpim,w3);
#endif
  *ret = tmpre + I*tmpim;
 
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    *ret *= plaq_weight[ix*16+mu*4+nu];
  }
#endif
   
}



double complex spatial_plaquette_sblk()
{
	static double complex pa;
	static double complex r0;
	
_OMP_PRAGMA ( single )
    {
      r0 = 0.;
    }

	_PIECE_FOR(&glattice,ixp)
	{
		
		_SITE_FOR_SUM(&glattice,ixp,ix,r0)
		{
			double complex tmp;

			cplaq_sblk(&tmp,ix,2,1);
			r0 += tmp;

			cplaq_sblk(&tmp,ix,3,1);
			r0 += tmp;

			cplaq_sblk(&tmp,ix,3,2);
			r0 += tmp;
		}
	}


_OMP_PRAGMA ( single )
    {
	pa=r0;

    }

	global_sum((double*)&pa,2);

#ifdef BC_T_OPEN
		pa[k] /= NG*GLB_VOLUME*(GLB_T-1)/GLB_T;
#else
		pa /= NG*GLB_VOLUME;
#endif
return pa;
}
