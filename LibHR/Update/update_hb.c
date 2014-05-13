/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
 *
 * File update.c
 *
 * Update programs
 *
 *******************************************************************************/

#define PROJECT_INTERVAL 10

#include "suN.h"
#include "utils.h"
#include "global.h"
#include "update.h"
#include "communications.h"
#define PI 3.141592653589793238462643383279502884197
#include <math.h>
#include <stdlib.h>

static suNg v __attribute__ ((aligned (16)));
static int * dyn_gauge=NULL;

void project_gauge_field(void)
{
  _DECLARE_INT_ITERATOR(ix);

  _MASTER_FOR(&glattice,ix) {
    project_to_suNg(pu_gauge(ix,0));
    project_to_suNg(pu_gauge(ix,1));
    project_to_suNg(pu_gauge(ix,2));
    project_to_suNg(pu_gauge(ix,3));
  }
} 

#if defined(BASIC_SF) || defined(ROTATED_SF)
static void g_up_Dirichlet_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(T-1,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

#if defined(BASIC_SF) || defined(ROTATED_SF) || defined(BC_T_MIXED)
static void g_dn_Dirichlet_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(0,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(1,ix,iy,iz);
	  dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

#if defined(BC_T_OPEN) || defined(BC_T_MIXED)
static void g_up_open_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(T-1,ix,iy,iz);
	  dyn_gauge[index*4]=0;
	}
  }
}
#endif

#if defined(BC_T_OPEN) 
static void g_dn_open_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(0,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

static void init_hb_boundary() {
  dyn_gauge =(int*)malloc(sizeof(int)*glattice.gsize_gauge*4);
  int i;
  for(i=0;i<glattice.gsize_gauge*4;i++) dyn_gauge[i]=1;
#if defined(BASIC_SF) || defined(ROTATED_SF)
  g_up_Dirichlet_BCs();
  g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_MIXED
  g_up_open_BCs();
  g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_OPEN
  g_up_open_BCs();
  g_dn_open_BCs();
#endif
}






static void update_all(double beta,int type)
{
  int mu,i,j;
  static int count=PROJECT_INTERVAL;

  if (count>=PROJECT_INTERVAL) {
    project_gauge_field();
    count=0;
  }
  ++count;
   
  for(mu=0;mu<4;mu++){
    start_gf_sendrecv(u_gauge);
    for(j=glat_even.master_start[0];j<=glat_even.master_end[0];j++){
      if(dyn_gauge[j*4+mu]!=0){
      staples(j,mu,&v);
      cabmar(beta,pu_gauge(j,mu),&v,type);
      }
    }
    complete_gf_sendrecv(u_gauge);
    for(i=1;i<glat_even.local_master_pieces;i++)
      for(j=glat_even.master_start[i];j<=glat_even.master_end[i];j++){
	if(dyn_gauge[j*4+mu]!=0){
	  staples(j,mu,&v);
	  cabmar(beta,pu_gauge(j,mu),&v,type);
	}
      }
  }
   
  for(mu=0;mu<4;mu++){
    start_gf_sendrecv(u_gauge);
    for(j=glat_odd.master_start[0];j<=glat_odd.master_end[0];j++){
      if(dyn_gauge[j*4+mu]!=0){
	staples(j,mu,&v);
	cabmar(beta,pu_gauge(j,mu),&v,type);
      }
    }
    complete_gf_sendrecv(u_gauge);
    for(i=1;i<glat_odd.local_master_pieces;i++)
      for(j=glat_odd.master_start[i];j<=glat_odd.master_end[i];j++){
	if(dyn_gauge[j*4+mu]!=0){
	  staples(j,mu,&v);
	  cabmar(beta,pu_gauge(j,mu),&v,type);
	}
      }
  }

} 


void update(double beta,int nhb,int nor)
{
  if(dyn_gauge==NULL ) init_hb_boundary();
  int n;
   
  for (n=0;n<nhb;n++){     
    update_all(beta,0);
  }

  for (n=0;n<nor;n++){     
    update_all(beta,1);
  }

 
}

