
/******************************************************************************
 *
 * File check11.c
 *
 * Consistency checks on the programs in the module linalg
 *
 * Author: luigi del debbio <luigi.del.debbio@ed.ac.uk>
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "linear_algebra.h"
#include "representation.h"
#include "global.h"
#include "logger.h"
#include "utils.h"
#include "io.h"

#include "communications.h"

#define MAX_ROTATE 50

static complex v[25];
static double EPSILON=1.e-12;
static spinor_field *ppk[5];

static int initr=0;
static const suNf_spinor s0={{{{{0.0f}}}}};
static suNf_spinor *psi;


static void alloc_ws_rotate(void)
{
  psi=calloc(MAX_ROTATE,sizeof(suNf_spinor));
 
  error((psi==NULL),1,"alloc_ws_rotate [linalg.c]",
	"Unable to allocate workspace");

  initr=1;
}

static void rotate_ptr(int n,spinor_field *pkk[],complex v[])
{
  if (initr==0)
    alloc_ws_rotate();
  
  error((n<1)||(n>MAX_ROTATE),1,"rotate [eva.c]",
        "Parameter n is out of range");
  
  for(int i=0; i<n; i++)
    error((*pkk)->type!=(*(pkk+i))->type,1,"not available", "Spinors don't match!");
  
#undef _OMP_PRAGMA //This doesn't works with multiple threads
#define _OMP_PRAGMA(s)
  
  _MASTER_FOR(pkk[0]->type,ix)
  {
    for (int k=0;k<n;k++)
    {
      suNf_spinor *pk=&(psi[k]);
      suNf_spinor *pj=_FIELD_AT(pkk[0],ix);
      complex *z=&v[k];
      
      _spinor_mulc_f(*pk,*z,*pj);
      
      for (int j=1;j<n;j++)
	    {
	      pj=_FIELD_AT(pkk[j],ix);
	      z+=n;
        
	      _spinor_mulc_add_assign_f(*pk,*z,*pj);
	    }
    }
    
    for (int k=0;k<n;k++)
      *_FIELD_AT(pkk[k],ix)=psi[k];
  }
}

static void project(spinor_field *pk,spinor_field *pl)
{
  complex sp;

  sp.re=-spinor_field_prod_re_f(pl,pk);
  sp.im=-spinor_field_prod_im_f(pl,pk);

  spinor_field_mulc_add_assign_f(pk,sp,pl);
}   

static double normalize(spinor_field *ps)
{
  double r,ri;

  r=spinor_field_sqnorm_f(ps);
  r=sqrt(r);
  error(r<EPSILON,1,"normalize [eva.c]","vector has vanishing norm");

  ri=1.0/r;
  spinor_field_mul_f(ps,ri,ps);

  return (double)(r);
}

static complex sp(spinor_field *pk,spinor_field *pl)
{

  complex *rpk,*rpl,z;
  
  double x=0.0;
  double y=0.0;
  
  _TWO_SPINORS_FOR_SUM(pk,pl,x,y) {
    for (int i=0;i<(4*NF);i++)
    {
      complex *rpk = (complex*)_SPINOR_PTR(pk) + i;
      complex *rpl = (complex*)_SPINOR_PTR(pl) + i;
      x+=(double)((*rpk).re*(*rpl).re+(*rpk).im*(*rpl).im);
      y+=(double)((*rpk).re*(*rpl).im-(*rpk).im*(*rpl).re);
      //rpk+=1; //?? why these increment
      //rpl+=1; //??
    }
  }
  
  z.re=x;
  z.im=y;
  
#ifdef WITH_MPI
  global_sum((double*)&z,2);
#endif
  
  return z;
}


int main(int argc,char *argv[])
{
  int i,j;
  double r;
  double rd,zsqd;
  double d,dmax;
  complex w;
  complex zd,wd;
  spinor_field *ws;
  spinor_field *pk,*pl;
  spinor_field *tmp;

  char pame[256];

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(pame,">out_%d",PID); logger_stdout(pame);
    sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
  }

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input(glb_var.read,"test_input");

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  /* setup random numbers */
  read_input(rlx_var.read,"test_input");
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */


  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  lprintf("CPTEST",0,"spinor gsize=%d\n",glattice.gsize_spinor);
  lprintf("CPTEST",0,"spinor nbuffers=%d\n",glattice.nbuffers_spinor);
  lprintf("CPTEST",0,"spinor ncopies=%d\n",glattice.ncopies_spinor);
  lprintf("CPTEST",0,"gauge gsize=%d\n",glattice.gsize_gauge);
  lprintf("CPTEST",0,"gauge nbuffers=%d\n",glattice.nbuffers_gauge);
  lprintf("CPTEST",0,"gauge ncopies=%d\n",glattice.ncopies_gauge);
  lprintf("CPTEST",0,"lmp=%d\n",glattice.local_master_pieces);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif
  start_gf_sendrecv(u_gauge);

  represent_gauge_field();

       
   
  lprintf("LA TEST",0,"Consistency of the programs in the module linalg\n");
  lprintf("LA TEST",0,"------------------------------------------------\n");   

   
  represent_gauge_field();

  tmp=alloc_spinor_field_f(1, &glattice);
  ws=alloc_spinor_field_f(10, &glattice);


  for (i=0;i<10;i++)
    gaussian_spinor_field(&ws[i]);

  dmax=0.0;
   
  for (i=0;i<10;i++)
    {
      pk=&ws[i];
      pl=&ws[9-i];
      w=sp(pk,pl);
      
      zd=spinor_field_prod_f(pk,pl);
      rd=spinor_field_sqnorm_f(pk)*spinor_field_sqnorm_f(pl);
      d=((zd.re-(double)w.re)*(zd.re-(double)w.re)+
         (zd.im-(double)w.im)*(zd.im-(double)w.im));
      d=sqrt(d/rd);
      if (d>dmax)
	dmax=d;

      rd=spinor_field_prod_re_f(pk,pl);
      d=fabs(zd.re/rd-1.0);
      if (d>dmax)
	dmax=d;

      zd=spinor_field_prod_f(pk,pk);
      rd=spinor_field_sqnorm_f(pk);
      
      d=fabs(zd.im/rd);
      if (d>dmax)
	dmax=d;

      d=fabs(zd.re/rd-1.0f);
      if (d>dmax)
	dmax=d;
    }
  lprintf("LA TEST",0,"Check of spinor_field_prod, spinor_field_prod_re\n");
  lprintf("LA TEST",0,"and spinor_field_sqnorm: %.2e\n\n",dmax);

   
  dmax=0.0;
  zd.re= 0.345;
  zd.im=-0.876;
  zsqd=zd.re*zd.re+zd.im*zd.im;
   
  for (i=0;i<9;i++)
    {
      pk=&ws[i];
      pl=&ws[i+1];
      
      wd=spinor_field_prod_f(pk,pl);
      rd=spinor_field_sqnorm_f(pk)+zsqd*spinor_field_sqnorm_f(pl)
	+2.0*(zd.re*wd.re-zd.im*wd.im);
      
      spinor_field_mulc_add_assign_f(pk,zd,pl);

      d=fabs(rd/spinor_field_sqnorm_f(pk)-1.0);
      if (d>dmax)
	dmax=d;
    }
  lprintf("LA TEST",0,"Consistency of spinor_prod, norm_square\n");
  lprintf("LA TEST",0,"and mulc_spinor_add: %.2e\n\n",dmax);
  for (i=0;i<10;i++)
    gaussian_spinor_field(&ws[i]);

  dmax=0.0;
   
  for (i=0;i<10;i++)
    {
      pk=&ws[i];
      
      if (i>0)
	{
	  pl=&ws[i-1];
	  project(pk,pl);
	  zd=spinor_field_prod_f(pk,pl);
         
	  d=(fabs(zd.re)+
	     fabs(zd.im))/
            sqrt(spinor_field_sqnorm_f(pk));
         
	  if (d>dmax)
            dmax=d;
	}
      
      normalize(pk);
      rd=spinor_field_sqnorm_f(pk);
      
      d=fabs(rd-1.0f);
      if (d>dmax)
	dmax=d;
    }
   
  lprintf("LA TEST",0,"Consistency of spinor_prod, norm_square,\n");
  lprintf("LA TEST",0,"normalize and project: %.2e\n\n",dmax);
     
   
  for (i=0;i<5;i++)
    {
      pk=&ws[i];
      pl=&ws[i+5];
       
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
       
      for (j=0;j<5;j++)
	{
	  v[5*i+j].re=0.1234f*(double)(i^2)-0.8976f*(double)(j);
	  v[5*i+j].im=0.2231f*(double)(i)+0.9922f*(double)(j^2);
	}
       
      ppk[i]=pl;
    }
   
  rotate_ptr(5,ppk,v);
  dmax=0.0;
   
  for (i=5;i<10;i++)
    {
      pk=&ws[i];
      
      for (j=0;j<5;j++)
	{
	  zd.re=-(double)v[5*j+(i-5)].re;
	  zd.im=-(double)v[5*j+(i-5)].im;
         
	  pl=&ws[j];
	  spinor_field_mulc_add_assign_f(pk,zd,pl);
	}
      
      rd=spinor_field_sqnorm_f(pk);
      
      d=fabs(rd);
      if (d>dmax)
	dmax=d;
    }
   
  dmax/=spinor_field_sqnorm_f(&ws[0]);
  dmax=sqrt(dmax);
   
  lprintf("LA TEST",0,"Consistency of mulc_spinor_add\n");
  lprintf("LA TEST",0,"and rotate: %.2e\n\n",dmax);
 

  dmax=0.0;
  
  for (i=0;i<5;i++)
    {
      pk=&ws[i];
      pl=&ws[9-i];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      spinor_field_g5_f(tmp,pk);
      spinor_field_g5_f(pk,tmp);
      
      zd.re=-1.0;
      zd.im=0.0;
      
      spinor_field_mulc_add_assign_f(pl,zd,pk);
      r=spinor_field_sqnorm_f(pl)/spinor_field_sqnorm_f(pk);
      d=sqrt(r);
      if (d>dmax)
	dmax=d;
      
      gaussian_spinor_field(pl);
      zd=spinor_field_prod_f(pk,pl);
      spinor_field_g5_f(pk,pk);
      spinor_field_g5_f(pl,pl);
      wd=spinor_field_prod_f(pk,pl);
      
      d=(fabs(zd.re-wd.re)+fabs(zd.im-wd.im))/
	(fabs(zd.re)+fabs(zd.im));
      if (d>dmax)
	dmax=d;
    }
  
  lprintf("LA TEST",0,"Check of spinor_field_g5_f: %.2e\n\n",dmax);
  
  
  dmax=0.0;
  
  for (i=0;i<5;i++)
    {
      pk=&ws[i];
      pl=&ws[9-i];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      d=-2.5;
      spinor_field_lc1_f(d,pk,pl);
      
      zd.re=1.5;
      zd.im=0.0;
      spinor_field_mulc_add_assign_f(pk,zd,pl);
      d=spinor_field_sqnorm_f(pk)/spinor_field_sqnorm_f(pl);
      
      if (d>dmax)
	dmax=d;
    }
  
  lprintf("LA TEST",0,"Check of lc1: %.2e\n\n",dmax);
  
  dmax=0.0;
  
  for (i=0;i<5;i++)
    {
      pk=&ws[i];
      pl=&ws[9-i];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      d=1.0;
      r=2.5;
      spinor_field_lc2_f(d,r,pk,pl);
      
      zd.re=-3.5;
      zd.im=0.0;
      spinor_field_mulc_add_assign_f(pk,zd,pl);
      d=spinor_field_sqnorm_f(pk)/spinor_field_sqnorm_f(pl);
      
      if (d>dmax)
	dmax=d;
    }
  
  lprintf("LA TEST",0,"Check of lc2: %.2e\n\n",dmax);
  
  dmax=0.0;
  
  for (i=0;i<5;i++)
    {
      pk=&ws[i];
      pl=&ws[9-i];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      d=3.5;
      r=-1.5;
      spinor_field_lc3_f(d,r,pk,pl,pk);
      
      zd.re=-1.0;
      zd.im=0.0;
      spinor_field_mulc_add_assign_f(pk,zd,pl);
      d=spinor_field_sqnorm_f(pk)/spinor_field_sqnorm_f(pl);
      
      if (d>dmax)
	dmax=d;
    }
  
  lprintf("LA TEST",0,"Check of lc3: %.2e\n\n",dmax);
  
  finalize_process();

  return 0;

}
