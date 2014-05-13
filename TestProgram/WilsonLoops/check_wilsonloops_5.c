/*******************************************************************************
*
* Check Wilson loops(r=1,t=1) = plaquette
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "geometry.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "linear_algebra.h"
#include "communications.h"
#include "observables.h"
#include "error.h"


extern struct {
  int c[3];
  int* path;
  int length;
  int** perm;
  int nperms;
  int nsteps;
} WL_path[256];



int main(int argc,char *argv[])
{
  setup_process(&argc,&argv);
  
  logger_setlevel(0,10000); /* log all */
  if (PID!=0) { logger_disable(); }
  logger_map("DEBUG","debug");
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 


  
  read_input(glb_var.read,"test_input");
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }
  
  geometry_mpi_eo();
  /* setup random numbers */
  read_input(rlx_var.read,"test_input");


  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"The lattice global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"The lattice borders are (%d,%d,%d,%d)\n",T_BORDER,X_BORDER,Y_BORDER,Z_BORDER);
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  
  u_gauge=alloc_gfield(&glattice);
  
  suNg_field* h_gauge=alloc_gfield(&glattice);
  suNg* poly=amalloc(sizeof(suNg)*X*Y*Z,ALIGN);
  
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  
  double st_plaq_ave;
  double st_plaq[3];
  
  for(int k=0;k<3;k++) {
    st_plaq[k]=0.;
    for(int t=0;t<T;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      int i=ipt(t,x,y,z);
      st_plaq[k] += plaq(i,k+1,0);
    }
    st_plaq[k] /= GLB_T*GLB_X*GLB_Y*GLB_Z;
  }
  global_sum(st_plaq,3);
  st_plaq_ave=(st_plaq[0]+st_plaq[1]+st_plaq[2])/3.;

  for(int k=0;k<3;k++)
    lprintf("MAIN",0,"st_plaq[%d] = %e\n",k,st_plaq[k]);
  
  WL_initialize();
  
  int c[3]={1,0,0};
  WL_load_path(c,1);
  
  WL_Hamiltonian_gauge(h_gauge,u_gauge);
  
  WL_broadcast_polyakov(poly,h_gauge);

  double** WL;
  WL=amalloc(sizeof(double*),ALIGN);
  WL[0]=amalloc(sizeof(double)*GLB_T,ALIGN);

  double** tmp;
  tmp=amalloc(sizeof(double*),ALIGN);
  tmp[0]=amalloc(sizeof(double)*GLB_T,ALIGN);
  
  for(int t=0;t<GLB_T;t++) WL[0][t]=0.;
  
  int counter=0;
  for(int p=0;p<WL_path[0].nperms;p++) {
    for(int w=0;w<8;w++) {
      int sign[3];
      sign[0]=(w%2==0)?1:-1;
      sign[1]=((w/2)%2==0)?1:-1;
      sign[2]=((w/4)%2==0)?1:-1;
      if(WL_path[0].c[0]==0 && sign[0]==-1) continue;
      if(WL_path[0].c[1]==0 && sign[1]==-1) continue;
      if(WL_path[0].c[2]==0 && sign[2]==-1) continue;
    
      WL_correlators(tmp,h_gauge,poly,WL_path[0].nsteps,WL_path[0].path,WL_path[0].length,WL_path[0].perm[p],sign);
      
      lprintf("MAIN",0,"sign=(%d,%d,%d) perm=(%d,%d,%d) WL = %e\n",sign[0],sign[1],sign[2],
        WL_path[0].perm[p][0],WL_path[0].perm[p][1],WL_path[0].perm[p][2],
        tmp[0][1]);
      
      for(int t=0;t<GLB_T;t++) WL[0][t]+=tmp[0][t];
      counter++;
    }
  }

  for(int t=0;t<GLB_T;t++) WL[0][t]/=counter;

  double err=fabs(WL[0][1]-st_plaq_ave);
  
  
  lprintf("MAIN",0,"Checking that Wilson loops(r=1,t=1) = plaquette\n");
  lprintf("MAIN",0,"\n");


  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",err);
  lprintf("MAIN",0,"(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN",0,"\n");
  if(err<1.e-15)
    lprintf("MAIN",0,"check_wilsonloops_5 ... OK\n");
  else
    lprintf("MAIN",0,"check_wilsonloops_5 ... FAILED\n");


  
  WL_free();
  
  free_gfield(u_gauge);
  free_gfield(h_gauge);
  afree(poly);

  afree(WL[0]);
  afree(WL);
  afree(tmp[0]);
  afree(tmp);

  finalize_process();
  exit(0);
}
