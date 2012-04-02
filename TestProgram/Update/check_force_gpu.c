/******************************************************************************
*
* Test of modules
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
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
#include "gpu.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

#define _F_DIR0(u,chi1,chi2)				      \
_vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
_suNf_multiply(p.c[0],*(pu_gauge_f(x,0)),ptmp);		      \
_vector_add_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
_suNf_multiply(p.c[1],*(pu_gauge_f(x,0)),ptmp);		      \
_vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
_vector_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
_suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
_vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
_suNf_multiply(p.c[0],*(pu_gauge_f(x,1)),ptmp);		      \
_vector_i_add_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
_suNf_multiply(p.c[1],*(pu_gauge_f(x,1)),ptmp);		      \
_vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
_vector_i_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
_suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
_vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
_suNf_multiply(p.c[0],*(pu_gauge_f(x,2)),ptmp);		      \
_vector_sub_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
_suNf_multiply(p.c[1],*(pu_gauge_f(x,2)),ptmp);		      \
_vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
_vector_add_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
_suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
_vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
_suNf_multiply(p.c[0],*(pu_gauge_f(x,3)),ptmp);		      \
_vector_i_sub_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
_suNf_multiply(p.c[1],*(pu_gauge_f(x,3)),ptmp);		      \
_vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
_vector_i_add_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
_suNf_FMAT((u),p)


static void flip_mom(suNg_av_field  *momenta)
{
    
    suNg_algebra_vector  *dptr;
    
    _DECLARE_INT_ITERATOR(ix);
    
    geometry_descriptor *gd=momenta->type;
    
    int dx;
    
    _MASTER_FOR(gd,ix) {
        for(dx=0;dx <4 ; dx++)
        {
            dptr=(suNg_algebra_vector*)(momenta->ptr+4*ix+dx);
            
            _algebra_vector_mul_g(*dptr,-1.0,*dptr);
            
        }
    }
}

static double mom_sq(suNg_av_field  *momenta)
{
    
    suNg_algebra_vector  *dptr;
    
    _DECLARE_INT_ITERATOR(ix);
    
    geometry_descriptor *gd=momenta->type;
    
    int dx;
    double norm2=0.;
    
    _MASTER_FOR(gd,ix) {
        for(dx=0;dx <4 ; dx++)
        {
            double tmp;
            dptr=(suNg_algebra_vector*)(momenta->ptr+4*ix+dx);
            _algebra_vector_sqnorm_g(tmp,*dptr);
            norm2+=tmp;
        }
    }
    
    return norm2;
}



int main(int argc,char *argv[])
{
    _DECLARE_INT_ITERATOR(x);
    static suNg_algebra_vector f;
    static suNf_vector ptmp;
    static suNf_spinor p;

   char pame[256];
   int i;
    
      static suNf_FMAT s1;
   spinor_field *Xs, *Ys;
    suNg_av_field *force;
    
   float elapsed;
   gpu_timer t1;

   force_hmc_par par;
    double dt;
   
   /* setup process id and communications */
   setup_process(&argc,&argv);
   
   /* logger setup */
   logger_setlevel(0,10000); /* log all */
   //   logger_setlevel(0,10); /* log all */
   logger_map("DEBUG","debug");
#ifdef WITH_MPI
   sprintf(pame,">out_%d",PID); logger_stdout(pame);
   sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
#endif
   
   lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
   
   /* read input file */
   read_input(glb_var.read,"test_input");
   rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed);
   
   
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
   
   u_gauge=alloc_gfield(&glattice);
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
   u_gauge_f=alloc_gfield_f(&glattice);
#endif

   u_gauge_flt=alloc_gfield_flt(&glattice);
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
   u_gauge_f_flt=alloc_gfield_f_flt(&glattice);
#endif

   
   lprintf("MAIN",0,"Generating a random gauge field... ");
   fflush(stdout);
   random_u(u_gauge);
   gfield_copy_to_gpu(u_gauge);
   assign_ud2u_cpu();
   gfield_copy_to_gpu_flt(u_gauge_flt);
   lprintf("MAIN",0,"done.\n");
   
   lprintf("MAIN",0,"Representing gauge field... ");
   represent_gauge_field();
   //gfield_copy_to_gpu_f(u_gauge_f);
   represent_gauge_field_flt();
   //gfield_copy_to_gpu_f_flt(u_gauge_f_flt);
   lprintf("MAIN",0,"done.\n");

   Xs=alloc_spinor_field_f(2,&glattice);
   Ys=Xs+1;
    force = alloc_suNg_av_field(&glattice);
    
   gaussian_spinor_field(Xs);
   spinor_field_copy_to_gpu_f(Xs);
    gaussian_spinor_field(Ys);
    spinor_field_copy_to_gpu_f(Ys);
    zero_momenta(force);

    lprintf("TEST",0,"Init: Mom sq=%e\n",mom_sq(force));
    
    
    par.hasenbusch=0;
    dt=0.1;
   
    _PIECE_FOR(&glattice,x) { 
        _SITE_FOR(&glattice,x) {
            
            for (int mu=0; mu<4; ++mu) {
                int y;
                suNf_spinor *chi1, *chi2;
                _suNf_FMAT_zero(s1);
                switch (mu) {
                    case 0:
                        y=iup(x,0);
                        chi1=_FIELD_AT(Xs,x);
                        chi2=_FIELD_AT(Ys,y);
                        _F_DIR0(s1,chi1,chi2);
                        chi1=_FIELD_AT(Ys,x);
                        chi2=_FIELD_AT(Xs,y);
                        _F_DIR0(s1,chi1,chi2);
                        break;
                    case 1:
                        y=iup(x,1);
                        chi1=_FIELD_AT(Xs,x);
                        chi2=_FIELD_AT(Ys,y);
                        _F_DIR1(s1,chi1,chi2);
                        chi1=_FIELD_AT(Ys,x);
                        chi2=_FIELD_AT(Xs,y);
                        _F_DIR1(s1,chi1,chi2);
                        break;
                    case 2:
                        y=iup(x,2);
                        chi1=_FIELD_AT(Xs,x);
                        chi2=_FIELD_AT(Ys,y);
                        _F_DIR2(s1,chi1,chi2);
                        chi1=_FIELD_AT(Ys,x);
                        chi2=_FIELD_AT(Xs,y);
                        _F_DIR2(s1,chi1,chi2);
                        break;
                    default: /* DIR 3 */
                        y=iup(x,3);
                        chi1=_FIELD_AT(Xs,x);
                        chi2=_FIELD_AT(Ys,y);
                        _F_DIR3(s1,chi1,chi2);
                        chi1=_FIELD_AT(Ys,x);
                        chi2=_FIELD_AT(Xs,y);
                        _F_DIR3(s1,chi1,chi2);
                }
                
                _algebra_project(f,s1);
    
#ifdef UPDATE_EO
                if(par.hasenbusch != 1) {
                    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-dt*(_REPR_NORM2/_FUND_NORM2),f);
                } else {
                    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-par.bD*dt*(_REPR_NORM2/_FUND_NORM2),f);
                }
#else
                if(par.hasenbusch != 1) {
                    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),dt*(_REPR_NORM2/_FUND_NORM2),f);	
                } else {
                    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),par.bD*dt*(_REPR_NORM2/_FUND_NORM2),f);
                }
#endif

            }
        }
    }

    lprintf("TEST",0,"Before flip: Mom sq=%e\n",mom_sq(force));
    
    flip_mom(force);

    lprintf("TEST",0,"After flip: Mom sq=%e\n",mom_sq(force));

    
    force_hmc_gpu(force,Xs,Ys,dt,&par);
    
    lprintf("TEST",0,"Mom sq=%e\n",mom_sq(force));

    
   /* TEST CG_M */

/*
    
    t1 = gpuTimerStart();   
   cgiters = cg_mshift(&par, M, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M.dbl(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("CG TEST",0,"test cg[%d] = %e (req. %e)\n",i,tau,par.err2);
   }
   lprintf("CG TEST",0,"time: %1.10gms\n",elapsed);

   lprintf("CG TEST",0,"\n\nTesting CG multishift with single precision acceleration\n");
   lprintf("CG TEST",0,"------------------------------------------------------------\n");
   
   t1 = gpuTimerStart();
   cgiters = cg_mshift_flt(&par, M, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M.dbl(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("CG TEST",0,"test cg[%d] = %e (req. %e)\n",i,tau,par.err2);
     spinor_field_zero_f(&res[i]);
   }
   lprintf("CG TEST",0,"time: %1.10gms\n",elapsed);


   lprintf("CG TEST",0,"\n\nTesting CG multishift version 2 with single precision acceleration\n");
   lprintf("CG TEST",0,"------------------------------------------------------------\n");

   t1 = gpuTimerStart();
   cgiters = cg_mshift_flt2(&par, M, s1, res);
   elapsed = gpuTimerStop(t1);
   lprintf("CG TEST",0,"Converged in %d iterations\n",cgiters);
   for(i=0;i<par.n;++i){
     M.dbl(s2,&res[i]);
     spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
     spinor_field_sub_assign_f(s2,s1);
     tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
     lprintf("CG TEST",0,"test cg[%d] = %e (req. %e)\n",i,tau,par.err2);
     spinor_field_zero_f(&res[i]);
   }
   lprintf("CG TEST",0,"time: %1.10gms\n",elapsed);

   lprintf("CG TEST",0,"DONE!\n");
 
 */

   free_spinor_field_f(Xs);
   finalize_process();

   exit(0);
}
