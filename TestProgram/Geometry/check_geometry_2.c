/*******************************************************************************
*
* NOCOMPILE = !WITH_MPI
*
* Check communication of gauge field in T direction
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"
#include "setup.h"
#include "random.h"
#include "suN.h"

#define _print_mat(a) lprintf("TEST",1,"(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",creal((a).c[0]),creal((a).c[1]),creal((a).c[2]),creal((a).c[3]),creal((a).c[4]),creal((a).c[5]),creal((a).c[6]),creal((a).c[7]),creal((a).c[8]), \
  cimag((a).c[0]), cimag((a).c[1]), cimag((a).c[2]), cimag((a).c[3]), cimag((a).c[4]), cimag((a).c[5]), cimag((a).c[6]), cimag((a).c[7]), cimag((a).c[8]))


int main(int argc, char *argv[])
{
    setup_process(&argc, &argv);

    if (!(NP_T>1)) {
        lprintf("TEST",1,"This test only works with NP_T>1. Exiting now...\n");
        return 0;
    }

    setup_gauge_fields();

    unit_u(u_gauge);
    complete_gf_sendrecv(u_gauge);

    double p=avr_plaquette();
    lprintf("TEST",1,"Plaquette= %lf\n",p);


    suNg unity;
    
    for (int x=0; x<X; x++) 
    for (int y=0; y<Y; y++) 
    for (int z=0; z<Z; z++) { 
        int ix=0,iy,iz;
        double d=0.;
        suNg *U;

        ix=ipt(T,x,y,z);
        iy=ipt(T-1,x,y,z);
        iz=iup(iy,0);
        U=_4FIELD_AT(u_gauge,ix,0);
        _suNg_unit(unity);
        _suNg_sub_assign(unity,*U);
        _suNg_sqnorm(d,unity);
        if(d>1.e-10) {
            _print_mat(*U);
            lprintf("TEST",1,"Norm [%d,%d,%d,%d]=%g ix=%d iz=%d\n",T,x,y,z,d,ix,iz);
        }
        if(ix!=iz) lprintf("ERROR",1,"ipt and iup incompatible: [%d,%d,%d,%d] ix=%d iz=%d\n",T,x,y,z,ix,iz);

        ix=ipt(-1,x,y,z);
        iy=ipt(0,x,y,z);
        iz=idn(iy,0);
        U=_4FIELD_AT(u_gauge,ix,0);
        _suNg_unit(unity);
        _suNg_sub_assign(unity,*U);
        _suNg_sqnorm(d,unity);
        if(d>1.e-10) {
            _print_mat(*U);
            lprintf("TEST",1,"Norm [%d,%d,%d,%d]=%g ix=%d iz=%d\n",-1,x,y,z,d,ix,iz);
        }
        if(ix!=iz) lprintf("ERROR",1,"ipt and idn incompatible: [%d,%d,%d,%d] ix=%d iz=%d\n",-1,x,y,z,ix,iz);

    }

    /* close communications */
    finalize_process();

    return 0;
}
