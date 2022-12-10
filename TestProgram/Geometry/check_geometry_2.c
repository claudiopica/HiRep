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
#include "linear_algebra.h"

#define _print_mat(a) lprintf("TEST",1,"(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",creal((a).c[0]),creal((a).c[1]),creal((a).c[2]),creal((a).c[3]),creal((a).c[4]),creal((a).c[5]),creal((a).c[6]),creal((a).c[7]),creal((a).c[8]), \
  cimag((a).c[0]), cimag((a).c[1]), cimag((a).c[2]), cimag((a).c[3]), cimag((a).c[4]), cimag((a).c[5]), cimag((a).c[6]), cimag((a).c[7]), cimag((a).c[8]))

int isNotZero(double zero, double tolerance) {
    return (fabs(zero)>tolerance);
}

int testUfieldIsUnit(int ix, int dir){
    suNg unity;
    double d;
    suNg *U=_4FIELD_AT(u_gauge,ix,dir);
    _suNg_unit(unity);
    _suNg_sub_assign(unity,*U);
    _suNg_sqnorm(d,unity);
    if(isNotZero(d,1.e-10)) {
        // _print_mat(*U);
        return 1;
    }
    return 0;

}

int main(int argc, char *argv[])
{
    int errors=0;
    setup_process(&argc, &argv);
    atexit(&finalize_process);

    setup_gauge_fields();

    lprintf("test",1,"GLATTICE\n");
    print_gd(&glattice);

    unit_u(u_gauge);
    complete_gf_sendrecv(u_gauge);

    //test plauette value
    double p=avr_plaquette();
    lprintf("TEST",1,"Plaquette= %lf\n",p);
    if (isNotZero(1.-p, 1.e-10)) {
        errors++;
        lprintf("ERROR",1,"Plaquette value is wrong\n");
    }

    // T dir
    if (NP_T>1) {
        for (int x=0; x<X; x++) 
        for (int y=0; y<Y; y++) 
        for (int z=0; z<Z; z++) { 
            int ix=0,iy,iz;
            ix=ipt(T,x,y,z);
            iy=ipt(T-1,x,y,z);
            iz=iup(iy,0);
            int err = testUfieldIsUnit(ix,0);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",T,x,y,z,ix,0);
            if(ix!=iz) lprintf("ERROR",1,"ipt and iup incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",T,x,y,z,ix,iz,0);
            errors += err;

            ix=ipt(-1,x,y,z);
            iy=ipt(0,x,y,z);
            iz=idn(iy,0);
            err = testUfieldIsUnit(ix,0);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",-1,x,y,z,ix,0);
            if(ix!=iz) lprintf("ERROR",1,"ipt and idn incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",-1,x,y,z,ix,iz,0);
            errors += err;
        }
    }

    // X dir
    if (NP_X>1) {
        for (int t=0; t<T; t++) 
        for (int y=0; y<Y; y++) 
        for (int z=0; z<Z; z++) { 
            int ix=0,iy,iz;
            ix=ipt(t,X,y,z);
            iy=ipt(t,X-1,y,z);
            iz=iup(iy,1);
            int err = testUfieldIsUnit(ix,1);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",t,X,y,z,ix,1);
            if(ix!=iz) lprintf("ERROR",1,"ipt and iup incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",t,X,y,z,ix,iz,1);
            errors += err;

            ix=ipt(t,-1,y,z);
            iy=ipt(t,0,y,z);
            iz=idn(iy,1);
            err = testUfieldIsUnit(ix,1);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",t,-1,y,z,ix,1);
            if(ix!=iz) lprintf("ERROR",1,"ipt and idn incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",t,-1,y,z,ix,iz,1);
            errors += err;
        }
    }

    // Y dir
    if (NP_Y>1) {
        for (int t=0; t<T; t++) 
        for (int x=0; x<X; x++) 
        for (int z=0; z<Z; z++) { 
            int ix=0,iy,iz;
            ix=ipt(t,x,Y,z);
            iy=ipt(t,x,Y-1,z);
            iz=iup(iy,2);
            int err = testUfieldIsUnit(ix,2);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",t,x,Y,z,ix,2);
            if(ix!=iz) lprintf("ERROR",1,"ipt and iup incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",t,x,Y,z,ix,iz,2);
            errors += err;

            ix=ipt(t,x,-1,z);
            iy=ipt(t,x,0,z);
            iz=idn(iy,2);
            err = testUfieldIsUnit(ix,2);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",t,x,-1,z,ix,2);
            if(ix!=iz) lprintf("ERROR",1,"ipt and idn incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",t,x,-1,z,ix,iz,2);
            errors += err;
        }
    }

    // Z dir
    if (NP_Z>1) {
        for (int t=0; t<T; t++) 
        for (int x=0; x<X; x++) 
        for (int y=0; y<Y; y++) { 
            int ix=0,iy,iz;
            ix=ipt(t,x,y,Z);
            iy=ipt(t,x,y,Z-1);
            iz=iup(iy,3);
            int err = testUfieldIsUnit(ix,3);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",t,x,y,Z,ix,3);
            if(ix!=iz) lprintf("ERROR",1,"ipt and iup incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",t,x,y,Z,ix,iz,3);
            errors += err;

            ix=ipt(t,x,y,-1);
            iy=ipt(t,x,y,0);
            iz=idn(iy,3);
            err = testUfieldIsUnit(ix,3);
            if (err) lprintf("ERROR",1,"U not unit @ [%d,%d,%d,%d] ix=%d dir=%d\n",t,x,y,-1,ix,3);
            if(ix!=iz) lprintf("ERROR",1,"ipt and idn incompatible: [%d,%d,%d,%d] ix=%d iz=%d dir=%d\n",t,x,y,-1,ix,iz,3);
            errors += err;
        }
    }

    lprintf("TEST",1,"Errors = %d\n",errors);
    return errors;
}
