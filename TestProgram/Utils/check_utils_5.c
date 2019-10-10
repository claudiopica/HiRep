/*******************************************************************************
*
* Check of the shift field routine.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "utils.h"
#include "update.h"
#include "observables.h"
#include "random.h"
#include "logger.h"
#include "communications.h"
#include "representation.h"
#include "setup.h"
#include "memory.h"
#include "hr_complex.h"
#include "linear_algebra.h"

int main(int argc, char *argv[])
{

    int return_value = 0;
    double tmp[4];
    int shift[4];
    double plaq[2];
    spinor_field *ss, *ss_e, *ss_o;
    double complex sobs[7];

    setup_process(&argc, &argv);

    setup_gauge_fields();
    ss = alloc_spinor_field_f(2, &glattice);
    ss_e = alloc_spinor_field_f(2, &glat_even);
    ss_o = alloc_spinor_field_f(2, &glat_odd);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_gf_sendrecv(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the plaquette... ");
    plaq[0] = avr_plaquette();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Generating and applying a random even shift to the gauge field... ");
    ranlxd(tmp, 4);
    bcast(tmp, 4);

    shift[0] = 2 * (int)(tmp[0] * GLB_T);
    shift[1] = 2 * (int)(tmp[1] * GLB_X);
    shift[2] = 2 * (int)(tmp[2] * GLB_Y);
    shift[3] = 2 * (int)(tmp[3] * GLB_Z);

    shift_fields(shift, NULL, u_gauge, NULL, u_gauge);

    lprintf("MAIN", 0, "done.\n\n");
    lprintf("MAIN", 0, "Shift = (%d,%d,%d,%d)\n ", shift[0], shift[1], shift[2], shift[3]);

    lprintf("MAIN", 0, "Measuring the plaquette on the shifted field... ");
    plaq[1] = avr_plaquette();
    lprintf("MAIN", 0, "done.\n\n");

    plaq[1] -= plaq[0];

    sqrt(plaq[1] * plaq[1]);

    lprintf("MAIN", 0, "Checking invariance of the plaquette 0pp on shifted configurations.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", sqrt(plaq[1] * plaq[1]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (sqrt(plaq[1] * plaq[1]) > 10.e-14)
        return_value++;

    lprintf("MAIN", 0, "Generating two random spinor fields on the whole lattice... ");

    gaussian_spinor_field(ss);
    gaussian_spinor_field(ss + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the norms and dot products of the fields... ");
    sobs[0] = spinor_field_sqnorm_f(ss);
    sobs[1] = spinor_field_sqnorm_f(ss + 1);
    sobs[2] = spinor_field_prod_f(ss, ss + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Applying the random shift to the spinor fields... ");
    shift_fields(shift, ss, NULL, ss, NULL);
    sobs[3] = spinor_field_sqnorm_f(ss);
    sobs[4] = spinor_field_prod_f(ss, ss + 1);
    shift_fields(shift, ss + 1, NULL, ss + 1, NULL);
    sobs[5] = spinor_field_sqnorm_f(ss + 1);
    sobs[6] = spinor_field_prod_f(ss, ss + 1);

    lprintf("MAIN", 0, "done.\n\n");

    sobs[3] = (sobs[3] - sobs[0]) / GLB_VOLUME;
    sobs[5] = (sobs[5] - sobs[1]) / GLB_VOLUME;
    sobs[6] = (sobs[6] - sobs[2]) / GLB_VOLUME;
    sobs[4] = (sobs[4] - sobs[2]) / GLB_VOLUME;

    lprintf("MAIN", 0, "Checking invariance of all the combinations of scalar products before and after the shift.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[3]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[4]));
    lprintf("MAIN", 0, "(should be different from zero)\n\n");
    if (cabs(sobs[4]) < 10.e-4)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[5]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[6]) > 10.e-14)
        return_value++;

    lprintf("MAIN", 0, "Generating two random spinor fields on the even lattice... ");

    gaussian_spinor_field(ss_e);
    gaussian_spinor_field(ss_e + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the norms and dot products field... ");
    sobs[0] = spinor_field_sqnorm_f(ss_e);
    sobs[1] = spinor_field_sqnorm_f(ss_e + 1);
    sobs[2] = spinor_field_prod_f(ss_e, ss_e + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Applying the random shift to the spinor field... ");
    if ((shift[0] + shift[1] + shift[2] + shift[3]) % 2 == 0)
    {
        shift_fields(shift, ss_e, NULL, ss_e, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_e);
        sobs[4] = spinor_field_prod_f(ss_e, ss_e + 1);
        shift_fields(shift, ss_e + 1, NULL, ss_e + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_e + 1);
        sobs[6] = spinor_field_prod_f(ss_e, ss_e + 1);
    }
    else
    {
        shift_fields(shift, ss_e, NULL, ss_e, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_e);
        sobs[4] = spinor_field_prod_f(ss_e, ss_e + 1);
        shift_fields(shift, ss_e + 1, NULL, ss_e + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_e + 1);
        sobs[6] = spinor_field_prod_f(ss_e, ss_e + 1);
    }
    lprintf("MAIN", 0, "done.\n\n");

    sobs[3] = (sobs[3] - sobs[0]) / GLB_VOLUME;
    sobs[5] = (sobs[5] - sobs[1]) / GLB_VOLUME;
    sobs[6] = (sobs[6] - sobs[2]) / GLB_VOLUME;
    sobs[4] = (sobs[4] - sobs[2]) / GLB_VOLUME;

    lprintf("MAIN", 0, "Checking invariance of all the combinations of scalar products before and after the shift.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[3]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[4]));
    lprintf("MAIN", 0, "(should be different from zero)\n\n");
    if (cabs(sobs[4]) < 10.e-4)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[5]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[6]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Generating two random spinor fields on the odd lattice... ");

    gaussian_spinor_field(ss_o);
    gaussian_spinor_field(ss_o + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the norms and dot products field... ");
    sobs[0] = spinor_field_sqnorm_f(ss_o);
    sobs[1] = spinor_field_sqnorm_f(ss_o + 1);
    sobs[2] = spinor_field_prod_f(ss_o, ss_o + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Applying the random shift to the spinor field... ");
    if ((shift[0] + shift[1] + shift[2] + shift[3]) % 2 == 0)
    {
        shift_fields(shift, ss_o, NULL, ss_o, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_o);
        sobs[4] = spinor_field_prod_f(ss_o, ss_o + 1);
        shift_fields(shift, ss_o + 1, NULL, ss_o + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_o + 1);
        sobs[6] = spinor_field_prod_f(ss_o, ss_o + 1);
    }
    else
    {
        shift_fields(shift, ss_o, NULL, ss_o, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_o);
        sobs[4] = spinor_field_prod_f(ss_o, ss_o + 1);
        shift_fields(shift, ss_o + 1, NULL, ss_o + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_o + 1);
        sobs[6] = spinor_field_prod_f(ss_o, ss_o + 1);
    }
    lprintf("MAIN", 0, "done.\n\n");

    sobs[3] = (sobs[3] - sobs[0]) / GLB_VOLUME;
    sobs[5] = (sobs[5] - sobs[1]) / GLB_VOLUME;
    sobs[6] = (sobs[6] - sobs[2]) / GLB_VOLUME;
    sobs[4] = (sobs[4] - sobs[2]) / GLB_VOLUME;

    lprintf("MAIN", 0, "Checking invariance of all the combinations of scalar products before and after the shift.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[3]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[4]));
    lprintf("MAIN", 0, "(should be different from zero)\n\n");
    if (cabs(sobs[4]) < 10.e-4)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[5]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[6]) > 10.e-14)
        return_value++;

    lprintf("MAIN", 0, "Generating and applying a random odd shift to the gauge field... ");
    ranlxd(tmp, 4);
    bcast(tmp, 4);

    shift[0] = 2 * (int)(tmp[0] * GLB_T);
    shift[1] = 2 * (int)(tmp[1] * GLB_X) + 1;
    shift[2] = 2 * (int)(tmp[2] * GLB_Y);
    shift[3] = 2 * (int)(tmp[3] * GLB_Z);

    shift_fields(shift, NULL, u_gauge, NULL, u_gauge);

    lprintf("MAIN", 0, "done.\n\n");
    lprintf("MAIN", 0, "Shift = (%d,%d,%d,%d)\n ", shift[0], shift[1], shift[2], shift[3]);

    lprintf("MAIN", 0, "Measuring the plaquette on the shifted field... ");
    plaq[1] = avr_plaquette();
    lprintf("MAIN", 0, "done.\n\n");

    plaq[1] -= plaq[0];

    sqrt(plaq[1] * plaq[1]);

    lprintf("MAIN", 0, "Checking invariance of the plaquette 0pp on shifted configurations.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", sqrt(plaq[1] * plaq[1]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (sqrt(plaq[1] * plaq[1]) > 10.e-14)
        return_value++;

    lprintf("MAIN", 0, "Generating two random spinor fields on the whole lattice... ");

    gaussian_spinor_field(ss);
    gaussian_spinor_field(ss + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the norms and dot products of the fields... ");
    sobs[0] = spinor_field_sqnorm_f(ss);
    sobs[1] = spinor_field_sqnorm_f(ss + 1);
    sobs[2] = spinor_field_prod_f(ss, ss + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Applying the random shift to the spinor fields... ");
    shift_fields(shift, ss, NULL, ss, NULL);
    sobs[3] = spinor_field_sqnorm_f(ss);
    sobs[4] = spinor_field_prod_f(ss, ss + 1);
    shift_fields(shift, ss + 1, NULL, ss + 1, NULL);
    sobs[5] = spinor_field_sqnorm_f(ss + 1);
    sobs[6] = spinor_field_prod_f(ss, ss + 1);

    lprintf("MAIN", 0, "done.\n\n");

    sobs[3] = (sobs[3] - sobs[0]) / GLB_VOLUME;
    sobs[5] = (sobs[5] - sobs[1]) / GLB_VOLUME;
    sobs[6] = (sobs[6] - sobs[2]) / GLB_VOLUME;
    sobs[4] = (sobs[4] - sobs[2]) / GLB_VOLUME;

    lprintf("MAIN", 0, "Checking invariance of all the combinations of scalar products before and after the shift.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[3]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[4]));
    lprintf("MAIN", 0, "(should be different from zero)\n\n");
    if (cabs(sobs[4]) < 10.e-4)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[5]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[6]) > 10.e-14)
        return_value++;

    lprintf("MAIN", 0, "Generating two random spinor fields on the even lattice... ");

    gaussian_spinor_field(ss_e);
    gaussian_spinor_field(ss_e + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the norms and dot products field... ");
    sobs[0] = spinor_field_sqnorm_f(ss_e);
    sobs[1] = spinor_field_sqnorm_f(ss_e + 1);
    sobs[2] = spinor_field_prod_f(ss_e, ss_e + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Applying the random shift to the spinor field... ");
    if ((shift[0] + shift[1] + shift[2] + shift[3]) % 2 == 0)
    {
        shift_fields(shift, ss_e, NULL, ss_e, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_e);
        sobs[4] = spinor_field_prod_f(ss_e, ss_e + 1);
        shift_fields(shift, ss_e + 1, NULL, ss_e + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_e + 1);
        sobs[6] = spinor_field_prod_f(ss_e, ss_e + 1);
    }
    else
    {
        shift_fields(shift, ss_e, NULL, ss_o, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_o);
        sobs[4] = spinor_field_prod_f(ss_o, ss_o + 1);
        shift_fields(shift, ss_e + 1, NULL, ss_o + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_o + 1);
        sobs[6] = spinor_field_prod_f(ss_o, ss_o + 1);
    }
    lprintf("MAIN", 0, "done.\n\n");

    sobs[3] = (sobs[3] - sobs[0]) / GLB_VOLUME;
    sobs[5] = (sobs[5] - sobs[1]) / GLB_VOLUME;
    sobs[6] = (sobs[6] - sobs[2]) / GLB_VOLUME;
    sobs[4] = (sobs[4] - sobs[2]) / GLB_VOLUME;

    lprintf("MAIN", 0, "Checking invariance of all the combinations of scalar products before and after the shift.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[3]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[4]));
    lprintf("MAIN", 0, "(should be different from zero)\n\n");
    if (cabs(sobs[4]) < 10.e-4)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[5]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[6]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Generating two random spinor fields on the odd lattice... ");

    gaussian_spinor_field(ss_o);
    gaussian_spinor_field(ss_o + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Measuring the norms and dot products field... ");
    sobs[0] = spinor_field_sqnorm_f(ss_o);
    sobs[1] = spinor_field_sqnorm_f(ss_o + 1);
    sobs[2] = spinor_field_prod_f(ss_o, ss_o + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Applying the random shift to the spinor field... ");
    if ((shift[0] + shift[1] + shift[2] + shift[3]) % 2 == 0)
    {
        shift_fields(shift, ss_o, NULL, ss_o, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_o);
        sobs[4] = spinor_field_prod_f(ss_o, ss_o + 1);
        shift_fields(shift, ss_o + 1, NULL, ss_o + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_o + 1);
        sobs[6] = spinor_field_prod_f(ss_o, ss_o + 1);
    }
    else
    {
        shift_fields(shift, ss_o, NULL, ss_e, NULL);
        sobs[3] = spinor_field_sqnorm_f(ss_e);
        sobs[4] = spinor_field_prod_f(ss_e, ss_e + 1);
        shift_fields(shift, ss_o + 1, NULL, ss_e + 1, NULL);
        sobs[5] = spinor_field_sqnorm_f(ss_e + 1);
        sobs[6] = spinor_field_prod_f(ss_e, ss_e + 1);
    }
    lprintf("MAIN", 0, "done.\n\n");

    sobs[3] = (sobs[3] - sobs[0]) / GLB_VOLUME;
    sobs[5] = (sobs[5] - sobs[1]) / GLB_VOLUME;
    sobs[6] = (sobs[6] - sobs[2]) / GLB_VOLUME;
    sobs[4] = (sobs[4] - sobs[2]) / GLB_VOLUME;

    lprintf("MAIN", 0, "Checking invariance of all the combinations of scalar products before and after the shift.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[3]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[4]));
    lprintf("MAIN", 0, "(should be different from zero)\n\n");
    if (cabs(sobs[4]) < 10.e-4)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[5]) > 10.e-14)
        return_value++;
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", cabs(sobs[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(sobs[6]) > 10.e-14)
        return_value++;

    finalize_process();
    return return_value;
}
