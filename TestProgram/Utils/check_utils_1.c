/*******************************************************************************
*
* Check of the workspace routines
*
*******************************************************************************/

#include "libhr.h"

#define nwrk 5

int main(int argc, char *argv[]) {
    int return_value = 0, i;
    suNg_field *wrk[nwrk + 1];
    int *iup_wrk[nwrk + 1];
    int *idn_wrk[nwrk + 1];
    int idx_wrk[nwrk + 1];
    hr_complex plaq[nwrk + 1], test;
    double dop;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Requesting and initializing %d  workspace gauge fields\n", nwrk);

    for (i = 0; i < nwrk; i++) {
        idx_wrk[i] = reserve_wrk_space_with_pointers(wrk + i, iup_wrk + i, idn_wrk + i);
        random_u(wrk[i]);
    }
    lprintf("MAIN", 0, "done.\n\n");

    for (i = 0; i < nwrk; i++) {
        set_wrk_space(i);
        plaq[i] = avr_plaquette_wrk();
    }

    lprintf("MAIN", 0, "Copying wrkfield %d on wrkfield %d\n", 0, nwrk - 2);

    copy_suNg_field(wrk[nwrk - 2], wrk[0]);

    set_wrk_space(idx_wrk[nwrk - 2]);
    test = avr_plaquette_wrk();
    test -= plaq[0];

    dop = sqrt(_complex_prod_re(test, test));

    lprintf("MAIN", 0, "Checking gauge plaquette of the wrkfield %d and wrkfield %d.\n ", 0, nwrk - 2);
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
    if (dop > 10.e-14) { return_value++; }
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Releasing wrkfield %d\n", nwrk - 2);
    release_wrk_space(idx_wrk[nwrk - 2]);
    test = avr_plaquette_wrk();
    test -= plaq[0];

    dop = sqrt(_complex_prod_re(test, test));
    lprintf("MAIN", 0, "Checking that wrkfield %d memory is still avaialable.\n ", nwrk - 2);
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
    if (dop > 10.e-14) { return_value++; }
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Requesting and initializing 2 additional workspace gauge fields\n");
    idx_wrk[nwrk - 2] = reserve_wrk_space_with_pointers(wrk + nwrk - 2, iup_wrk + nwrk - 2, idn_wrk + nwrk - 2);
    random_u(wrk[nwrk - 2]);
    idx_wrk[nwrk] = reserve_wrk_space_with_pointers(wrk + nwrk, iup_wrk + nwrk, idn_wrk + nwrk);
    random_u(wrk[nwrk]);
    lprintf("MAIN", 0, "done.\n\n");

    test = avr_plaquette_wrk();
    test -= plaq[0];

    dop = sqrt(_complex_prod_re(test, test));
    lprintf("MAIN", 0, "Checking that wrkfield %d has been updated.\n ", nwrk - 2);
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
    lprintf("MAIN", 0, "(should be greater 1*10^(-6) or so)\n");
    if (dop < 10.e-6) { return_value++; }
    lprintf("MAIN", 0, "done.\n\n");

    int id;

    lprintf("MAIN", 0, "Checking that the default geometry pointer of the wrkfield %d are identical to the default ones.\n",
            nwrk - 2);
    for (id = 0; id < 4 * glattice.gsize_gauge; id++) {
        if (iup_wrk[nwrk - 2][id] != iup[id] || idn_wrk[nwrk - 2][id] != idn[id]) {
            lprintf("MAIN", 0, "Different pointers in %d\n", id);
            return_value++;
        }
    }
    lprintf("MAIN", 0, "done.\n\n");

    finalize_process();
    return return_value;
}
