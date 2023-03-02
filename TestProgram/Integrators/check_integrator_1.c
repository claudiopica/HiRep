/****************************************************************************

*
*******************************************************************************/
#include "check_integrator_1.h"

#define SCALING_RANGE 20

/* flow control variable */
hmc_flow flow = init_hmc_flow(flow);

int main(int argc, char *argv[]) {
    int return_value = 0;
    double res, elapsed_lf, elapsed_o2, elapsed_o4;
    /* setup process communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    /* Init Monte Carlo */

    init_mc_ghmc(&flow, get_input_filename());

    integrator_par *pint = flow.hmc_v->hmc_p.integrator;
    do {
        lprintf("Main", 0, "number of steps per monomial level %ld\n", (long int)pint->nsteps);
        pint = pint->next;
    } while (pint != NULL);

    lprintf("MAIN", 0, "Initial plaquette: %1.16e\n", avr_plaquette());

    double rr, rr1;
    struct timeval start, end, etime; /* //for trajectory timing */
    gettimeofday(&start, 0);
    rr = integrate_ghmc(0, &(flow.hmc_v->hmc_p));

    rr1 = integrate_ghmc(1, &(flow.hmc_v->hmc_p));
    gettimeofday(&end, 0);

    timeval_subtract(&etime, &end, &start);
    res = fabs(rr1 - rr) / fabs(rr1);

    lprintf(
        "MAIN", 0,
        "Checking that the user defined action doesn't change on reusing the integrator twice from the same starting fields %.2e\n",
        res);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n\n");

    if (res > 1.0e-12) { return_value++; }

    for (int int_type = 0; int_type < 3; int_type++) {
        set_integrator_type(&(flow.hmc_v->hmc_p), int_type);

        gettimeofday(&start, 0);
        rr = integrate_ghmc(1, &(flow.hmc_v->hmc_p));

        rr1 = integrate_ghmc(1, &(flow.hmc_v->hmc_p));
        gettimeofday(&end, 0);

        timeval_subtract(&etime, &end, &start);
        res = fabs(rr1 - rr) / fabs(rr1);

        lprintf("MAIN", 0,
                "Checking that the action doesn't change on reusing the integrator twice from the same starting fields %.2e\n",
                res);
        lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n\n");

        if (res > 1.0e-12) { return_value++; }
    }
    if (1 == 1) {
        double rt0[SCALING_RANGE], rt1[SCALING_RANGE], rt2[SCALING_RANGE];

        for (int i = 0; i < SCALING_RANGE; i++) {
            set_integrator_type(&(flow.hmc_v->hmc_p), LEAPFROG);
            set_first_integrator_nsteps(&(flow.hmc_v->hmc_p), i * 10 + 20);
            gettimeofday(&start, 0);
            rt0[i] = integrate_ghmc(1, &(flow.hmc_v->hmc_p));
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            elapsed_lf = etime.tv_sec * 1000. + etime.tv_usec * 0.001;

            set_integrator_type(&(flow.hmc_v->hmc_p), O2MN);
            gettimeofday(&start, 0);
            rt1[i] = integrate_ghmc(1, &(flow.hmc_v->hmc_p));
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            elapsed_o2 = etime.tv_sec * 1000. + etime.tv_usec * 0.001;

            set_integrator_type(&(flow.hmc_v->hmc_p), O4MN);
            gettimeofday(&start, 0);
            rt2[i] = integrate_ghmc(1, &(flow.hmc_v->hmc_p));
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            elapsed_o4 = etime.tv_sec * 1000. + etime.tv_usec * 0.001;

            lprintf(
                "MAIN", 0,
                "Delta H for first level monomial nsteps: %d dt: %lf LF: %1.16e %1.16e msec  O2: %1.16e %1.16e msec  O4: %1.16e %1.16e msec\n\n",
                i * 10 + 20, 1.0 / ((double)(i * 10 + 20)), fabs(rt0[i]), elapsed_lf, fabs(rt1[i]), elapsed_o2, fabs(rt2[i]),
                elapsed_o4);
        }

        double scaling_deviation1;
        double scaling_deviation2;
        double dt, ldt = 1.0 / ((SCALING_RANGE - 1) * 10 + 20);
        for (int i = 0; i < SCALING_RANGE - 1; i++) {
            dt = 1.0 / (i * 10 + 20);
            scaling_deviation1 = fabs(1.0 - rt0[i] * rt1[SCALING_RANGE - 1] / rt1[i] / rt0[SCALING_RANGE - 1]);
            scaling_deviation2 =
                fabs(1.0 - dt * dt * rt1[i] * rt2[SCALING_RANGE - 1] / rt2[i] / rt1[SCALING_RANGE - 1] / ldt / ldt);
            lprintf(
                "MAIN", 0,
                "Scaling deviation for LF/O2 (normalized at the smallest dt) at dt=%lf %1.16e \n(should be around 1*10^(-2) or so)\n\n",
                dt, scaling_deviation1);
            lprintf(
                "MAIN", 0,
                "Scaling deviation for dt^2*O2/O4 (normalized at the smallest dt) at dt=%lf %1.16e \n(should be around 1*10^(-2) or so)\n\n",
                dt, scaling_deviation2);

            if (scaling_deviation1 > 0.5 || scaling_deviation2 > 0.5) { return_value++; }
        }
    }

    /* close communications */
    finalize_process();

    return return_value;
}
