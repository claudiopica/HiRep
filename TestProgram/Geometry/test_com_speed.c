/*******************************************************************************
*
* Test the comuunication speed
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
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
#include "spectrum.h"
#include "setup.h"
#include "random.h"
#define repts 500

int main(int argc, char *argv[])
{
    int i, j, k;

    /* setup process communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    init_BCs(NULL);

    random_u(u_gauge);

    struct timeval start, end, etime; /* //for trajectory timing */
    lprintf("MAIN", 0, "Starting the measure of time for the complete sendrecv (The measure is optimal when the timing is null)\n");

    double tmp = 0;
    double elapsed;
    long int flopcounter = 10;

    //warmup
    for (k = 0; k < 20; k++)
    {
        start_gf_sendrecv(u_gauge);

        complete_gf_sendrecv(u_gauge);
    }

    for (k = 0; k < 20; k++)
    {
        elapsed = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            start_gf_sendrecv(u_gauge);
            for (j = 0; j < flopcounter; j++)
                tmp = tmp + 1.0;

            gettimeofday(&start, 0);
            complete_gf_sendrecv(u_gauge);
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            elapsed += etime.tv_sec * 1000 + etime.tv_usec * 0.001;
        }
        elapsed /= repts;
        lprintf("MAIN", 0, "Complete send_receive for a bulk of %ld flop ops: time in [%lf millisec]\n", flopcounter, elapsed);
        flopcounter *= 2;
    }

    /* close communications */
    finalize_process();

    return 0;
}
