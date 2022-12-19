/*******************************************************************************
 * NOCOMPILE= !WITH_MPI
 * Test the communication speed
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
#include "representation.h"
#define repts 100

int flopmin;
int flopmax;
int intervals = 10;
int df;

#if defined(REPR_ADJOINT)
int flopsite = 8 * NF * (7 + 8 * NF);
#else
int flopsite = 8 * NF * (7 + 16 * NF);
#endif

int singlestep;

static void random_g(suNg_field *g)
{
    _MASTER_FOR(&glattice, ix)
    {
        random_suNg(_FIELD_AT(g, ix));
    }
}

static void transform(suNg_field *gtransf, suNg_field *gfield)
{
    _MASTER_FOR(&glattice, ix)
    {
        for (int mu = 0; mu < 4; mu++)
        {
            int iy = iup(ix, mu);
            suNg *u = _4FIELD_AT(gfield, ix, mu);
            suNg v;
            _suNg_times_suNg_dagger(v, *u, *_FIELD_AT(gtransf, iy));
            _suNg_times_suNg(*u, *_FIELD_AT(gtransf, ix), v);
        }
    }

    start_gf_sendrecv(gfield);
    complete_gf_sendrecv(gfield);
}

static void transform_s(suNg_field *gfield, spinor_field *in)
{
    suNf_vector tmp1, tmp2, tmp3, tmp4;
    suNf gfx;
    _MASTER_FOR(&glattice, ix)
    {
        suNf_spinor *s = _FIELD_AT(in, ix);

        _group_represent2(&gfx, _FIELD_AT(gfield, ix));

        _suNf_multiply(tmp1, gfx, s->c[0]);
        _suNf_multiply(tmp2, gfx, s->c[1]);
        _suNf_multiply(tmp3, gfx, s->c[2]);
        _suNf_multiply(tmp4, gfx, s->c[3]);

        s->c[0] = tmp1;
        s->c[1] = tmp2;
        s->c[2] = tmp3;
        s->c[3] = tmp4;
    }
}

int main(int argc, char *argv[])
{
    int i, j, k, l;
    suNg_field *g;
    spinor_field *s1;
    double norm;

    /* setup process communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    init_BCs(NULL);

    random_u(u_gauge);

    g = alloc_gtransf(&glattice);
    random_g(g);
    start_gt_sendrecv(g);
    complete_gt_sendrecv(g);

    s1 = alloc_spinor_field_f(1, &glattice);

    gaussian_spinor_field(s1);

    struct timeval start, end, etime; /* //for trajectory timing */

    double tmp = 0;
    double opt_trick = 0.;
    double elapsed, elapsed_mean, elapsed_var;
    double msf, mgf, mbulk;

    lprintf("MAIN", 0, "Reference: the number of flop per site in the actual setup for the application of the plain Dirac operator to a full spinor is %ld \n", flopsite);
    lprintf("MAIN", 0, "The number of flop in your bulk in the actual setup for the application of the plain Dirac operator is %ld \n", flopsite * (glattice.master_end[0] - glattice.master_start[0] + 1));

    singlestep = flopsite * (glattice.master_end[0] - glattice.master_start[0] + 1) / 100;

    //warmup
    for (k = 0; k < 20; k++)
    {
        start_gf_sendrecv(u_gauge);

        complete_gf_sendrecv(u_gauge);
    }

    MPI_Barrier(GLB_COMM);
    elapsed = 0.;
    elapsed_mean = 0.;
    elapsed_var = 0.;

    for (k = 0; k < repts; k++)
    {
        gettimeofday(&start, 0);

        start_gf_sendrecv(u_gauge);

        complete_gf_sendrecv(u_gauge);

        gettimeofday(&end, 0);
        transform(g, u_gauge);

        timeval_subtract(&etime, &end, &start);
        elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
        elapsed_mean += elapsed;
        elapsed_var += elapsed * elapsed;
    }

    elapsed_mean /= repts;
    elapsed_var /= repts;
    elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));
    lprintf("MAIN", 0, "Full Send Recv Reference +-+-+-+-+-+-+\n");

    lprintf("MAIN", 0, "start_gf_send +complete_gf_recv timing in [%lf millisec] +- [%lf millisec]\n", elapsed_mean, elapsed_var);
    MPI_Barrier(GLB_COMM);
    mgf = elapsed_mean;

    elapsed = 0.;
    elapsed_mean = 0.;
    elapsed_var = 0.;

    for (k = 0; k < repts; k++)
    {
        gettimeofday(&start, 0);

        start_sf_sendrecv(s1);

        complete_sf_sendrecv(s1);

        gettimeofday(&end, 0);
        transform_s(g, s1);

        timeval_subtract(&etime, &end, &start);
        elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
        elapsed_mean += elapsed;
        elapsed_var += elapsed * elapsed;
    }

    elapsed_mean /= repts;
    elapsed_var /= repts;
    elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

    lprintf("MAIN", 0, "start_sf_send +complete_gf_recv timing in [%lf millisec] +- [%lf millisec]\n", elapsed_mean, elapsed_var);
    MPI_Barrier(GLB_COMM);
    msf = elapsed_mean;

    elapsed = 0.;
    elapsed_mean = 0.;
    elapsed_var = 0.;
    for (i = 0; i < repts; i++)
    {
        tmp = 0;
        gettimeofday(&start, 0);

        for (l = 0; l < singlestep / 2; l++)
            tmp = tmp + l * 0.00021234512;

        gettimeofday(&end, 0);
        timeval_subtract(&etime, &end, &start);
        elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
        elapsed_mean += elapsed;
        elapsed_var += elapsed * elapsed;
        opt_trick += tmp;
    }
    norm = 1.3432 * opt_trick;
    opt_trick /= norm;

    elapsed_mean /= repts;
    mbulk = elapsed_mean;

    flopmin = (int)(0.1 * msf / mbulk);

    if (flopmin > 20)
        flopmin = 20;

    flopmax = (int)(4.0 * mgf / mbulk);
    df = (flopmax - flopmin) / intervals;

    int *steps = malloc((intervals + 1) * sizeof(int));
    int id = 0;
    for (k = 0; k < intervals; k++)
    {
        steps[id] = flopmin + df * id;
        if (flopmin + df * (id + 1) > 100 && flopmin + df * (id) < 100)
        {
            id++;
            steps[id] = 100;
        }
        id++;
    }
    lprintf("MAIN", 0, "Bulk timing +-+-+-+-+-+-+\n");

    for (k = 0; k < intervals + 1; k++)
    {
        elapsed = 0.;
        elapsed_mean = 0.;
        elapsed_var = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            gettimeofday(&start, 0);
            for (j = 0; j < steps[k]; j++)
                for (l = 0; l < singlestep / 2; l++)
                    tmp = tmp + l * 0.00021234512;

            gettimeofday(&end, 0);
            MPI_Barrier(GLB_COMM);
            timeval_subtract(&etime, &end, &start);
            elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
            elapsed_mean += elapsed;
            elapsed_var += elapsed * elapsed;
            opt_trick += tmp;
        }
        elapsed_mean /= repts;
        elapsed_var /= repts;
        elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

        lprintf("MAIN", 0, "%ld flop ops: time in [%lf millisec] +- [%lf millisec]\n", steps[k] * singlestep, elapsed_mean, elapsed_var);
    }
    MPI_Barrier(GLB_COMM);
    opt_trick /= norm;
    lprintf("MAIN", 0, "MPI Barrier timing +-+-+-+-+-+-+\n");

    for (k = 0; k < intervals + 1; k++)
    {
        elapsed = 0.;
        elapsed_mean = 0.;
        elapsed_var = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            for (j = 0; j < steps[k]; j++)
                for (l = 0; l < singlestep / 2; l++)
                    tmp = tmp + l * 0.00021234512;

            gettimeofday(&start, 0);
            MPI_Barrier(GLB_COMM);
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
            elapsed_mean += elapsed;
            elapsed_var += elapsed * elapsed;
            opt_trick += tmp;
        }
        elapsed_mean /= repts;
        elapsed_var /= repts;
        elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

        lprintf("MAIN", 0, "MPI_Barrier after %ld flop ops: time in [%lf millisec] +- [%lf millisec]\n", steps[k] * singlestep, elapsed_mean, elapsed_var);
    }
    MPI_Barrier(GLB_COMM);

    opt_trick /= norm;
    lprintf("MAIN", 0, "SF timing complete_sf_send_recv +-+-+-+-+-+-+\n");

    for (k = 0; k < intervals + 1; k++)
    {
        elapsed = 0.;
        elapsed_mean = 0.;
        elapsed_var = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            start_sf_sendrecv(s1);
            for (j = 0; j < steps[k]; j++)
                for (l = 0; l < singlestep / 2; l++)
                    tmp = tmp + l * 0.0000021234512;

            gettimeofday(&start, 0);
            complete_sf_sendrecv(s1);
            gettimeofday(&end, 0);
            transform_s(g, s1);

            timeval_subtract(&etime, &end, &start);
            elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
            elapsed_mean += elapsed;
            elapsed_var += elapsed * elapsed;
            opt_trick += tmp;
        }
        elapsed_mean /= repts;
        elapsed_var /= repts;
        elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

        lprintf("MAIN", 0, "complete_sf_send_recv for a bulk of %ld flop ops: time in [%lf millisec] +- [%lf millisec]\n", steps[k] * singlestep, elapsed_mean, elapsed_var);
    }

    opt_trick /= norm;

    MPI_Barrier(GLB_COMM);
    lprintf("MAIN", 0, "SF timing start_sf_send + bulk  + complete_sf_recv +-+-+-+-+-+-+\n");

    for (k = 0; k < intervals + 1; k++)
    {
        elapsed = 0.;
        elapsed_mean = 0.;
        elapsed_var = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            gettimeofday(&start, 0);
            start_sf_sendrecv(s1);
            for (j = 0; j < steps[k]; j++)
                for (l = 0; l < singlestep / 2; l++)
                    tmp = tmp + l * 0.00021234512;

            complete_sf_sendrecv(s1);
            gettimeofday(&end, 0);
            transform_s(g, s1);

            timeval_subtract(&etime, &end, &start);
            elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
            elapsed_mean += elapsed;
            elapsed_var += elapsed * elapsed;
            opt_trick += tmp;
        }
        elapsed_mean /= repts;
        elapsed_var /= repts;
        elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

        lprintf("MAIN", 0, "start_sf_send + bulk of %ld flop ops + complete_sf_recv : time in [%lf millisec] +- [%lf millisec]\n", steps[k] * singlestep, elapsed_mean, elapsed_var);
    }

    opt_trick /= norm;

    MPI_Barrier(GLB_COMM);
    lprintf("MAIN", 0, "Gauge timing complete_gf_send_recv +-+-+-+-+-+-+\n");

    for (k = 0; k < intervals + 1; k++)
    {
        elapsed = 0.;
        elapsed_mean = 0.;
        elapsed_var = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            start_gf_sendrecv(u_gauge);
            for (j = 0; j < steps[k]; j++)
                for (l = 0; l < singlestep / 2; l++)
                    tmp = tmp + l * 0.00021234512;

            gettimeofday(&start, 0);
            complete_gf_sendrecv(u_gauge);
            gettimeofday(&end, 0);
            transform(g, u_gauge);

            timeval_subtract(&etime, &end, &start);
            elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
            elapsed_mean += elapsed;
            elapsed_var += elapsed * elapsed;
            opt_trick += tmp;
        }
        elapsed_mean /= repts;
        elapsed_var /= repts;
        elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

        lprintf("MAIN", 0, "complete_gf_send_recv for a bulk of %ld flop ops: time in [%lf millisec] +- [%lf millisec]\n", steps[k] * singlestep, elapsed_mean, elapsed_var);
    }
    opt_trick /= norm;

    MPI_Barrier(GLB_COMM);
    lprintf("MAIN", 0, "Gauge timing start_gf_send + bulk  + complete_gf_recv +-+-+-+-+-+-+\n");

    for (k = 0; k < intervals + 1; k++)
    {
        elapsed = 0.;
        elapsed_mean = 0.;
        elapsed_var = 0.;
        for (i = 0; i < repts; i++)
        {
            tmp = 0;
            gettimeofday(&start, 0);
            start_gf_sendrecv(u_gauge);
            for (j = 0; j < steps[k]; j++)
                for (l = 0; l < singlestep / 2; l++)
                    tmp = tmp + l * 0.00021234512;

            complete_gf_sendrecv(u_gauge);
            gettimeofday(&end, 0);
            transform(g, u_gauge);

            timeval_subtract(&etime, &end, &start);
            elapsed = etime.tv_sec * 1000 + etime.tv_usec * 0.001;
            elapsed_mean += elapsed;
            elapsed_var += elapsed * elapsed;
            opt_trick += tmp;
        }
        elapsed_mean /= repts;
        elapsed_var /= repts;
        elapsed_var = sqrt((elapsed_var - elapsed_mean * elapsed_mean) / (repts - 1));

        lprintf("MAIN", 0, "start_gf_send + bulk of %ld flop ops + complete_gf_recv : time in [%lf millisec] +- [%lf millisec]\n", steps[k] * singlestep, elapsed_mean, elapsed_var);
    }

    lprintf("MAIN", 0, "check number %1.10e (must be differnt from zero)\n ", opt_trick);

    /* close communications */
    finalize_process();

    return 0;
}
