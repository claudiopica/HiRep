/******************************************************************************
* File: check_smeared_stoch_sources.c
*
* Checks of smeared stochastic sources
*
* Author: Pietro Butti
* Started: Sept. 2025
******************************************************************************/

#include "libhr.h"

/* Mesons parameters */
typedef struct input_mesons {
    char mstring[256];

    /* for the reading function */
    input_record_t read[5];
    double csw;
    double alpha; //PIETRO//
    int nhits_2pt;

} input_mesons;

#define init_input_mesons(varname)                                                                                \
    {                                                                                                             \
        .read = {                                                                                                 \
            { "fermion mass", "mes:mass = %s", STRING_T, (varname).mstring },                                     \
            { "csw coefficient", "mes:csw = %lf", DOUBLE_T, &(varname).csw },                                     \
            { "smearing alpha", "mes:alpha = %lf", DOUBLE_T, &(varname).alpha },                                  \
            { "number of noisy sources per cnfg for 2pt fn", "mes:nhits_2pt = %d", INT_T, &(varname).nhits_2pt }, \
            { NULL, NULL, INT_T, NULL }                                                                           \
        }                                                                                                         \
    }

input_glb glb_ip = init_input_glb(glb_ip);
input_mesons mes_ip = init_input_mesons(mes_ip);

void create_alltimes_pt_source_rand(spinor_field *source, int *locs) {
    int t, beta, col, lidx, ix;
    double ran;

    int Xsrc, Ysrc, Zsrc;

    for (beta = 0; beta < 4; ++beta) {
#ifdef WITH_GPU
        zero_spinor_field_cpu(&source[beta]);
#endif
        zero_spinor_field(&source[beta]);
    }

    if (PID == 0) {
        // Generate random spatial location
        for (t = 0; t < GLB_T; t++) {
            locs[4 * t] = t;
            ranlxd(&ran, 1);
            locs[1 + 4 * t] = (int)(ran * GLB_X);
            ranlxd(&ran, 1);
            locs[2 + 4 * t] = (int)(ran * GLB_Y);
            ranlxd(&ran, 1);
            locs[3 + 4 * t] = (int)(ran * GLB_Z);
        }
    } // end rank 0 global task
#ifdef WITH_MPI
    // Broadcast info to every processes
    MPI_Bcast(locs, 4 * GLB_T, MPI_INT, 0, GLB_COMM);
#endif

    for (t = 0; t < T; t++) {
        Xsrc = locs[1 + 4 * (zerocoord[0] + t)];
        Ysrc = locs[2 + 4 * (zerocoord[0] + t)];
        Zsrc = locs[3 + 4 * (zerocoord[0] + t)];

        if ((zerocoord[1] <= Xsrc && Xsrc < zerocoord[1] + X) && (zerocoord[2] <= Ysrc && Ysrc < zerocoord[2] + Y) &&
            (zerocoord[3] <= Zsrc && Zsrc < zerocoord[3] + Z)) {
            // Select point
            ix = ipt(t, Xsrc - zerocoord[1], Ysrc - zerocoord[2], Zsrc - zerocoord[3]);

            // Fill source
            for (col = 0; col < NF; col++) {
                for (beta = 0; beta < 4; beta++) {
                    lidx = beta + col * 4;
                    _FIELD_AT(&source[lidx], ix)->c[beta].c[col] = 1.;
                }
            }
        }
    }

    // -----------

    for (beta = 0; beta < 4; ++beta) {
#ifdef WITH_GPU
        copy_to_gpu(source + beta);
#endif
        start_sendrecv_spinor_field(source + beta);
        complete_sendrecv_spinor_field(source + beta);
    }
}

void create_alltimes_pt_source_loc(spinor_field *source, int Xsrc, int Ysrc, int Zsrc) {
    int t, beta, col, lidx, ix;

    for (beta = 0; beta < 4; ++beta) {
#ifdef WITH_GPU
        zero_spinor_field_cpu(&source[beta]);
#endif
        zero_spinor_field(&source[beta]);
    }

    // Check whether (t,xs,ys,zs) is in this process
    if (zerocoord[1] <= Xsrc && Xsrc < zerocoord[1] + X && zerocoord[2] <= Ysrc && Ysrc < zerocoord[2] + Y &&
        zerocoord[3] <= Zsrc && Zsrc < zerocoord[3] + Z) {
        for (t = 0; t < T; t++) {
            // Select point
            ix = ipt(t, Xsrc - zerocoord[1], Ysrc - zerocoord[2], Zsrc - zerocoord[3]);

            // Fill source
            for (col = 0; col < NF; col++) {
                for (beta = 0; beta < 4; beta++) {
                    lidx = beta + col * 4;
                    _FIELD_AT(&source[lidx], ix)->c[beta].c[col] = 1.;
                }
            }
        }
    }

    for (beta = 0; beta < 4; ++beta) {
#ifdef WITH_GPU
        copy_to_gpu(source + beta);
#endif
        start_sendrecv_spinor_field(source + beta);
        complete_sendrecv_spinor_field(source + beta);
    }
}

double smear_factor(int vecp[3], double alpha) {
    double cosp = cos(2. * M_PI * (double)vecp[0] / (double)GLB_X) + cos(2. * M_PI * (double)vecp[1] / (double)GLB_Y) +
                  cos(2. * M_PI * (double)vecp[2] / (double)GLB_Z);
    double z = 1. / (1. + 6. * alpha) * (1. + 2. * alpha * cosp);

    return z;
}

hr_complex exp_i_p_x(int p0, int p1, int p2, int x0, int x1, int x2) {
    return (cos(2. * M_PI / (double)GLB_X * (double)(p0 * x0)) + I * sin(2. * M_PI / (double)GLB_X * (double)(p0 * x0))) *
           (cos(2. * M_PI / (double)GLB_Y * (double)(p1 * x1)) + I * sin(2. * M_PI / (double)GLB_Y * (double)(p1 * x1))) *
           (cos(2. * M_PI / (double)GLB_Z * (double)(p2 * x2)) + I * sin(2. * M_PI / (double)GLB_Z * (double)(p2 * x2)));
}

int check_alltime_rand(spinor_field *source, hr_complex *acc, int *vecp, int *locs, double alpha) {
    int t, ix, iy, iz, ii, col, beta;

    hr_complex pw, tmp, accumulator;
    // hr_complex f_p = smear_factor(vecp, alpha);
    spinor_field *lsource;

    for (t = 0; t < T; t++) {
        // Fourier transfrom + spin/color sum in a block
        accumulator = 0. + I * 0.;
        for (ix = 0; ix < X; ix++) {
            for (iy = 0; iy < Y; iy++) {
                for (iz = 0; iz < Z; iz++) {
                    // Select point
                    ii = ipt(t, ix, iy, iz);

                    // Compute plane waves
                    pw = exp_i_p_x(-vecp[0], -vecp[1], -vecp[2], ix + zerocoord[1] - locs[1 + 4 * (zerocoord[0] + t)],
                                   iy + zerocoord[2] - locs[2 + 4 * (zerocoord[0] + t)],
                                   iz + zerocoord[3] - locs[3 + 4 * (zerocoord[0] + t)]);

                    // Cycle over spin/color components
                    for (int ls = 0; ls < NF * 4; ls++) {
                        lsource = source + ls;

                        for (col = 0; col < NF; col++) {
                            for (beta = 0; beta < 4; beta++) {
                                // int lidx = beta + col * 4;
                                // Accumulate
                                tmp = pw * _FIELD_AT(lsource, ii)->c[beta].c[col];
                                // if (creal(tmp) * creal(tmp) + cimag(tmp) * cimag(tmp) > 1.e-12) {
                                //     lprintf("not zero", 0, "point %d %d %d %d %d %d value %.10e %.10e\n", ix, iy, iz, t, col,
                                //             beta, creal(tmp), cimag(tmp));
                                // }

                                accumulator += tmp;
                            }
                        }
                    }
                }
            }
        }
        // Accumulate
        acc[t + zerocoord[0]] = accumulator;
    }

    global_sum((double *)acc, 2 * GLB_T);
    return 0;
}

int main(int argc, char *argv[]) {
    int n_sources;
    int beta, col, lidx;

    hr_complex *acc;

    logger_map("DEBUG", "debug");
    logger_setlevel(0, 200);

    /* setup process id and communications */
    setup_process(&argc, &argv);
    // read_input(pars_ip.read, get_input_filename());
    read_input(mes_ip.read, get_input_filename());

    setup_gauge_fields();
    unit_u(u_gauge);
    represent_gauge_field();

    lprintf("CORR", 0, "smearing parameter     : alpha = %e \n", mes_ip.alpha);

    n_sources = NF * 4;
    // spinor_field *source = alloc_spinor_field(n_sources, &glattice);
    // spinor_field *source_pt = alloc_spinor_field(n_sources, &glattice);

    // // Create one source per spin and color index
    // create_alltimes_pt_source_loc(source_pt, xs, ys, zs);

    int *src_loc = malloc(4 * GLB_T * sizeof(int));
    spinor_field *src = alloc_spinor_field(n_sources, &glattice);
    spinor_field *src_pt = alloc_spinor_field(n_sources, &glattice);

    create_alltimes_pt_source_rand(src_pt, src_loc);
    for (int t = 0; t < GLB_T; t++) {
        lprintf("MAIN", 0, "t=%i [%i, %i, %i]\n", src_loc[4 * t], src_loc[1 + 4 * t], src_loc[2 + 4 * t], src_loc[3 + 4 * t]);
    }

    // Smear the source
    for (col = 0; col < NF; col++) {
        for (beta = 0; beta < 4; beta++) {
            lidx = beta + col * 4;
            // Fphi_cpu_(&source[lidx], &source_pt[lidx], mes_ip.alpha);
            gaussian_smearing(&src[lidx], &src_pt[lidx], u_gauge_f, mes_ip.alpha);
        }
    }

    // Perform this to a representative of momenta -----------------------
    int moms[3 * 4] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1 };
    hr_complex f_p;

    for (int imom = 0; imom < 4; imom++) {
        int vecp[3] = { moms[3 * imom], moms[1 + 3 * imom], moms[2 + 3 * imom] };
        f_p = smear_factor(vecp, mes_ip.alpha) * (double)n_sources;

        lprintf("TEST", 0, "========= p = [%i, %i, %i], f(p) = %e + i %e ========= \n", vecp[0], vecp[1], vecp[2], creal(f_p),
                cimag(f_p));

        // Run test
        acc = (hr_complex *)malloc(GLB_T * sizeof(hr_complex));
        for (int t=0; t<GLB_T; t++) {
            acc[t] = 0.;
        }

        // Fourier transform in each block
        check_alltime_rand(src, acc, vecp, src_loc, mes_ip.alpha);

        for (int t=0; t<GLB_T; t++) {
            lprintf("TEST", 0, " source @ [t=%i, random]  res = %e %e\n", t, creal((acc[t]-f_p)/f_p), cimag((acc[t]-f_p)/f_p));
        }


    }

    // -------------------------------------------------------------------


    free_spinor_field(src);
    free_spinor_field(src_pt);

    finalize_process();

    return 0;
}