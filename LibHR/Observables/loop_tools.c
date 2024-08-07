/******************************************************************************
 *
 * Methods to compute disconnected loops
 * Copyright (c) 2014, R. Arthur, V. Drach, A. Hietanen
 * All rights reserved.
 * 
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 *
 *******************************************************************************/

#include "observables.h"
#include "libhr_core.h"
#include "io.h"
#include "utils.h"
#include "memory.h"
#include "Update/representation.h"
#include "Update/avr_plaquette.h"
#include "inverters.h"

#if (!defined(BC_T_SF_ROTATED) && !defined(BC_T_SF) && !defined(FERMION_THETA))

#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

void measure_bilinear_loops_4spinorfield(spinor_field *prop, spinor_field *source, int src_id, int tau, int col, int eo,
                                         storage_switch swc, data_storage_array **ret) {
    hr_complex **corr;
    double *corr_re[16];
    double *corr_im[16];
    suNf_spin_matrix sma, smb, sm1;
    int NGamma = 16;
    int ix, t, x, y, z, tc, beta;
    int iGamma;
    int tau_min = 0;
    int tau_max = GLB_T;
    hr_complex tr;

    struct timeval start, end, etime;

    gettimeofday(&start, 0);

    for (int i = 0; i < NGamma; i++) {
        corr_re[i] = (double *)calloc(GLB_T, sizeof(double));
        corr_im[i] = (double *)calloc(GLB_T, sizeof(double));
    }

    corr = (hr_complex **)malloc(sizeof(hr_complex *) * NGamma);
    for (int i = 0; i < NGamma; i++) {
        corr[i] = (hr_complex *)calloc(GLB_T, sizeof(hr_complex));
    }

    if (tau != -1) {
        tau_min = tau;
        tau_max = tau + 1;
    }

    for (t = 0; t < T; t++) {
        tc = (zerocoord[0] + t + GLB_T) % GLB_T;
        /* loop on spatial volume */
        for (x = 0; x < X; x++) {
            for (y = 0; y < Y; y++) {
                for (z = 0; z < Z; z++) {
                    ix = ipt(t, x, y, z);

                    for (beta = 0; beta < 4; beta++) {
                        _spinmatrix_assign_row(sma, *_FIELD_AT(&source[beta], ix), beta);
                        _spinmatrix_assign_row(smb, *_FIELD_AT(&prop[beta], ix), beta);
                    }

                    _g5_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[0][tc], tr);

                    _g1_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[1][tc], tr);
                    _g2_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[2][tc], tr);
                    _g3_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[3][tc], tr);

                    _g5g0_spinmatrix(sm1, smb); // minus sign is missing
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[4][tc], tr);

                    _g0g1_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[5][tc], tr);
                    _g0g2_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[6][tc], tr);
                    _g0g3_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[7][tc], tr);
                    _spinmatrix_mul_trace(tr, sma, smb);
                    _complex_add_assign(corr[8][tc], tr);

                    _g5g1_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[9][tc], tr);
                    _g5g2_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[10][tc], tr);
                    _g5g3_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[11][tc], tr);

                    _g0_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[12][tc], tr);
                    _g5g0g1_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[13][tc], tr);
                    _g5g0g2_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[14][tc], tr);
                    _g5g0g3_spinmatrix(sm1, smb);
                    _spinmatrix_mul_trace(tr, sma, sm1);
                    _complex_add_assign(corr[15][tc], tr);
                }
            }
        }

    } /* end loop t */

    for (t = 0; t < T; t++) {
        tc = (zerocoord[0] + t + GLB_T) % GLB_T;
        for (iGamma = 0; iGamma < NGamma; iGamma++) {
            if (iGamma != 0 && iGamma != 1 && iGamma != 2 && iGamma != 3 && iGamma != 8 && iGamma != 12) {
                /* multiply by -i to match old convention [and to have hermitean operators ] */
                corr[iGamma][tc] = -I * corr[iGamma][tc];
            }
            corr_re[iGamma][tc] = creal(corr[iGamma][tc]);
            corr_im[iGamma][tc] = cimag(corr[iGamma][tc]);
        }
    }

    /* print into output file */
    for (iGamma = 0; iGamma < NGamma; iGamma++) {
        global_sum(corr_re[iGamma], GLB_T);
        global_sum(corr_im[iGamma], GLB_T);
        global_sum((double *)(corr[iGamma]), 2 * GLB_T);
    }

    if (src_id == 0 && tau == -1 && col == -1 && eo == -1) {
        lprintf("CORR", 0, " Output format [t iGamma iSource Re Im] \n"); // no dilution
    }
    if (src_id == 0 && tau == 0 && col == -1 && eo == -1) {
        lprintf("CORR", 0, "Output format  [t iGamma iSource Re Im] \n"); // time dilution
    }
    if (src_id == 0 && tau == 0 && col == 0 && eo == -1) {
        lprintf("CORR", 0, "Output format  [t iGamma iSource col Re Im] \n"); // time + col dilution
    }
    if (src_id == 0 && tau == 0 && col == 0 && eo == 0) {
        lprintf("CORR", 0, "Output format  [t iGamma iSource col eo Re Im] \n"); // time + col + eo dilution
    }
    if (src_id == 0 && tau == -1 && col == 0 && eo == -1) {
        lprintf("CORR", 0, "Output format  [t iGamma iSource col Re Im] \n"); //   col dilution
    }
    if (src_id == 0 && tau == -1 && col == 0 && eo == 0) {
        lprintf("CORR", 0, "Output format  [t iGamma iSource col eo Re Im] \n"); // col + eo dilution
    }

    if (col == -1 && eo == -1) {
        for (t = tau_min; t < tau_max; ++t) {
            for (iGamma = 0; iGamma < NGamma; iGamma++) {
                lprintf("CORR", 0, "%i %i %i %3.10e %3.10e \n", t, iGamma, src_id, corr_re[iGamma][t], corr_im[iGamma][t]);
                if (swc == STORE) {
                    int idx[4] = { src_id, iGamma, t, 0 };
                    *data_storage_element(*ret, 0, idx) = corr_re[iGamma][t];
                    idx[3] = 1;
                    *data_storage_element(*ret, 0, idx) = corr_im[iGamma][t];
                }
            }
        }
    }

    if (col != -1 && eo == -1) {
        for (t = tau_min; t < tau_max; ++t) {
            for (iGamma = 0; iGamma < NGamma; iGamma++) {
                lprintf("CORR", 0, "%i %i %i %i %3.10e %3.10e \n", t, iGamma, src_id, col, corr_re[iGamma][t],
                        corr_im[iGamma][t]);
                if (swc == STORE) {
                    int idx[5] = { src_id, col, iGamma, t, 0 };
                    *data_storage_element(*ret, 0, idx) = corr_re[iGamma][t];
                    idx[4] = 1;
                    *data_storage_element(*ret, 0, idx) = corr_im[iGamma][t];
                }
            }
        }
    }
    if (col == -1 && eo != -1) {
        for (t = tau_min; t < tau_max; ++t) {
            for (iGamma = 0; iGamma < NGamma; iGamma++) {
                lprintf("CORR", 0, "%i %i %i %i %3.10e %3.10e \n", t, iGamma, src_id, eo, corr_re[iGamma][t],
                        corr_im[iGamma][t]);
                if (swc == STORE) {
                    int idx[5] = { src_id, eo, iGamma, t, 0 };
                    *data_storage_element(*ret, 0, idx) = corr_re[iGamma][t];
                    idx[4] = 1;
                    *data_storage_element(*ret, 0, idx) = corr_im[iGamma][t];
                }
            }
        }
    }
    if (col != -1 && eo != -1) {
        for (t = tau_min; t < tau_max; ++t) {
            for (iGamma = 0; iGamma < NGamma; iGamma++) {
                lprintf("CORR", 0, "%i %i %i %i %i %3.10e %3.10e \n", t, iGamma, src_id, col, eo, corr_re[iGamma][t],
                        corr_im[iGamma][t]);
                if (swc == STORE) {
                    int idx[6] = { src_id, eo, col, iGamma, t, 0 };
                    *data_storage_element(*ret, 0, idx) = corr_re[iGamma][t];
                    idx[5] = 1;
                    *data_storage_element(*ret, 0, idx) = corr_im[iGamma][t];
                }
            }
        }
    }

    fflush(stdout);
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    if (tau != -1) {
        lprintf("TIMING", 0, "Contractions for source %d done and time %d [%ld sec %ld usec]\n", src_id, tau, etime.tv_sec,
                etime.tv_usec);
    }
    if (tau == -1) {
        lprintf("TIMING", 0, "Contractions for source %d done and all timeslices [%ld sec %ld usec]\n", src_id, etime.tv_sec,
                etime.tv_usec);
    }
}

void measure_loops(double *m, int nhits, int conf_num, double precision, int source_type, int n_mom, storage_switch swc,
                   data_storage_array **ret) {
    int k, l;
    int n_spinor;
    int eo, tau, col;
    struct timeval start, end, etime;

    if (source_type == 0) { lprintf("CORR", 0, "Pure volume source  will be used  \n"); }
    if (source_type == 1) { lprintf("CORR", 0, "Gauge fixed source  with time and spin dilution will be used \n"); }
    if (source_type == 2) { lprintf("CORR", 0, "Time and spin dilution  will be used \n"); }
    if (source_type == 3) { lprintf("CORR", 0, "Time, spin and color dilution  will be used \n"); }
    if (source_type == 4) { lprintf("CORR", 0, "Time, spin , color and eo dilution  will be used \n"); }
    if (source_type == 5) { lprintf("CORR", 0, "Spin , color and eo dilution  will be used \n"); }

    gettimeofday(&start, 0);
    init_propagator_eo(1, m, precision);

    spinor_field *source;
    spinor_field *prop;
    suNg_field *u_gauge_old = NULL;

    if (source_type == 0) {
        source = alloc_spinor_field(1, &glattice);
        prop = alloc_spinor_field(1, &glattice);
#ifdef WITH_GPU
        zero_spinor_field_cpu(prop);
#endif
        zero_spinor_field(prop);
    } else {
        source = alloc_spinor_field(4, &glattice);
        prop = alloc_spinor_field(4, &glattice);
        for (int i = 0; i < 4; i++) {
#ifdef WITH_GPU
            zero_spinor_field_cpu(prop + i);
#endif
            zero_spinor_field(prop + i);
        }
    }

    if (swc == STORE && *ret == NULL) {
        if (source_type == 0) {
            int idx[5] = { nhits, pow(n_mom, 3), 16, GLB_T, 2 };
            *ret = allocate_data_storage_array(1);
            allocate_data_storage_element(*ret, 0, 5, idx); // ( 1 ) * (nhits*nmom^3*ngamma*GLB_T * 2 reals )
        } else if (source_type == 2) {
            int idx[4] = { nhits, 16, GLB_T, 2 };
            *ret = allocate_data_storage_array(1);
            allocate_data_storage_element(*ret, 0, 4, idx); // ( 1 ) * (nhits*ngamma*GLB_T * 2 reals )
        } else if (source_type == 1 || source_type == 3) {
            int idx[5] = { nhits, NF, 16, GLB_T, 2 };
            *ret = allocate_data_storage_array(1);
            allocate_data_storage_element(*ret, 0, 5, idx); // ( 1 ) * (nhits*NF*ngamma*GLB_T * 2 reals )
        } else if (source_type == 4) {
            int idx[6] = { nhits, 2, NF, 16, GLB_T, 2 };
            *ret = allocate_data_storage_array(1);
            allocate_data_storage_element(*ret, 0, 6, idx); // ( 1 ) * (nhits*2(eo)*NF*ngamma*GLB_T * 2 reals )
        } else if (source_type == 5) {
            int idx[6] = { nhits, 2, NF, 16, GLB_T, 2 };
            *ret = allocate_data_storage_array(1);
            allocate_data_storage_element(*ret, 0, 6, idx); // ( 1 ) * (nhits*2(eo)*NF*ngamma*GLB_T * 2 reals )
        } else {
            error(1, 1, "measure_loops [loop_tools.c]", "Source_type not implemented");
        }
    }
    if (source_type == 1) {
        u_gauge_old = alloc_suNg_field(&glattice);
        copy_suNg_field(u_gauge_old, u_gauge);
#ifdef WITH_GPU
        zero_spinor_field_cpu(prop);
#endif
        zero_spinor_field(prop);
        //Fix the Gauge
        double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
                              1.8, //overrelax
                              10000, //maxit
                              1e-12, //tolerance
                              u_gauge //gauge
        );
        lprintf("GFWALL", 0, "Gauge fixed action  %1.6f\n", act);
        double p2 = calc_plaq(u_gauge);
        lprintf("TEST", 0, "fixed_gauge plaq %1.6f\n", p2);
        full_plaquette();
        represent_gauge_field();
    }

    for (k = 0; k < nhits; k++) {
        lprintf("MAIN", 0, "k = %d/%d source_type=%d\n", k + 1, nhits, source_type);
        if (source_type == 0) /* generation of a volume source with Z2xZ2 noise */
        {
            create_z2_volume_source(source);
            calc_propagator(prop, source, 1); // No dilution
#ifdef WITH_GPU
            copy_from_gpu(prop);
#endif
            lprintf("CORR", 0, "Start to perform the contractions ... \n");
            measure_bilinear_loops_spinorfield(prop, source, k, n_mom, swc, ret);
            lprintf("CORR", 0, "Contraction done\n");
        }

        if (source_type == 1) // experimental Gauge Fixed Wall sources
        {
            //error(1, 1, "measure_loops [loop_tools.c]", "Source_type ==1 (gauge fixed wall sources) is broken and untested.");

            for (tau = 0; tau < GLB_T; tau++) {
                for (l = 0; l < NF; l++) {
                    create_gauge_fixed_wall_source(source, tau, l);
                    calc_propagator(prop, source, 4); //4 for spin dilution
                    create_point_source(source, tau, l); //to get the contraction right
                    //measure_mesons(discon_correlators, prop, source, 1, 0);
#ifdef WITH_GPU
                    for (int beta = 0; beta < 4; beta++) {
                        copy_from_gpu(prop + beta);
                    }
#endif
                    measure_bilinear_loops_4spinorfield(prop, source, k, tau, l, -1, swc, ret);
                }
            }

            //This gets the norm of the 2pt wrong by a factor GLB_VOL3 but the norm of the disconnected right
            //print_mesons(discon_correlators, 1.0, conf_num, 1, m, GLB_T, 1, "DISCON_GFWALL");

        } /* gfwall */

        if (source_type == 2) {
            for (tau = 0; tau < GLB_T; ++tau) {
                create_diluted_source_equal_atau(source, tau);
                calc_propagator(prop, source, 4); //4 for spin dilution
#ifdef WITH_GPU
                for (int beta = 0; beta < 4; beta++) {
                    copy_from_gpu(prop + beta);
                }
#endif
                measure_bilinear_loops_4spinorfield(prop, source, k, tau, -1, -1, swc, ret);
            }

        } /* time + spin dilution  */

        if (source_type == 3) {
            for (tau = 0; tau < GLB_T; ++tau) {
                for (col = 0; col < NF; ++col) {
                    create_diluted_source_equal_atau_col(source, tau, col);
                    calc_propagator(prop, source, 4); //4 for spin dilution
#ifdef WITH_GPU
                    for (int beta = 0; beta < 4; beta++) {
                        copy_from_gpu(prop + beta);
                    }
#endif
                    measure_bilinear_loops_4spinorfield(prop, source, k, tau, col, -1, swc, ret);
                }
            }
        } /* time  + spin + color dilution  */

        if (source_type == 4) {
            n_spinor = 4;
            for (tau = 0; tau < GLB_T; ++tau) {
                for (col = 0; col < NF; ++col) {
                    for (eo = 0; eo < 2; ++eo) {
                        create_diluted_source_equal_atau_col(source, tau, col);
                        zero_even_or_odd_site_spinorfield(source, n_spinor, eo);
                        calc_propagator(prop, source, 4); //4 for spin dilution
#ifdef WITH_GPU
                        for (int beta = 0; beta < 4; beta++) {
                            copy_from_gpu(prop + beta);
                        }
#endif
                        measure_bilinear_loops_4spinorfield(prop, source, k, tau, col, eo, swc, ret);
                    }
                }
            }
        } /* time + spin + color +eo  dilution  */

        if (source_type == 5) {
            n_spinor = 4;
            for (col = 0; col < NF; ++col) {
                for (eo = 0; eo < 2; ++eo) {
                    create_noise_source_equal_col_dil(source, col);
                    zero_even_or_odd_site_spinorfield(source, n_spinor, eo); //set even or odd site to zero
                    calc_propagator(prop, source, 4); //4 for spin dilution
#ifdef WITH_GPU
                    for (int beta = 0; beta < 4; beta++) {
                        copy_from_gpu(prop + beta);
                    }
#endif
                    measure_bilinear_loops_4spinorfield(prop, source, k, -1, col, eo, swc, ret);
                }
            }
        } /* volume source + spin + color + eo  dilution  */
    }
    if (u_gauge_old != NULL) {
        copy_suNg_field(u_gauge, u_gauge_old);
        represent_gauge_field();
        free_suNg_field(u_gauge_old);
    }
    free_spinor_field(source);
    free_spinor_field(prop);
    free_propagator_eo();
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("TIMING", 0, "Sources generation, invert and contract for %i sources done [%ld sec %ld usec]\n", nhits,
            etime.tv_sec, etime.tv_usec);
}

void measure_bilinear_loops_spinorfield(spinor_field *prop, spinor_field *source, int src_id, int n_mom, storage_switch swc,
                                        data_storage_array **ret) {
    int px, py, pz, ip;
    int NGamma = 16;
    int n_mom_tot = n_mom * n_mom * n_mom;
    double pdotx;
    hr_complex phase;
    hr_complex ***corr;
    double *corr_re[n_mom_tot][16];
    double *corr_im[n_mom_tot][16];

    int pt[4];
    pt[0] = pt[1] = pt[2] = pt[3] = 0;

    int i, j, ix, t, x, y, z, tc;
    corr = (hr_complex ***)malloc(sizeof(hr_complex **) * n_mom_tot);

    for (i = 0; i < n_mom_tot; i++) {
        corr[i] = (hr_complex **)malloc(sizeof(hr_complex *) * NGamma);
    }
    for (i = 0; i < n_mom_tot; i++) {
        for (j = 0; j < NGamma; j++) {
            corr[i][j] = (hr_complex *)calloc(GLB_T, sizeof(hr_complex));
        }
    }

    int offset = 0;
    suNf_spinor tmp_spinor;
    hr_complex tmp;
    struct timeval start, end, etime;

    gettimeofday(&start, 0);

    for (j = 0; j < n_mom_tot; ++j) {
        for (i = 0; i < NGamma; ++i) {
            //corr[j][i]=(hr_complex*) malloc(sizeof(hr_complex)*size);
            corr_re[j][i] = (double *)calloc(GLB_T, sizeof(double));
            corr_im[j][i] = (double *)calloc(GLB_T, sizeof(double));
        }
    }

    for (t = 0; t < T; t++) {
        /* offset set to zero here */
        tc = (zerocoord[0] + t + GLB_T) % GLB_T + offset;
        /* loop on spatial volume */
        for (x = 0; x < X; x++) {
            for (y = 0; y < Y; y++) {
                for (z = 0; z < Z; z++) {
                    ix = ipt(t, x, y, z);
                    ip = 0;
                    for (px = 0; px < n_mom; ++px) {
                        for (py = 0; py < n_mom; ++py) {
                            for (pz = 0; pz < n_mom; ++pz) {
                                pdotx = 2.0 * PI *
                                        (((double)px) * (x + zerocoord[1] - pt[1]) / GLB_X +
                                         ((double)py) * (y + zerocoord[2] - pt[2]) / GLB_Y +
                                         ((double)pz) * (z + zerocoord[3] - pt[3]) / GLB_Z);
                                phase = cexp(I * pdotx);

                                _spinor_g5_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][0][tc], phase, tmp);

                                _spinor_g1_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][1][tc], phase, tmp);

                                _spinor_g2_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][2][tc], phase, tmp);

                                _spinor_g3_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][3][tc], phase, tmp);

                                _spinor_g0g5_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][4][tc], phase, tmp);

                                _spinor_g0g1_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][5][tc], phase, tmp);

                                _spinor_g0g2_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][6][tc], phase, tmp);

                                _spinor_g0g3_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][7][tc], phase, tmp);

                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), *_FIELD_AT(prop, ix));
                                _complex_mul_assign(corr[ip][8][tc], phase, tmp);

                                _spinor_g5g1_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][9][tc], phase, tmp);

                                _spinor_g5g2_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][10][tc], phase, tmp);

                                _spinor_g5g3_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][11][tc], phase, tmp);

                                _spinor_g0_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][12][tc], phase, tmp);

                                _spinor_g5g0g1_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][13][tc], phase, tmp);

                                _spinor_g5g0g2_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][14][tc], phase, tmp);

                                _spinor_g5g0g3_f(tmp_spinor, *_FIELD_AT(prop, ix));
                                _spinor_prod_f(tmp, *_FIELD_AT(source, ix), tmp_spinor);
                                _complex_mul_assign(corr[ip][15][tc], phase, tmp);

                                ip = ip + 1;
                            }
                        }
                    }

                } /* end loop t */
            }
        }
    } /* end loop px,py,pz */

    int iGamma;
    for (t = 0; t < T; t++) {
        tc = (zerocoord[0] + t + GLB_T) % GLB_T + 0 * GLB_T + offset;
        for (j = 0; j < n_mom_tot; ++j) {
            for (iGamma = 0; iGamma < NGamma; iGamma++) {
                if (iGamma != 0 && iGamma != 1 && iGamma != 2 && iGamma != 3 && iGamma != 8 && iGamma != 12) {
                    /* multiply by -i by convention [to have hermitean operators ] */
                    corr[j][iGamma][tc] = (-I) * corr[j][iGamma][tc];
                }

                corr_re[j][iGamma][tc] = creal(corr[j][iGamma][tc]);
                corr_im[j][iGamma][tc] = cimag(corr[j][iGamma][tc]);
            }
        }
    }

    /* print into output file */
    for (j = 0; j < n_mom_tot; ++j) {
        for (iGamma = 0; iGamma < NGamma; iGamma++) {
            global_sum(corr_re[j][iGamma], GLB_T);
            global_sum(corr_im[j][iGamma], GLB_T);
            global_sum((double *)(corr[j][iGamma]), 2 * GLB_T);
        }
    }

    if (src_id == 0) { lprintf("CORR", 0, "loops for one noise vector, all local bilinear [t iGamma iSource Re Im] \n"); }

    for (t = 0; t < GLB_T; ++t) {
        for (iGamma = 0; iGamma < NGamma; iGamma++) {
            ip = 0;
            for (px = 0; px < n_mom; ++px) {
                for (py = 0; py < n_mom; ++py) {
                    for (pz = 0; pz < n_mom; ++pz) {
                        lprintf("CORR", 0, "%i %i %i %i %i %i %3.10e %3.10e \n", t, iGamma, src_id, px, py, pz,
                                corr_re[ip][iGamma][t], corr_im[ip][iGamma][t]);

                        if (swc == STORE) {
                            int idx[5] = { src_id, ip, iGamma, t, 0 };
                            *data_storage_element(*ret, 0, idx) = corr_re[ip][iGamma][t];
                            idx[4] = 1;
                            *data_storage_element(*ret, 0, idx) = corr_im[ip][iGamma][t];
                        }
                    }
                }
            }
            //	out_corr[ip][iGamma][t] += creal(corr[ip][iGamma][t]) + I * cimag(corr[ip][iGamma][t]);
            ip = ip + 1;
        }
    }
    fflush(stdout);
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("TIMING", 0, "Contractions for source %i done [%ld sec %ld usec]\n", src_id, etime.tv_sec, etime.tv_usec);

    for (i = 0; i < n_mom_tot; i++) {
        for (j = 0; j < NGamma; j++) {
            free(corr[i][j]);
        }
    }

    for (i = 0; i < n_mom_tot; i++) {
        free(corr[i]);
    }

    free(corr);

    for (j = 0; j < n_mom_tot; ++j) {
        for (i = 0; i < NGamma; ++i) {
            free(corr_re[j][i]);
            free(corr_im[j][i]);
        }
    }
}

#endif