/******************************************************************************
*
* File check_sources.c
*
* Checks of the stochastic sources with 4 spinor fields
*
* Author: Vincent Drach
*
******************************************************************************/

#include "libhr.h"

/* Mesons parameters */
typedef struct input_mesons {
    char mstring[256];

    /* for the reading function */
    input_record_t read[7];
    double precision;
    int nhits;
    int source_type;
    int n_mom;
    double csw;

} input_mesons;

#define init_input_mesons(varname)                                                            \
    {                                                                                         \
        .read = {                                                                             \
            { "quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring },      \
            { "inverter precision", "disc:precision = %lf", DOUBLE_T, &(varname).precision }, \
            { "number of inversions per cnfg", "disc:nhits = %d", INT_T, &(varname).nhits },  \
            { "maximum component of momentum", "disc:n_mom = %d", INT_T, &(varname).n_mom },  \
            { NULL, NULL, INT_T, NULL }                                                       \
        }                                                                                     \
    }

input_glb glb_ip = init_input_glb(glb_ip);
input_mesons mes_ip = init_input_mesons(mes_ip);

char char_t[100];
FILE *fp;
char path[1035];
int i;

static hr_complex DeltaKronecker(int l, int j) {
    if (l == j) {
        return 1.;
    } else {
        return 0.;
    }
}
static hr_complex is_eo(int t, int x, int y, int z, int eo) {
    if (((zerocoord[0] + t + zerocoord[1] + x + zerocoord[2] + y + zerocoord[3] + z) & 1) == eo) {
        return 1.;
    } else {
        return 0.;
    }
}

static hr_complex average_complex(hr_complex a[], int n) {
    hr_complex sum;
    _complex_0(sum);
    for (int l = 0; l < n; l++) {
        sum += a[l];
    }
    return sum / n;
}

static hr_complex sd(hr_complex a[], int n) {
    hr_complex mean;
    mean = average_complex(a, n);

    hr_complex sum = 0;
    for (int l = 0; l < n; l++) {
        sum += creal(a[l] - mean) * creal(a[l] - mean) + I * cimag(a[l] - mean) * cimag(a[l] - mean);
    }

    // return the standard deviation of the real and imaginary part in a complex number;
    return sqrt(creal(sum) / (n - 1)) + I * sqrt(cimag(sum) / (n - 1));
}

static int test_passed(double mean, double stdev) {
    if (stdev < 0.) { return 1; }

    if (fabs(mean) < 1e-10 && fabs(stdev) < 1e-10) { return 0; }
    if (mean < 0. && mean + stdev > 0.) { return 0; }
    if (mean >= 0. && mean - stdev < 0.) { return 0; }

    return 1;
}

int main(int argc, char *argv[]) {
    int return_value = 0;
    int x, y, z, t, ix;
    int rand_block[4];
    int shift[4];
    int b1, b2, c1, c2;
    int col, tau, eo, j;
    int k, counter;
    hr_complex av_, sd_;
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary

    logger_map("DEBUG", "debug");
    logger_setlevel(0, 200);

    /* setup process id and communications */
    setup_process(&argc, &argv);
    read_input(mes_ip.read, get_input_filename());

    tau = rand() % GLB_T;
    bcast_int(&tau, 1);

    lprintf("MAIN", 0, "This code test the following functions:\n");
    lprintf("MAIN", 0, "create_diluted_source_equal_atau(source,tau)\n");
    lprintf("MAIN", 0, "create_diluted_source_equal_atau_col(source,tau,col)\n");
    lprintf("MAIN", 0, "create_diluted_source_equal_atau_col(source,tau,col) + eo \n");

    lprintf("MAIN", 0, "Number of noise vector : nhits = %i \n", mes_ip.nhits);
    lprintf("MAIN", 0, "Chosen random timeslice = %i \n", tau);

    double norm;

    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *source_shifted = alloc_spinor_field(4, &glattice);

    lprintf("CORR", 0, "Testing spinorfield timeslice source without dilution .... \n");

    for (i = 0; i < mes_ip.nhits; i++) {
        create_diluted_source_equal_atau(source, tau);

        for (j = 0; j < 4; ++j) {
            norm = sqnorm_spinor_field(&source[j]);

            if (fabs(norm / ((double)NF * GLB_X * GLB_Y * GLB_Z) - 1.) > 1e-14) { return_value += 1; }
        } // loop on j
    }
    global_sum_int(&return_value, 1);
    if (return_value == 0) {
        lprintf("MAIN", 0, "test norm passed\n ");
    } else {
        lprintf("MAIN", 0, "test norm FAILED: return_value= %d\n ", return_value);
    }

    lprintf("CORR", 0, "Testing spinorfield timeslice source with spin dilution .... \n");
    for (col = 0; col < NF; ++col) {
        for (i = 0; i < mes_ip.nhits; i++) {
            create_diluted_source_equal_atau_col(source, tau, col);

            for (j = 0; j < 4; ++j) {
                norm = sqnorm_spinor_field(&source[j]);

                if (fabs(norm / ((double)GLB_X * GLB_Y * GLB_Z) - 1.) > 1e-14) { return_value += 1; }
            }
        }
    }
    global_sum_int(&return_value, 1);
    if (return_value == 0) {
        lprintf("MAIN", 0, "test norm passed\n ");
    } else {
        lprintf("MAIN", 0, "test norm FAILED: return_value= %d\n ", return_value);
    }

    lprintf("CORR", 0, "Testing spinorfield timeslice source with spin dilution and eo dilution .... \n");
    // spin dilution only.
    for (col = 0; col < NF; ++col) {
        for (i = 0; i < mes_ip.nhits; i++) {
            for (eo = 0; eo < 2; ++eo) {
                create_diluted_source_equal_atau_col(source, tau, col);
                zero_even_or_odd_site_spinorfield(source, 4, eo);

                for (j = 0; j < 4; ++j) {
                    norm = sqnorm_spinor_field(&source[j]);
                    if (fabs(norm / (0.5 * (double)GLB_X * GLB_Y * GLB_Z) - 1.) > 1e-14) { return_value += 1; }
                } // loop on j
            } // loop on eo.
        }
    }
    global_sum_int(&return_value, 1);
    if (return_value == 0) {
        lprintf("MAIN", 0, "test norm passed\n ");
    } else {
        lprintf("MAIN", 0, "test norm FAILED: return_value= %d\n ", return_value);
    }

    hr_complex tmp_vec[mes_ip.nhits];

    for (int ll = 0; ll < 10 * MPI_WORLD_SIZE; ll++) {
        if (ll == 0) {
            rand_block[0] = 0;
            rand_block[1] = 0;
            rand_block[2] = 0;
            rand_block[3] = 0;
        } else {
            rand_block[0] = rand() % (NP_T);
            rand_block[1] = rand() % (NP_X);
            rand_block[2] = rand() % (NP_Y);
            rand_block[3] = rand() % (NP_Z);
        }
        bcast_int(rand_block, 4);
        col = rand() % NF;
        bcast_int(&col, 1);

        int tmp_index = rand_block[0] + rand_block[1] + rand_block[2] + rand_block[3];
        shift[0] = rand_block[0] * T;
        shift[1] = rand_block[1] * X;
        shift[2] = rand_block[2] * Y;
        shift[3] = rand_block[3] * Z;

        lprintf("MAIN", 0, "Random block of coord: (%d %d %d %d) \n", rand_block[0], rand_block[1], rand_block[2],
                rand_block[3]);
        lprintf("MAIN", 0, "shift: (%d %d %d %d) \n", shift[0], shift[1], shift[2], shift[3]);

        lprintf("MAIN", 0, "Now testing create_diluted_source_equal_atau \n");

        for (i = 0; i < mes_ip.nhits; i++) {
            _complex_0(tmp_vec[i]);
        }

        counter = 0;
        for (i = 0; i < mes_ip.nhits; i++) {
            create_diluted_source_equal_atau(source, tau);
            for (k = 0; k < 4; k++) {
                shift_fields(shift, &source[k], NULL, &source_shifted[k], NULL);
            }

            for (x = 0; x < X; x += 1) {
                for (y = 0; y < Y; y += 1) {
                    for (z = 0; z < Z; z += 1) {
                        for (t = 0; t < T; t += 1) {
                            ix = ipt(t, x, y, z);
                            for (b1 = 0; b1 < 4; ++b1) {
                                for (c1 = 0; c1 < NF; ++c1) {
                                    for (b2 = 0; b2 < 4; ++b2) {
                                        for (c2 = 0; c2 < NF; ++c2) {
                                            for (k = 0; k < 4; k++) {
                                                tmp_vec[counter] +=
                                                    conj(_FIELD_AT(&source[k], ix)->c[b1].c[c1]) *
                                                        _FIELD_AT(&source_shifted[k], ix)->c[b2].c[c2] -
                                                    DeltaKronecker(tau - zerocoord[0], t) * DeltaKronecker(k, b2) *
                                                        DeltaKronecker(k, b1) * DeltaKronecker(tmp_index, 0) *
                                                        DeltaKronecker(b1, b2) *
                                                        DeltaKronecker(
                                                            c1, c2); // that quantity should be 0 in average over the noise.
                                            }
                                        }
                                    }
                                }
                            }
                        } //loop local volume
                    }
                }
            }
            counter += 1;
        } //loop nhits

        av_ = average_complex(tmp_vec, mes_ip.nhits);
        sd_ = sd(tmp_vec, mes_ip.nhits);

        if (test_passed(creal(av_), creal(sd_)) == 1 || test_passed(cimag(av_), cimag(sd_)) == 1) { return_value += 1; }

        if (test_passed(creal(av_), creal(sd_)) == 1 || test_passed(cimag(av_), cimag(sd_)) == 1) {
            lprintf("MAIN", 0, " failed : %d , %e(%e)   %e(%e)  \n", test_passed(creal(av_), creal(sd_)), creal(av_),
                    creal(sd_), cimag(av_), cimag(sd_));
        }

        lprintf("MAIN", 0, "Now testing create_diluted_source_equal_atau_col \n");

        for (i = 0; i < mes_ip.nhits; i++) {
            _complex_0(tmp_vec[i]);
        }

        counter = 0;

        for (i = 0; i < mes_ip.nhits; i++) {
            create_diluted_source_equal_atau_col(source, tau, col);

            for (k = 0; k < 4; k++) {
                shift_fields(shift, &source[k], NULL, &source_shifted[k], NULL);
            }

            for (x = 0; x < X; x += 1) {
                for (y = 0; y < Y; y += 1) {
                    for (z = 0; z < Z; z += 1) {
                        for (t = 0; t < T; t += 1) {
                            ix = ipt(t, x, y, z);
                            for (b1 = 0; b1 < 4; ++b1) {
                                for (c1 = 0; c1 < NF; ++c1) {
                                    for (b2 = 0; b2 < 4; ++b2) {
                                        for (c2 = 0; c2 < NF; ++c2) {
                                            for (k = 0; k < 4; k++) {
                                                tmp_vec[counter] +=
                                                    conj(_FIELD_AT(&source[k], ix)->c[b1].c[c1]) *
                                                        _FIELD_AT(&source_shifted[k], ix)->c[b2].c[c2] -
                                                    DeltaKronecker(col, c1) * DeltaKronecker(col, c2) *
                                                        DeltaKronecker(tau - zerocoord[0], t) * DeltaKronecker(k, b2) *
                                                        DeltaKronecker(k, b1) * DeltaKronecker(tmp_index, 0) *
                                                        DeltaKronecker(b1, b2) *
                                                        DeltaKronecker(
                                                            c1, c2); // that quantity should be 0 in average over the noise.
                                            }
                                        }
                                    }
                                }
                            }
                        } //loop local volume
                    }
                }
            }
            counter += 1;
        } //loop nhits

        av_ = average_complex(tmp_vec, mes_ip.nhits);
        sd_ = sd(tmp_vec, mes_ip.nhits);

        if (test_passed(creal(av_), creal(sd_)) == 1 || test_passed(cimag(av_), cimag(sd_)) == 1) { return_value += 1; }

        if (test_passed(creal(av_), creal(sd_)) == 1 || test_passed(cimag(av_), cimag(sd_)) == 1) {
            lprintf("MAIN", 0, " failed : %d , %e(%e)   %e(%e)  \n", test_passed(creal(av_), creal(sd_)), creal(av_),
                    creal(sd_), cimag(av_), cimag(sd_));
        }

        lprintf("MAIN", 0, "Now testing create_diluted_source_equal_atau_col + eo \n");

        for (eo = 0; eo < 2; ++eo) {
            for (i = 0; i < mes_ip.nhits; i++) {
                _complex_0(tmp_vec[i]);
            }

            counter = 0;
            for (i = 0; i < mes_ip.nhits; i++) {
                create_diluted_source_equal_atau_col(source, tau, col);
                zero_even_or_odd_site_spinorfield(source, 4, eo);

                for (k = 0; k < 4; k++) {
                    shift_fields(shift, &source[k], NULL, &source_shifted[k], NULL);
                }

                for (x = 0; x < X; x += 1) {
                    for (y = 0; y < Y; y += 1) {
                        for (z = 0; z < Z; z += 1) {
                            for (t = 0; t < T; t += 1) {
                                ix = ipt(t, x, y, z);
                                for (b1 = 0; b1 < 4; ++b1) {
                                    for (c1 = 0; c1 < NF; ++c1) {
                                        for (b2 = 0; b2 < 4; ++b2) {
                                            for (c2 = 0; c2 < NF; ++c2) {
                                                for (k = 0; k < 4; k++) {
                                                    tmp_vec[counter] +=
                                                        conj(_FIELD_AT(&source[k], ix)->c[b1].c[c1]) *
                                                            _FIELD_AT(&source_shifted[k], ix)->c[b2].c[c2] -
                                                        is_eo(t, x, y, z, eo) * DeltaKronecker(col, c1) *
                                                            DeltaKronecker(col, c2) * DeltaKronecker(tau - zerocoord[0], t) *
                                                            DeltaKronecker(k, b2) * DeltaKronecker(k, b1) *
                                                            DeltaKronecker(tmp_index, 0) * DeltaKronecker(b1, b2) *
                                                            DeltaKronecker(
                                                                c1, c2); // that quantity should be 0 in average over the noise.
                                                }
                                            }
                                        }
                                    }
                                }
                            } //loop local volume
                        }
                    }
                }
                counter += 1;
            } //loop nhits

            av_ = average_complex(tmp_vec, mes_ip.nhits);
            sd_ = sd(tmp_vec, mes_ip.nhits);

            if (test_passed(creal(av_), creal(sd_)) == 1 || test_passed(cimag(av_), cimag(sd_)) == 1) { return_value += 1; }

            if (test_passed(creal(av_), creal(sd_)) == 1 || test_passed(cimag(av_), cimag(sd_)) == 1) {
                lprintf("MAIN", 0, " failed : %d , %e(%e)   %e(%e)  \n", test_passed(creal(av_), creal(sd_)), creal(av_),
                        creal(sd_), cimag(av_), cimag(sd_));
            }

        } // loop eo

    } //loop on ll (the number of shifts)

    global_sum_int(&return_value, 1);
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);

    free_spinor_field(source);
    free_spinor_field(source_shifted);

    finalize_process();

    return return_value;
}
