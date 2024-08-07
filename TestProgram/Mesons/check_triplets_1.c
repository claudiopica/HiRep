
/******************************************************************************
 *
 * File check_triplets_1.c
 *
 * Check of the triplet mesons (free case)
 *
 * Author: Agostino Patella
 *
 ******************************************************************************/

#include "libhr.h"
#include <string.h>

//
//#error "Old version of Mesons, it should be updated"

enum MesonT { A = 0, Pi, Rho, B, Pi2, Rho2, Xt, Yt };
static int gid[8] = { 1, -1, -1, -1, 1, 1, 1, -1 }; /* gid = tr \bar{Gamma} Gamma / 4 */
static int g0[8] = { 1, 1, 1, -1, -1, -1, 1, -1 }; /* g0 = tr Gamma Gamma / 4 */
static int g1[8] = { 1, 1, -1, -1, 1, -1, -1, +1 }; /* g1 = tr gamma_1 \bar{Gamma} gamma_1 Gamma / 4 */
static int g2[8] = { 1, 1, 1, 1, 1, 1, -1, -1 }; /* g2 = tr gamma_2 \bar{Gamma} gamma_2 Gamma / 4 */
static int g3[8] = { 1, 1, 1, 1, 1, 1, -1, -1 }; /* g3 = tr gamma_3 \bar{Gamma} gamma_3 Gamma / 4 */
static int gb[8] = { 1, 1, 1, -1, -1, -1, 1, -1 };
static char nameT[8][256] = { "a", "pi", "rho", "b", "pi2", "rho2", "forbidden triplet 0+-", "forbidden triplet 1++" };

/* Mesons parameters */
typedef struct input_mesons {
    char mstring[256];
    double csw;
    /* for the reading function */
    input_record_t read[3];

} input_mesons;

#define init_input_mesons(varname)                                            \
    {                                                                         \
        .read = {                                                             \
            { "fermion mass", "mes:mass = %s", STRING_T, (varname).mstring }, \
            { "csw", "mes:csw = %lf", DOUBLE_T, &(varname).csw },             \
            { NULL, NULL, INT_T, NULL }                                       \
        }                                                                     \
    }

static double mass;
void free_correlators(double **triplets);

input_glb glb_ip = init_input_glb(glb_ip);
input_mesons mes_ip = init_input_mesons(mes_ip);

int main(int argc, char *argv[]) {
    int i, t;
    int return_value = 0;
    int g[4];
    double **ex_triplets;
    double **pta_triplets;
    char pame[256];
    double tol = 1.e-10;
    spinor_field **pta_qprop = NULL;

    /* setup process id and communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    read_input(mes_ip.read, get_input_filename());

    strcpy(pame, mes_ip.mstring);
    mass = atof(strtok(pame, ";"));

    lprintf("MAIN", 0, "mes:masses = %f\n", mass);
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&mes_ip.csw);
#endif

    unit_u(u_gauge);
    represent_gauge_field();
#ifdef REPR_FUNDAMENTAL
    apply_BCs_on_represented_gauge_field(); //This is a trick: the BCs are not applied in the case the REPR is fundamental because represent_gauge field assumes that the right BCs are already applied on the fundamental field!
#endif

    ex_triplets = (double **)malloc(sizeof(double *) * 8);
    for (i = 0; i < 8; i++) {
        ex_triplets[i] = (double *)calloc(GLB_T, sizeof(double *));
    }
    pta_triplets = (double **)malloc(sizeof(double *) * 8);
    for (i = 0; i < 8; i++) {
        pta_triplets[i] = (double *)calloc(GLB_T, sizeof(double *));
    }

    pta_qprop = (spinor_field **)malloc(sizeof(spinor_field *));
    pta_qprop[0] = alloc_spinor_field(4 * NF, &glattice);

    /* CALCOLO ESPLICITO */
    free_correlators(ex_triplets);

    /* MESONI CON PROPAGATORE POINT-TO-ALL */

    g[0] = g[1] = g[2] = g[3] = 0;
    pta_qprop_QMR_eo(g, pta_qprop, 1, &mass, 1e-28);

#ifdef WITH_GPU
    for (int k = 0; k < 4 * NF; k++) {
        copy_from_gpu(pta_qprop[0] + k);
    }
#endif

    id_correlator(pta_triplets[A], g[0], pta_qprop[0]);
    g0_correlator(pta_triplets[Xt], g[0], pta_qprop[0]);
    g5_correlator(pta_triplets[Pi], g[0], pta_qprop[0]);
    g0g5_correlator(pta_triplets[Pi2], g[0], pta_qprop[0]);
    g1_correlator(pta_triplets[Rho], g[0], pta_qprop[0]);
    g0g1_correlator(pta_triplets[Rho2], g[0], pta_qprop[0]);
    g5g1_correlator(pta_triplets[Yt], g[0], pta_qprop[0]);
    g0g5g1_correlator(pta_triplets[B], g[0], pta_qprop[0]);

    /* STAMPA */
    lprintf("TEST", 0, "\nANALITICO\tPOINT-TO-ALL\tERROR (must be less than %e)\n", tol);
    for (i = 0; i < 8; i++) {
        lprintf("TEST", 0, "TRIPLET CORRELATOR %s\n", nameT[i]);
        for (t = 0; t < GLB_T; t++) {
            lprintf("TEST", 0, "%e\t%e\t%e\n", ex_triplets[i][t], pta_triplets[i][t],
                    fabs(ex_triplets[i][t] - pta_triplets[i][t]));
            if (fabs(ex_triplets[i][t] - pta_triplets[i][t]) > tol) { return_value += 1; }
        }
    }

    global_sum_int(&return_value, 1);
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);

    for (i = 0; i < 8; i++) {
        free(pta_triplets[i]);
        free(ex_triplets[i]);
    }
    free(pta_triplets);
    free(ex_triplets);
    free(pta_qprop);

    finalize_process();
    return return_value;
}

double square(double x) {
    return x * x;
}

hr_complex ev(double k[4]) {
    hr_complex z;
    z = mass + 4.0 + cos((2.0 * M_PI * k[0]) / GLB_T) + cos((2.0 * M_PI * k[1]) / GLB_X) + cos((2.0 * M_PI * k[2]) / GLB_Y) +
        cos((2.0 * M_PI * k[3]) / GLB_Z) +
        I * sqrt(square(sin((2.0 * M_PI * k[0]) / GLB_T)) + square(sin((2.0 * M_PI * k[1]) / GLB_X)) +
                 square(sin((2.0 * M_PI * k[2]) / GLB_Y)) + square(sin((2.0 * M_PI * k[3]) / GLB_Z)));
    return z;
}

void free_correlators(double **triplets) {
    double A2[GLB_T], B2[4][GLB_T];
    hr_complex AA[GLB_T], BB[4][GLB_T];
    hr_complex tmp, eit;
    double norm2, z;
    int i, j, t;
    double k[4];
    double sigma[4] = { 0., 0., 0., 0. };
#ifdef BC_T_ANTIPERIODIC
    sigma[0] = .5;
#endif
#ifdef BC_X_ANTIPERIODIC
    sigma[1] = .5;
#endif
#ifdef BC_Y_ANTIPERIODIC
    sigma[2] = .5;
#endif
#ifdef BC_Z_ANTIPERIODIC
    sigma[3] = .5;
#endif

    z = -4. * NF / square(GLB_T * GLB_X * GLB_Y * GLB_Z);

    lprintf("FREE", 0, "sigma = (%f,%f,%f,%f)\n", sigma[0], sigma[1], sigma[2], sigma[3]);
    for (t = 0; t < GLB_T; t++) {
        A2[t] = 0.;
        B2[0][t] = B2[1][t] = B2[2][t] = B2[3][t] = 0.;
    }

    for (k[1] = sigma[1]; k[1] < GLB_X + sigma[1] - .5; k[1] += 1.) {
        for (k[2] = sigma[2]; k[2] < GLB_Y + sigma[2] - .5; k[2] += 1.) {
            for (k[3] = sigma[3]; k[3] < GLB_Z + sigma[3] - .5; k[3] += 1.) {
                for (t = 0; t < GLB_T; t++) {
                    AA[t] = 0.;
                    BB[0][t] = 0.;
                    BB[1][t] = 0.;
                    BB[2][t] = 0.;
                    BB[3][t] = 0.;
                }
                for (k[0] = sigma[0]; k[0] < GLB_T + sigma[0] - .5; k[0] += 1.) {
                    tmp = ev(k);
                    norm2 = tmp * conj(tmp);

                    for (t = 0; t < GLB_T; t++) {
                        eit = cos((2.0 * M_PI * t * k[0]) / GLB_T) + I * sin((2.0 * M_PI * t * k[0]) / GLB_T);
                        AA[t] += creal(tmp) * eit / norm2;
                        BB[0][t] += sin((2.0 * M_PI * k[0]) / GLB_T) * eit / norm2;
                        for (j = 1; j < 4; j++) {
                            BB[j][t] += sin((2.0 * M_PI * k[j]) / GLB_X) * eit / norm2;
                        }
                    }
                }
                for (t = 0; t < GLB_T; t++) {
                    A2[t] += AA[t] * conj(AA[t]) * z;
                    for (j = 0; j < 4; j++) {
                        B2[j][t] += BB[j][t] * conj(BB[j][t]) * z;
                    }
                }
            }
        }
    }

    lprintf("FREE", 10, "A2\tB2[0]\tB2[1]\tB2[2]\tB[3]\n", A2[0]);
    for (t = 0; t < GLB_T; t++) {
        lprintf("FREE", 10, "%e\t%e\t%e\t%e\t%e\n", A2[t], B2[0][t], B2[1][t], B2[2][t], B2[3][t]);
    }

    for (i = 0; i < 8; i++) {
        for (t = 0; t < GLB_T; t++) {
            triplets[i][t] =
                gb[i] * (gid[i] * A2[t] - g0[i] * B2[0][t] - g1[i] * B2[1][t] - g2[i] * B2[2][t] - g3[i] * B2[3][t]);
        }
    }

    lprintf("FREE", 0, "Exact free correlators computed.\n");
}
