
/******************************************************************************
*
*
* File check_disc_5.c
*
* Check of the  disc loops (free case): discon volume source (type = 5)
*
* Author: Vincent Drach
*
* NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
* NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
*
******************************************************************************/

#include "libhr.h"
#include <string.h>

static hr_complex gid[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; /* gid = tr Gamma  */
static hr_complex g0[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; /*  g0 = tr gamma_0 Gamma */
static hr_complex g1[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; /*  g1 = tr gamma_1 Gamma */
static hr_complex g2[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; /*  g2 = tr gamma_2 Gamma */
static hr_complex g3[16] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; /*  g3 = tr gamma_3 Gamma */

char *mes_channel_names[16] = { "g5", "g1",     "g2",     "g3",     "-ig0g5", "-ig0g1",   "-ig0g2",   "-ig0g3",
                                "id", "-ig5g1", "-ig5g2", "-ig5g3", "g0",     "-ig5g0g1", "-ig5g0g2", "-ig5g0g3" };

#define mult_mat(r, A, B)                                \
    {                                                    \
        int _i, _j, _k;                                  \
        hr_complex wm[4][4];                             \
        for (_i = 0; _i < 4; _i++)                       \
            for (_j = 0; _j < 4; _j++) {                 \
                _complex_0(wm[_i][_j]);                  \
                for (_k = 0; _k < 4; _k++) {             \
                    wm[_i][_j] += A[_i][_k] * B[_k][_j]; \
                }                                        \
            }                                            \
        for (_i = 0; _i < 4; _i++)                       \
            for (_j = 0; _j < 4; _j++) {                 \
                r[_i][_j] = wm[_i][_j];                  \
            }                                            \
    }

#define mult_mat_minusI(r, A)               \
    {                                       \
        int _i, _j;                         \
        for (_i = 0; _i < 4; _i++)          \
            for (_j = 0; _j < 4; _j++) {    \
                r[_i][_j] = -I * A[_i][_j]; \
            }                               \
    }

#define set_zero_mat(A)                  \
    {                                    \
        int _i, _j;                      \
        for (_i = 0; _i < 4; _i++)       \
            for (_j = 0; _j < 4; _j++) { \
                A[_i][_j] = 0.;          \
            }                            \
    }
#define trace_mat(r, A)              \
    {                                \
        int _i;                      \
        _complex_0(r);               \
        for (_i = 0; _i < 4; _i++) { \
            r += A[_i][_i];          \
        }                            \
    }

/*VD: This is kept to to some debugging.*/
/*  static void print_mat(hr_complex mat[4][4], const char name[]) {
int i,j;
lprintf("MAIN",0,"%s = \n", name);
for(i=0; i<4; i++) {
lprintf("MAIN",0,"[ ");
for(j=0; j<4; j++) {
lprintf("MAIN",0,"(%.2f,%.2f) ",creal(mat[i][j]),cimag(mat[i][j]));
}
lprintf("MAIN",0,"]\n");
}
}*/
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
            { "Fermion mass", "disc:mass = %s", STRING_T, (varname).mstring },                \
            { "inverter precision", "disc:precision = %lf", DOUBLE_T, &(varname).precision }, \
            { "number of inversions per cnfg", "disc:nhits5 = %d", INT_T, &(varname).nhits }, \
            { "maximum component of momentum", "disc:n_mom = %d", INT_T, &(varname).n_mom },  \
            { "csw coefficient", "mes:csw = %lf", DOUBLE_T, &(varname).csw },                 \
            { NULL, NULL, INT_T, NULL }                                                       \
        }                                                                                     \
    }

static double mass;
void free_loops(hr_complex *loops);

input_mesons mes_ip = init_input_mesons(mes_ip);

char char_t[100];
FILE *fp;
char path[1035];
/*Gamma / 4 */
hr_complex get_gid(hr_complex Gamma[4][4]) {
    hr_complex r;
    trace_mat(r, Gamma);
    //printf("get_gid %f %f \n",creal(r)/4.,cimag(r)/4.);
    return r;
}
/* gmu = tr Gamma gamma_mu   */
hr_complex get_gmu(hr_complex Gamma[4][4], int mu) {
    int sign;
    hr_complex tmp[4][4];
    hr_complex r;
    hr_complex gmu[4][4];

    if (mu == 1) { g1_debug(gmu, &sign); }
    if (mu == 2) { g2_debug(gmu, &sign); }
    if (mu == 3) { g3_debug(gmu, &sign); }
    if (mu == 0) { g0_debug(gmu, &sign); }

    //print_mat(Gamma,"Gamma");
    //  print_mat(gmu,"gmu");

    mult_mat(tmp, Gamma, gmu);

    trace_mat(r, tmp);
    //printf("get_gmu mu=%d %f %f \n",mu,creal(r)/4.,cimag(r)/4.);
    return r;
}

int compare_disc(hr_complex *corr_ex, hr_complex *corr_num, char *name[16], double tol, double tol_rel_scalar) {
    int retval = 0;
    int nGamma = 16;
    for (int n = 0; n < nGamma; n++) {
        if (n == 8) {
            if (cabs(corr_ex[8] - corr_num[8]) / cabs(corr_ex[8]) > tol_rel_scalar) {
                lprintf("TEST", 0, "Mismatch %s, rel diff: %e, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name[n],
                        cabs(corr_ex[n] - corr_num[n]) / cabs(corr_ex[n]), creal(corr_num[n]), cimag(corr_num[n]),
                        creal(corr_ex[n]), cimag(corr_ex[n]));
                retval += 1;
            } else {
                lprintf("TEST", 0, "Match %s, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name[n], creal(corr_num[n]),
                        cimag(corr_num[n]), creal(corr_ex[n]), cimag(corr_ex[n]));
            }
        } else {
            if (cabs(corr_ex[n] - corr_num[n]) > tol) {
                lprintf("TEST", 0, "Mismatch %s, absolute diff: %e, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name[n],
                        cabs(corr_ex[n] - corr_num[n]), creal(corr_num[n]), cimag(corr_num[n]), creal(corr_ex[n]),
                        cimag(corr_ex[n]));
                retval += 1;
            } else {
                lprintf("TEST", 0, "Match %s, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name[n], creal(corr_num[n]),
                        cimag(corr_num[n]), creal(corr_ex[n]), cimag(corr_ex[n]));
            }
        }
    }
    return retval;
}

int main(int argc, char *argv[]) {
    int i, sign;
    hr_complex *ex_loops;
    char pame[256];
    int n_Gamma = 16;
    int source_type = 5;
    int return_value = 0;
    data_storage_array *out_corr = NULL;
    hr_complex *mean_loops;
    double abs_tol = 1e-1;
    double rel_tol_scalar_loop = 5.e-3;
    struct timeval start, end, etime;
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary

    hr_complex g[16][4][4];
    hr_complex tmp[4][4];
    g5_debug(g[0], &sign);
    g1_debug(g[1], &sign);
    g2_debug(g[2], &sign);
    g3_debug(g[3], &sign);
    g0g5_debug(tmp, &sign);
    mult_mat_minusI(g[4], tmp);
    g0g1_debug(tmp, &sign);
    mult_mat_minusI(g[5], tmp);
    g0g2_debug(tmp, &sign);
    mult_mat_minusI(g[6], tmp);
    g0g3_debug(tmp, &sign);
    mult_mat_minusI(g[7], tmp);
    id_debug(g[8], &sign);
    g5g1_debug(tmp, &sign);
    mult_mat_minusI(g[9], tmp);
    g5g2_debug(tmp, &sign);
    mult_mat_minusI(g[10], tmp);
    g5g3_debug(tmp, &sign);
    mult_mat_minusI(g[11], tmp);
    g0_debug(g[12], &sign);
    mult_mat(g[13], g[0], g[5]); //  -i g5g0g1
    mult_mat(g[14], g[0], g[6]); //  -i g5g0g2
    mult_mat(g[15], g[0], g[7]); //  -i g5g0g3

    for (i = 0; i < 16; i++) {
        gid[i] = get_gid(g[i]);
        g0[i] = get_gmu(g[i], 0);
        g1[i] = get_gmu(g[i], 1);
        g2[i] = get_gmu(g[i], 2);
        g3[i] = get_gmu(g[i], 3);
    }

    /* setup process id and communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    read_input(mes_ip.read, get_input_filename());

    strcpy(pame, mes_ip.mstring);
    mass = atof(strtok(pame, ";"));

    lprintf("MAIN", 0, "disc:mass = %f\n", mass);
    lprintf("MAIN", 0, "disc:nhits = %i\n", mes_ip.nhits);
    lprintf("MAIN", 0, "Inverter precision = %e\n", mes_ip.precision);
    lprintf("MAIN", 0, "Number of momenta = %d\n", mes_ip.n_mom);
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&mes_ip.csw);
#endif

    gettimeofday(&start, 0);

    unit_u(u_gauge);
    represent_gauge_field();
#ifdef REPR_FUNDAMENTAL
    apply_BCs_on_represented_gauge_field(); //This is a trick: the BCs are not applied in the case the REPR is fundamental because represent_gauge field assumes that the right BCs are already applied on the fundamental field!
#endif

    lprintf("MAIN", 0, "source type is fixed to 5: Spin , color and eo dilution \n");

    lprintf("MAIN", 0, "Measuring D(t) =  sum_x psibar(x) Gamma psi(x)\n");

    lprintf("MAIN", 0, "Zerocoord{%d,%d,%d,%d}\n", zerocoord[0], zerocoord[1], zerocoord[2], zerocoord[3]);

    error(!(GLB_X == GLB_Y && GLB_X == GLB_Z), 1, "main", "This test works only for GLB_X=GLB_Y=GLB_Z");

    lprintf("CORR", 0, "Number of noise vector : nhits = %i \n", mes_ip.nhits);

    measure_loops(&mass, mes_ip.nhits, 0, mes_ip.precision, source_type, mes_ip.n_mom, STORE, &out_corr);

    //stochastic & time average
    mean_loops = (hr_complex *)calloc(n_Gamma, sizeof(hr_complex));
    for (int k = 0; k < mes_ip.nhits; k++) {
        for (int eo = 0; eo < 2; eo++) {
            for (int col = 0; col < NF; col++) {
                for (int j = 0; j < n_Gamma; j++) {
                    for (int t = 0; t < GLB_T; t++) {
                        int idx_re[6] = { k, eo, col, j, t, 0 };
                        int idx_im[6] = { k, eo, col, j, t, 1 };

                        mean_loops[j] +=
                            (*data_storage_element(out_corr, 0, idx_re) + I * *data_storage_element(out_corr, 0, idx_im)) /
                            (mes_ip.nhits * GLB_T);
                    }
                }
            }
        }
    }
    //for(int j=0; j<n_Gamma; j++)   for(int t=0; t<GLB_T; t++) mean_loops[j] += out_corr[0][j][t]/(mes_ip.nhits*GLB_T);

    //  /* CALCOLO ESPLICITO */
    ex_loops = (hr_complex *)calloc(16, sizeof(hr_complex));
    free_loops(ex_loops);

    return_value += compare_disc(ex_loops, mean_loops, mes_channel_names, abs_tol, rel_tol_scalar_loop);

    global_sum_int(&return_value, 1);
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Configuration : analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

    finalize_process();
    return return_value;
}

double square(double x) {
    return x * x;
}

double denom(double k[4]) {
    double res;
    res = square(mass + 4.0 - cos((2.0 * M_PI * k[0]) / GLB_T) - cos((2.0 * M_PI * k[1]) / GLB_X) -
                 cos((2.0 * M_PI * k[2]) / GLB_Y) - cos((2.0 * M_PI * k[3]) / GLB_Z)) +
          square(sin((2.0 * M_PI * k[0]) / GLB_T)) + square(sin((2.0 * M_PI * k[1]) / GLB_X)) +
          square(sin((2.0 * M_PI * k[2]) / GLB_Y)) + square(sin((2.0 * M_PI * k[3]) / GLB_Z));
    return res;
}

void free_loops(hr_complex *loops) {
    double A, B[4];
    double tmp, norm;
    int i;
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

    lprintf("FREE", 0, "sigma = (%f,%f,%f,%f)\n", sigma[0], sigma[1], sigma[2], sigma[3]);
    norm = (double)NF / GLB_T;

    A = B[0] = B[1] = B[2] = B[3] = 0.;

    for (k[1] = sigma[1]; k[1] < GLB_X + sigma[1] - .5; k[1] += 1.) {
        for (k[2] = sigma[2]; k[2] < GLB_Y + sigma[2] - .5; k[2] += 1.) {
            for (k[3] = sigma[3]; k[3] < GLB_Z + sigma[3] - .5; k[3] += 1.) {
                for (k[0] = sigma[0]; k[0] < GLB_T + sigma[0] - .5; k[0] += 1.) {
                    tmp = denom(k);

                    A += (mass + 4.0 - cos((2.0 * M_PI * k[0]) / GLB_T) - cos((2.0 * M_PI * k[1]) / GLB_X) -
                          cos((2.0 * M_PI * k[2]) / GLB_Y) - cos((2.0 * M_PI * k[3]) / GLB_Z)) /
                         tmp;
                    B[0] += sin((2.0 * M_PI * k[0]) / GLB_T) / tmp;
                    B[1] += sin((2.0 * M_PI * k[1]) / GLB_X) / tmp;
                    B[2] += sin((2.0 * M_PI * k[2]) / GLB_Y) / tmp;
                    B[3] += sin((2.0 * M_PI * k[3]) / GLB_Z) / tmp;
                }
            }
        }
    }

    lprintf("FREE", 10, "A=%e\t B[0=]%e\t B[1]=%e\t B[2]=%e\t B[3]=%e\n", A, B[0], B[1], B[2], B[3]);

    for (i = 0; i < 16; i++) {
        loops[i] = (gid[i] * A - I * (g0[i] * B[0] + g1[i] * B[1] + g2[i] * B[2] + g3[i] * B[3])) *
                   norm; // = sum_vec{x} psi(x) Gamma psibar(x)
    }

    lprintf("FREE", 0, "Exact free correlators computed.\n");
}
