/******************************************************************************
 *
 *
 * File check_scattering_length.c
 * Checks of the pi pi scattering length calculations 
 * Author: Vincent Drach & Fernando Romero Lopez
 * NOCOMPILE= BC_X_ANTIPERIODIC
 * NOCOMPILE= BC_Y_ANTIPERIODIC
 * NOCOMPILE= BC_Z_ANTIPERIODIC
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= FERMION_THETA
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "scattering.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "clover_tools.h"
#include "setup.h"
#include "cinfo.c"
#include "clover_tools.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define SQR(A) ((A) * (A))
#define CMUL(a, b) \
    (a).re *= b;   \
    (a).im *= b

/// \cond
#define BASENAME(filename) (strrchr((filename), '/') ? strrchr((filename), '/') + 1 : filename)

#define corr_ind(px, py, pz, n_mom, tc, nm, cm) ((px) * (n_mom) * (n_mom)*GLB_T * (nm) + (py) * (n_mom)*GLB_T * (nm) + (pz)*GLB_T * (nm) + ((cm)*GLB_T) + (tc))
/// \endcond

/// \cond
#define INDEX(px, py, pz, n_mom, tc) ((px + n_mom) * (2 * n_mom + 1) * (2 * n_mom + 1) * (GLB_T) + (py + n_mom) * (2 * n_mom + 1) * (GLB_T) + (pz + n_mom) * (GLB_T) + (tc))
/// \endcond

/**
 * @brief Structure containing data from the input file relevant to scattering.
 */
typedef struct _input_scatt
{
    char mstring[256];
    double csw;
    double precision;
    int nhits;

    /* for the reading function */
    input_record_t read[7];

} input_scatt;

#define init_input_scatt(varname)                                                          \
    {                                                                                      \
        .read = {                                                                          \
            {"Fermion masses", "mes:mass = %s", STRING_T, (varname).mstring},              \
            {"csw", "mes:csw = %lf", DOUBLE_T, &(varname).csw},                            \
            {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision}, \
            {"number of inversions per cnfg", "I0:nhits = %d", INT_T, &(varname).nhits},   \
            {NULL, NULL, INT_T, NULL}                                                      \
        }                                                                                  \
    }

char input_filename[256] = "input_file";

input_scatt mes_var = init_input_scatt(mes_var);

typedef struct
{
    char string[256];
    char configlist[256];
    int t, x, y, z;
    int nc, nf;
    double b, m;
    int n;
    int type;
} filename_t;

/**
 * @brief Sets the configuration to unit gauge
 */
void unit_gauge(suNg_field *gauge)
{
    int mu;
    int x, y, z, t, ix;
    for (t = 0; t < T; t++)
        for (x = 0; x < X; x++)
            for (y = 0; y < Y; y++)
                for (z = 0; z < Z; z++)
                {
                    ix = ipt(t, x, y, z);
                    for (mu = 0; mu < 4; ++mu)
                    {
                        _suNg_unit(*_4FIELD_AT(gauge, ix, mu));
                    }
                }
    start_gf_sendrecv(gauge);
    complete_gf_sendrecv(gauge);
}

typedef struct fourvector
{
    double v[4];
} fourvec;

/**
 * @brief Adds two four-vectors together replacing the first one with the sum, v1 += v2
 */
static void iadd(fourvec *v1, fourvec *v2)
{
    for (int i = 0; i < 4; ++i)
    {
        v1->v[i] += v2->v[i];
    }
}

/**
 * @brief Multiply four-vector by a real number
 */
static void imul(fourvec *v1, double a)
{
    for (int i = 0; i < 4; ++i)
    {
        v1->v[i] *= a;
    }
}

/**
 * @brief Returns sum over phat^2 + m (part of the propagator)
 */
double f1(fourvec p, double m)
{
    double tmp = 0.0;
    int i;

    for (i = 0; i < 4; ++i)
    {
        tmp += sin(p.v[i] / 2) * sin(p.v[i] / 2);
    }
    return m + 2 * tmp;
}

/**
 * @brief Part of the propagator
 */
static double f2(fourvec v1, fourvec v2)
{
    int i;
    double result = 0.0;
    for (i = 0; i < 4; ++i)
    {
        result += sin(v1.v[i]) * sin(v2.v[i]);
    }
    return result;
}

/**
 * @brief Part of the propagator
 */
double b_mu(fourvec p1, int mu)
{
    return sin(p1.v[mu]);
}

/**
 * @brief Calculates analytic expression for pion 2-point function
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
double complex twopoint(fourvec p, double m, int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    double complex res;
    res = 0.;
    double tmp;

    for (q1 = 0; q1 < L; ++q1)
        for (q2 = 0; q2 < L; ++q2)
            for (q3 = 0; q3 < L; ++q3)
                for (q41 = 0; q41 < LT; ++q41)
                    for (q42 = 0; q42 < LT; ++q42)
                    {
#ifdef BC_T_PERIODIC
                        mom1 = (fourvec){{q1, q2, q3, ((double)q41) * L / LT}};
#elif BC_T_ANTIPERIODIC
                        mom1 = (fourvec){{q1, q2, q3, ((double)2 * q41 + 1) * L / (2 * LT)}};
#endif
                        iadd(&mom1, &p);
                        imul(&mom1, 2.0 * PI / L);
#ifdef BC_T_PERIODIC
                        mom2 = (fourvec){{q1, q2, q3, ((double)q42) * L / LT}};
#elif BC_T_ANTIPERIODIC
                        mom2 = (fourvec){{q1, q2, q3, ((double)2 * q42 + 1) * L / (2 * LT)}};
#endif
                        imul(&mom2, 2.0 * PI / L);

                        tmp = (f1(mom1, m) * f1(mom2, m) + f2(mom1, mom2)) / ((SQR(f1(mom1, m)) + f2(mom1, mom1)) * (SQR(f1(mom2, m)) + f2(mom2, mom2)));
                        res += tmp * cexp(I * (2.0 * PI / LT) * t * (q42 - q41));
                    }

    res = 4 * res / L / L / L / LT / LT;
    return res;
}

/**
 * @brief Calculates analytic expression for rho 2-point function with gamma3 at source and sink
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
double complex twopoint_rho(fourvec p, double m, int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    double complex res;
    res = 0.;
    double tmp;

    for (q1 = 0; q1 < L; ++q1)
        for (q2 = 0; q2 < L; ++q2)
            for (q3 = 0; q3 < L; ++q3)
                for (q41 = 0; q41 < LT; ++q41)
                    for (q42 = 0; q42 < LT; ++q42)
                    {
#ifdef BC_T_PERIODIC
                        mom1 = (fourvec){{q1, q2, q3, ((double)q41) * L / LT}};
#elif BC_T_ANTIPERIODIC
                        mom1 = (fourvec){{q1, q2, q3, ((double)2 * q41 + 1) * L / (2 * LT)}};
#endif
                        iadd(&mom1, &p);
                        imul(&mom1, 2.0 * PI / L);
#ifdef BC_T_PERIODIC
                        mom2 = (fourvec){{q1, q2, q3, ((double)q42) * L / LT}};
#elif BC_T_ANTIPERIODIC
                        mom2 = (fourvec){{q1, q2, q3, ((double)2 * q42 + 1) * L / (2 * LT)}};
#endif

                        imul(&mom2, 2.0 * PI / L);

                        tmp = (f1(mom1, m) * f1(mom2, m) + f2(mom1, mom2) - 2 * sin(mom1.v[2]) * sin(mom2.v[2])) / ((SQR(f1(mom1, m)) + f2(mom1, mom1)) * (SQR(f1(mom2, m)) + f2(mom2, mom2)));
                        res += tmp * cexp(I * (2.0 * PI / LT) * t * (q42 - q41));
                    }

    res = 4 * res / L / L / L / LT / LT;
    return res;
}

/**
 * @brief Calculates analytic expression for rho 2-point function with gamma1 at source and gamma2 at the sink
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
double complex twopoint_rho12(fourvec p, double m, int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    double complex res;
    res = 0.;
    double tmp;

    for (q1 = 0; q1 < L; ++q1)
        for (q2 = 0; q2 < L; ++q2)
            for (q3 = 0; q3 < L; ++q3)
                for (q41 = 0; q41 < LT; ++q41)
                    for (q42 = 0; q42 < LT; ++q42)
                    {
#ifdef BC_T_PERIODIC
                        mom1 = (fourvec){{q1, q2, q3, ((double)q41) * L / LT}};
#elif BC_T_ANTIPERIODIC
                        mom1 = (fourvec){{q1, q2, q3, ((double)2 * q41 + 1) * L / (2 * LT)}};
#endif
                        iadd(&mom1, &p);
                        imul(&mom1, 2.0 * PI / L);
#ifdef BC_T_PERIODIC
                        mom2 = (fourvec){{q1, q2, q3, ((double)q42) * L / LT}};
#elif BC_T_ANTIPERIODIC
                        mom2 = (fourvec){{q1, q2, q3, ((double)2 * q42 + 1) * L / (2 * LT)}};
#endif
                        imul(&mom2, 2.0 * PI / L);

                        tmp = (-(sin(mom1.v[0]) * sin(mom2.v[1]) + sin(mom1.v[1]) * sin(mom2.v[0]))) / ((SQR(f1(mom1, m)) + f2(mom1, mom1)) * (SQR(f1(mom2, m)) + f2(mom2, mom2)));
                        res += tmp * cexp(I * (2.0 * PI / LT) * t * (q42 - q41));
                    }

    res = 4 * res / L / L / L / LT / LT;
    return res;
}

#define Q(A, L) (2 * PI * (A) / (L))
#define FV(A, B)                                         \
    (fourvec)                                            \
    {                                                    \
        {                                                \
            Q(A##1, L), Q(A##2, L), Q(A##3, L), Q(B, LT) \
        }                                                \
    }

/**
 * @brief Calculates analytic expression for the rectangle pipi->pipi contraction.
 * @param px Momentum at the sink connected to the source
 * @param py Momentum at the source
 * @param pz Momentum at the sink on the vertex opposite to the source
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
double complex Rect(fourvec px, fourvec py, fourvec pz, double m, int L, int LT, int t)
{
    fourvec mom[4];
    int q11, q12, q13, q14, q24, q34, q44, i, j;
    double complex res;
    res = 0.;
    double numerator, denominator;
    double af1[4];
    double af2[4][4];

    for (q11 = 0; q11 < L; ++q11)
        for (q12 = 0; q12 < L; ++q12)
            for (q13 = 0; q13 < L; ++q13)
                for (q14 = 0; q14 < LT; ++q14)
                    for (q24 = 0; q24 < LT; ++q24)
                        for (q34 = 0; q34 < LT; ++q34)
                            for (q44 = 0; q44 < LT; ++q44)
                            {
#ifdef BC_T_PERIODIC
                                mom[0] = FV(q1, q14);
                                mom[1] = (fourvec){{q11, q12, q13, ((double)q24) * L / LT}};
                                iadd(&mom[1], &px);
                                imul(&mom[1], 2.0 * PI / L);
                                mom[2] = (fourvec){{q11, q12, q13, ((double)q34) * L / LT}};
                                iadd(&mom[2], &px);
                                iadd(&mom[2], &pz);
                                imul(&mom[2], 2.0 * PI / L);
                                mom[3] = (fourvec){{q11, q12, q13, ((double)q44) * L / LT}};
                                iadd(&mom[3], &px);
                                iadd(&mom[3], &py);
                                iadd(&mom[3], &pz);
                                imul(&mom[3], 2.0 * PI / L);
#elif BC_T_ANTIPERIODIC
                                mom[0] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q14 + 1) * PI / (LT)}};
                                mom[1] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q24 + 1) * PI / (LT)}};
                                mom[2] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q34 + 1) * PI / (LT)}};
                                mom[3] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q44 + 1) * PI / (LT)}};
#endif
                                denominator = 1.0;
                                for (i = 0; i < 4; ++i)
                                {
                                    af1[i] = f1(mom[i], m);
                                    for (j = 0; j < 4; ++j)
                                    {
                                        af2[i][j] = f2(mom[i], mom[j]);
                                    }
                                    denominator *= (SQR(af1[i]) + af2[i][i]);
                                }

                                numerator = af1[0] * af1[1] * af1[2] * af1[3] + af1[0] * af1[1] * af2[2][3] - af1[0] * af1[2] * af2[1][3] + af1[0] * af1[3] * af2[1][2] + af1[1] * af1[2] * af2[0][3] - af1[1] * af1[3] * af2[0][2] + af1[2] * af1[3] * af2[0][1] + af2[0][1] * af2[2][3] - af2[0][2] * af2[1][3] + af2[0][3] * af2[1][2];

                                res += cexp(I * (t * (q24 - q44)) * 2.0 * PI / LT) * numerator / denominator;
                            }

    res = 4 * res / L / L / L / LT / LT / LT / LT;
    return res;
}

complex double C(fourvec px, fourvec py, double m, int L, int LT, int t)
{
    fourvec mom[4];
    int q11, q12, q13, q14, q24, q34, q44, i, j;
    complex double res;
    res = 0.0 * I;
    double numerator, denominator;
    double af1[4];
    double af2[4][4];

    for (q11 = 0; q11 < L; ++q11)
        for (q12 = 0; q12 < L; ++q12)
            for (q13 = 0; q13 < L; ++q13)
                for (q14 = 0; q14 < LT; ++q14)
                    for (q24 = 0; q24 < LT; ++q24)
                        for (q34 = 0; q34 < LT; ++q34)
                            for (q44 = 0; q44 < LT; ++q44)
                            {
#ifdef BC_T_PERIODIC
                                mom[0] = FV(q1, q14);
                                mom[1] = (fourvec){{q11, q12, q13, ((double)q24) * L / LT}};
                                iadd(&mom[1], &px);
                                imul(&mom[1], 2.0 * PI / L);
                                mom[2] = (fourvec){{q11, q12, q13, ((double)q34) * L / LT}};
                                iadd(&mom[2], &px);
                                imul(&mom[2], 2.0 * PI / L);
                                mom[3] = (fourvec){{q11, q12, q13, ((double)q44) * L / LT}};
                                iadd(&mom[3], &px);
                                iadd(&mom[3], &py);
                                imul(&mom[3], 2.0 * PI / L);
#elif BC_T_ANTIPERIODIC
                                mom[0] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q14 + 1) * PI / (LT)}};
                                mom[1] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q24 + 1) * PI / (LT)}};
                                mom[2] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q34 + 1) * PI / (LT)}};
                                mom[3] = (fourvec){{q11 * 2.0 * PI / L, q12 * 2.0 * PI / L, q13 * 2.0 * PI / L, ((double)2 * q44 + 1) * PI / (LT)}};
#endif
                                denominator = 1.0;
                                for (i = 0; i < 4; ++i)
                                {
                                    af1[i] = f1(mom[i], m);
                                    for (j = 0; j < 4; ++j)
                                    {
                                        af2[i][j] = f2(mom[i], mom[j]);
                                    }
                                    denominator *= (SQR(af1[i]) + af2[i][i]);
                                }

                                numerator = af1[0] * af1[1] * af1[2] * af1[3] + af1[0] * af1[1] * af2[2][3] - af1[0] * af1[2] * af2[1][3] + af1[0] * af1[3] * af2[1][2] + af1[1] * af1[2] * af2[0][3] - af1[1] * af1[3] * af2[0][2] + af1[2] * af1[3] * af2[0][1] + af2[0][1] * af2[2][3] - af2[0][2] * af2[1][3] + af2[0][3] * af2[1][2];

                                res += cexp((double)(t * (-q14 + q24 - q34 + q44)) * 2.0 * I * PI / LT) * numerator / denominator;
                            }

    return 4 * res / L / L / L / L / L / L / LT / LT / LT / LT;
}
int compare_corr_disc(double complex *corr_num, char *name, double tol)
{
    int retval = 0;
    for (int t = 0; t < GLB_T; t++)
    {
        if (cabs(corr_num[t]) > tol)
        {
            lprintf("TEST", 0, "Mismatch %s, t=%d,  numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name, t, creal(corr_num[t]), cimag(corr_num[t]), 0., 0.);
            retval += 1;
        }
        else
        {
            lprintf("TEST", 0, "Match %s, t=%d, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name, t, creal(corr_num[t]), cimag(corr_num[t]), 0., 0.);
        }
    }

    return retval;
}

int compare_corr(double complex *corr_ex, double complex *corr_num, int tstart, char *name, double tol)
{
    int retval = 0;
    for (int t = tstart; t < GLB_T; t++)
    {
        if (cabs(corr_ex[t] - corr_num[t]) / cabs(corr_ex[t]) > tol)
        {
            lprintf("TEST", 0, "Mismatch %s, t=%d, relative diff: %e, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name, t, cabs(corr_ex[t] - corr_num[t]) / cabs(corr_ex[t]), creal(corr_num[t]), cimag(corr_num[t]), creal(corr_ex[t]), cimag(corr_ex[t]));
            retval += 1;
        }
        else
        {
            lprintf("TEST", 0, "Match %s, t=%d, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name, t, creal(corr_num[t]), cimag(corr_num[t]), creal(corr_ex[t]), cimag(corr_ex[t]));
        }
    }

    return retval;
}

int main(int argc, char *argv[])
{
    int return_value = 0;
    double m[256];
    int ncorr = 16;
    double tol = 1.5e-1;
    double tol_disc = 1.e-10;
    meson_observable **mo_arr;
    fourvec zero_p = (fourvec){{0, 0, 0, 0}};

    error(!(GLB_X == GLB_Y && GLB_X == GLB_Z), 1, "main", "This test works only for GLB_X=GLB_Y=GLB_Z");

    setup_process(&argc, &argv);

    setup_gauge_fields();

    read_input(glb_var.read, get_input_filename());
    read_input(mes_var.read, get_input_filename());
    read_input(rlx_var.read, get_input_filename());

    int numsources = mes_var.nhits;

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&mes_var.csw);
#endif
    m[0] = atof(mes_var.mstring);
    init_propagator_eo(1, m, mes_var.precision);

    lprintf("MAIN", 0, "mass is : %e\n", m[0]);
    lprintf("MAIN", 0, "Number of hits : %d \n", numsources);

    struct timeval start, end, etime;
    gettimeofday(&start, 0);

    unit_gauge(u_gauge);
    represent_gauge_field();
#ifdef REPR_FUNDAMENTAL
    apply_BCs_on_represented_gauge_field(); //This is a trick: the BCs are not applied in the case the REPR is fundamental because represent_gauge field assumes that the right BCs are already applied on the fundamental field!
#endif

    double complex Pistoch[GLB_T], Pitheo[GLB_T];
    double complex Rhostoch[GLB_T], Rhotheo[GLB_T];
    double complex Rstoch[GLB_T], Rtheo[GLB_T];
    double complex Raltstoch[GLB_T];
    double complex Dstoch[GLB_T], Dtheo[GLB_T];
    double complex Cstoch[GLB_T], Ctheo[GLB_T];
    double complex Vstoch[GLB_T], Vtheo[GLB_T];
    double complex Discstoch[GLB_T];

    for (int t = 0; t < GLB_T; t++)
    {
        Rstoch[t] = 0.0;
        Raltstoch[t] = 0.0;
        Dstoch[t] = 0.0;
        Cstoch[t] = 0.0;
        Vstoch[t] = 0.0;
        Pistoch[t] = 0.0;
        Rhostoch[t] = 0.0;
        Discstoch[t] = 0.0;
    }

    mo_arr = (meson_observable **)malloc(sizeof(meson_observable *) * ncorr);
    for (int i = 0; i < ncorr; i++)
        mo_arr[i] = (meson_observable *)malloc(sizeof(meson_observable));
    init_mo(mo_arr[0], "D", GLB_T);
    init_mo(mo_arr[1], "C", GLB_T);
    init_mo(mo_arr[2], "R", GLB_T);
    init_mo(mo_arr[3], "V", GLB_T);
    init_mo(mo_arr[4], "pi", GLB_T);
    init_mo(mo_arr[5], "rho11", GLB_T);
    init_mo(mo_arr[6], "rho12", GLB_T);
    init_mo(mo_arr[7], "rho13", GLB_T);
    init_mo(mo_arr[8], "rho21", GLB_T);
    init_mo(mo_arr[9], "rho22", GLB_T);
    init_mo(mo_arr[10], "rho23", GLB_T);
    init_mo(mo_arr[11], "rho31", GLB_T);
    init_mo(mo_arr[12], "rho32", GLB_T);
    init_mo(mo_arr[13], "rho33", GLB_T);
    init_mo(mo_arr[14], "Ralt", GLB_T);
    init_mo(mo_arr[15], "disc", GLB_T);

    for (int n = 0; n < numsources; n++)
    {
        if (n % 100 == 0)
            lprintf("MAIN", 0, "nhits: %d / %d \n/", n, numsources);
        measure_pion_scattering_I0(m, 1, mes_var.precision, NULL, NULL, 1, mo_arr);

        // now copy the content to a placeholder.
        for (int t = 0; t < GLB_T; t++)
        {
            Dstoch[t] += (mo_arr[0]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[0]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (numsources);
            Cstoch[t] += (mo_arr[1]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[1]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (numsources);
            Rstoch[t] += (mo_arr[2]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[2]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (numsources);
            Vstoch[t] += (mo_arr[3]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[3]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (numsources);
            Pistoch[t] += (mo_arr[4]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[4]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (numsources);
            Rhostoch[t] += (mo_arr[5]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[5]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (3 * numsources);
            Rhostoch[t] -= (mo_arr[9]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[9]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (3 * numsources); // there is an i in the interpolating field
            Rhostoch[t] += (mo_arr[13]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[13]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (3 * numsources);
            Discstoch[t] += (mo_arr[15]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[15]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]) / (numsources);
        }
        for (int i = 0; i < ncorr; i++)
            reset_mo(mo_arr[i]);

        //io4pt(R,0,n,path,"R2","test"); // to write the output to a file.
    }

    for (int t = 0; t < GLB_T; t++)
    {
        Dtheo[t] = cpow(NF * GLB_X * GLB_Y * GLB_Z * twopoint(zero_p, m[0], GLB_X, GLB_T, t), 2);
        Ctheo[t] = NF * cpow(GLB_X * GLB_Y * GLB_Z, 2) * C(zero_p, zero_p, m[0], GLB_X, GLB_T, t);
        Vtheo[t] = NF * GLB_X * GLB_Y * GLB_Z * twopoint(zero_p, m[0], GLB_X, GLB_T, 0);
        Pitheo[t] = NF * GLB_X * GLB_Y * GLB_Z * twopoint(zero_p, m[0], GLB_X, GLB_T, t);
        Rhotheo[t] = NF * GLB_X * GLB_Y * GLB_Z * twopoint_rho(zero_p, m[0], GLB_X, GLB_T, t); // VD why no NF ??
        Rtheo[t] = NF * GLB_X * GLB_Y * GLB_Z * Rect(zero_p, zero_p, zero_p, m[0], GLB_X, GLB_T, t);
    }

    //  for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"V analytical V %e + 1I %e  numerical %e + 1I %e  ratio (num/theo) %e + 1I %e\n" ,creal(Vtheo[i]),cimag(Vtheo[i]),creal(Vstoch[i]), cimag(Vstoch[i]),creal(Vstoch[i]/Vtheo[i]),cimag(Vstoch[i]/Vtheo[i]));
    //  for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"C analytical C %e + 1I %e  numerical %e + 1I %e  ratio (num/theo) %e + 1I %e\n" ,creal(Ctheo[i]),cimag(Ctheo[i]),creal(Cstoch[i]), cimag(Cstoch[i]),creal(Cstoch[i]/Ctheo[i]),cimag(Cstoch[i]/Ctheo[i]));
    //  for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"D analytical D %e + 1I %e  numerical %e + 1I %e  ratio (num/theo) %e + 1I %e\n" ,creal(Dtheo[i]),cimag(Dtheo[i]),creal(Dstoch[i]), cimag(Dstoch[i]),creal(Dstoch[i]/Dtheo[i]),cimag(Dstoch[i]/Dtheo[i]));
    //  for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"R analytical R %e + 1I %e  numerical %e + 1I %e  ratio (num/theo) %e + 1I %e\n" ,creal(Rtheo[i]),cimag(Rtheo[i]),creal(Rstoch[i]), cimag(Rstoch[i]),creal(Rstoch[i]/Rtheo[i]),cimag(Rstoch[i]/Rtheo[i]));

    // Now compare R and Ralt
    measure_pion_scattering_I0(m, numsources, mes_var.precision, NULL, NULL, 0, mo_arr); // now compute Ralt

    for (int t = 0; t < GLB_T; t++)
        Raltstoch[t] = (mo_arr[14]->corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo_arr[14]->corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]);
    for (int i = 0; i < GLB_T; ++i)
        lprintf("MAIN", 0, "R analytical Ralt %e + 1I %e  numerical %e + 1I %e  ratio (num/theo) %e + 1I %e\n", creal(Rtheo[i]), cimag(Rtheo[i]), creal(Raltstoch[i]), cimag(Raltstoch[i]), creal(Raltstoch[i] / Rtheo[i]), cimag(Raltstoch[i] / Rtheo[i]));

    return_value += compare_corr(Pitheo, Pistoch, 0, "Pi", tol);
    return_value += compare_corr(Rhotheo, Rhostoch, 0, "Rho (averaged over 11,22,33) ", tol);
    return_value += compare_corr(Vtheo, Vstoch, 0, "V", tol);
    return_value += compare_corr(Ctheo, Cstoch, 0, "C", tol);
    return_value += compare_corr(Dtheo, Dstoch, 0, "D", tol);
    return_value += compare_corr(Rtheo, Rstoch, 0, "R", tol);
    return_value += compare_corr(Rtheo, Raltstoch, 1, "Ralt", tol); // t=0 is incorrect for Ralt.
    return_value += compare_corr_disc(Discstoch, "disc", tol_disc);

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Configuration : analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

    global_sum_int(&return_value, 1);
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);

    for (int i = 0; i < ncorr; i++)
        free_mo(mo_arr[i]);

    finalize_process();

    return return_value;
}
