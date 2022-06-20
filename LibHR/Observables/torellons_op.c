
#include <stdlib.h>
#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "utils.h"
#include "glueballs.h"
#include <string.h>
#define ntors 29
static double PI = 3.141592653589793238462643383279502884197;
static double complex *tor_path_storage = NULL;
static double complex poly0(int in)
{
    return polyleg(in, 1)->tr;
}

static double complex poly1(int in)
{
    return polyleg(in, 2)->tr;
}

static double complex poly2(int in)
{
    return polyleg(in, 3)->tr;
}

static inline double diPoly_p_0_0_0_Ir_1_C_1_n_1(int idx)
{
    return +(16.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[0 + idx]) + (16.) * cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[1 + idx]) + (16.) * cimag(tor_path_storage[2 + idx]) * cimag(tor_path_storage[2 + idx]) + (16.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[0 + idx]) + (16.) * creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[1 + idx]) + (16.) * creal(tor_path_storage[2 + idx]) * creal(tor_path_storage[2 + idx]);
}

static double complex poly8(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly10(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly6(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly13(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly11(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly12(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly7(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly9(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static inline double diPoly_p_0_0_1_Ir_1_C_1_n_1(int idx)
{
    return +cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[8 + idx]) + cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[10 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[6 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[13 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[11 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[12 + idx]) + cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[7 + idx]) + cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[9 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[8 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[10 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[6 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[13 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[11 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[12 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[7 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[9 + idx]);
}

static double complex poly4(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly5(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static inline double diPoly_p_0_0_1_Ir_3_C_1_n_1(int idx)
{
    return +(-2.8284271247461901) * cimag(tor_path_storage[4 + idx]) * creal(tor_path_storage[1 + idx]) + (-2.8284271247461901) * cimag(tor_path_storage[5 + idx]) * creal(tor_path_storage[1 + idx]) + (2.8284271247461901) * cimag(tor_path_storage[1 + idx]) * creal(tor_path_storage[4 + idx]) + (2.8284271247461901) * cimag(tor_path_storage[1 + idx]) * creal(tor_path_storage[5 + idx]);
}

static double complex poly14(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static double complex poly15(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    wl = polyleg(in, 3);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static inline double diPoly_p_0_1_0_Ir_1_C_1_n_1(int idx)
{
    return +(4.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[14 + idx]) + (4.) * cimag(tor_path_storage[2 + idx]) * cimag(tor_path_storage[15 + idx]) + (4.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[14 + idx]) + (4.) * creal(tor_path_storage[2 + idx]) * creal(tor_path_storage[15 + idx]);
}

static double complex poly3(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly16(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly17(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly18(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static inline double diPoly_p_1_0_0_Ir_1_C_1_n_1(int idx)
{
    return +(2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[3 + idx]) + (2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[16 + idx]) + (2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[17 + idx]) + (2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[18 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[3 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[16 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[17 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[18 + idx]);
}

static double complex poly19(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static double complex poly21(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static double complex poly22(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static inline double diPoly_p_1_0_0_Ir_1_C_1_n_2(int idx)
{
    return +(2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[14 + idx]) + (2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[19 + idx]) + (2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[21 + idx]) + (2.) * cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[22 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[14 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[19 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[21 + idx]) + (2.) * creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[22 + idx]);
}

static double complex poly20(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    wl = polyleg(in, 3);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static inline double diPoly_p_1_0_1_Ir_1_C_1_n_1(int idx)
{
    return +cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[14 + idx]) + cimag(tor_path_storage[2 + idx]) * cimag(tor_path_storage[15 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[19 + idx]) + cimag(tor_path_storage[2 + idx]) * cimag(tor_path_storage[20 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[14 + idx]) + creal(tor_path_storage[2 + idx]) * creal(tor_path_storage[15 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[19 + idx]) + creal(tor_path_storage[2 + idx]) * creal(tor_path_storage[20 + idx]);
}

static double complex poly23(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static double complex poly24(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res1, res, wl->p[1]);
    _suNg_trace(p, res1);
    return p;
}

static inline double diPoly_p_1_1_0_Ir_1_C_1_n_1(int idx)
{
    return +cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[23 + idx]) + cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[8 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[24 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[6 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[23 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[8 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[24 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[6 + idx]);
}

static double complex poly26(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static double complex poly27(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    wl = polyleg(in, 2);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static double complex poly28(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static double complex poly25(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;
    wilson_lines *wl;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    wl = polyleg(in, 1);

    _suNg_times_suNg(res, res1, wl->p[0]);
    _suNg_trace(p, res);
    return p;
}

static inline double diPoly_p_1_1_0_Ir_1_C_1_n_2(int idx)
{
    return +cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[26 + idx]) + cimag(tor_path_storage[1 + idx]) * cimag(tor_path_storage[27 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[28 + idx]) + cimag(tor_path_storage[0 + idx]) * cimag(tor_path_storage[25 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[26 + idx]) + creal(tor_path_storage[1 + idx]) * creal(tor_path_storage[27 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[28 + idx]) + creal(tor_path_storage[0 + idx]) * creal(tor_path_storage[25 + idx]);
}

static int last_t = -10;
void request_space_tors_evaluation() { last_t = -10; }
static void eval_time_momentum_torellons(int t, int px, int py, int pz, double complex *np)
{
    int nnx, nny, nnz, idx = 0, in;
    double complex ce = I * 2.0 * PI / GLB_X;
    if (tor_path_storage == NULL)
    {
        tor_path_storage = malloc(ntors * X * Y * Z * sizeof(double complex));
        for (in = 0; in < ntors * X * Y * Z; in++)
            tor_path_storage[in] = 0.;
    };
    if (t != last_t)
    {
        last_t = t;
        for (nny = 0; nny < Y; nny++)
            for (nnz = 0; nnz < Z; nnz++)
                for (nnx = 0; nnx < X; nnx++)
                {
                    in = ipt(t, nnx, nny, nnz);
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    tor_path_storage[0 + idx] = poly0(in);
                    tor_path_storage[3 + idx] = poly3(in);
                    tor_path_storage[6 + idx] = poly6(in);
                    tor_path_storage[11 + idx] = poly11(in);
                    tor_path_storage[12 + idx] = poly12(in);
                    tor_path_storage[13 + idx] = poly13(in);
                    tor_path_storage[14 + idx] = poly14(in);
                    tor_path_storage[16 + idx] = poly16(in);
                    tor_path_storage[17 + idx] = poly17(in);
                    tor_path_storage[18 + idx] = poly18(in);
                    tor_path_storage[19 + idx] = poly19(in);
                    tor_path_storage[21 + idx] = poly21(in);
                    tor_path_storage[22 + idx] = poly22(in);
                    tor_path_storage[24 + idx] = poly24(in);
                    tor_path_storage[25 + idx] = poly25(in);
                    tor_path_storage[28 + idx] = poly28(in);
                };
        for (nnz = 0; nnz < Z; nnz++)
            for (nnx = 0; nnx < X; nnx++)
                for (nny = 0; nny < Y; nny++)
                {
                    in = ipt(t, nnx, nny, nnz);
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    tor_path_storage[1 + idx] = poly1(in);
                    tor_path_storage[4 + idx] = poly4(in);
                    tor_path_storage[5 + idx] = poly5(in);
                    tor_path_storage[7 + idx] = poly7(in);
                    tor_path_storage[8 + idx] = poly8(in);
                    tor_path_storage[9 + idx] = poly9(in);
                    tor_path_storage[10 + idx] = poly10(in);
                    tor_path_storage[23 + idx] = poly23(in);
                    tor_path_storage[26 + idx] = poly26(in);
                    tor_path_storage[27 + idx] = poly27(in);
                };
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    in = ipt(t, nnx, nny, nnz);
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    tor_path_storage[2 + idx] = poly2(in);
                    tor_path_storage[15 + idx] = poly15(in);
                    tor_path_storage[20 + idx] = poly20(in);
                };
    };
    if (px == 0 && py == 0 && pz == 0)
    {
        np[0] = 0.0;
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    np[0] += diPoly_p_0_0_0_Ir_1_C_1_n_1(idx);
                };
    };
    if (px == 0 && py == 0 && pz == 1)
    {
        np[1] = 0.0;
        np[2] = 0.0;
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    np[1] += cexp(ce * (double)(nnz)) * diPoly_p_0_0_1_Ir_1_C_1_n_1(idx);
                    np[2] += cexp(ce * (double)(nnz)) * diPoly_p_0_0_1_Ir_3_C_1_n_1(idx);
                };
    };
    if (px == 0 && py == 1 && pz == 0)
    {
        np[3] = 0.0;
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    np[3] += cexp(ce * (double)(nny)) * diPoly_p_0_1_0_Ir_1_C_1_n_1(idx);
                };
    };
    if (px == 1 && py == 0 && pz == 0)
    {
        np[4] = 0.0;
        np[5] = 0.0;
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    np[4] += cexp(ce * (double)(nnx)) * diPoly_p_1_0_0_Ir_1_C_1_n_1(idx);
                    np[5] += cexp(ce * (double)(nnx)) * diPoly_p_1_0_0_Ir_1_C_1_n_2(idx);
                };
    };
    if (px == 1 && py == 0 && pz == 1)
    {
        np[6] = 0.0;
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    np[6] += cexp(ce * (double)(nnx + nnz)) * diPoly_p_1_0_1_Ir_1_C_1_n_1(idx);
                };
    };
    if (px == 1 && py == 1 && pz == 0)
    {
        np[7] = 0.0;
        np[8] = 0.0;
        for (nnx = 0; nnx < X; nnx++)
            for (nny = 0; nny < Y; nny++)
                for (nnz = 0; nnz < Z; nnz++)
                {
                    idx = ntors * (nnx + X * (nny + Y * nnz));
                    np[7] += cexp(ce * (double)(nnx + nny)) * diPoly_p_1_1_0_Ir_1_C_1_n_1(idx);
                    np[8] += cexp(ce * (double)(nnx + nny)) * diPoly_p_1_1_0_Ir_1_C_1_n_2(idx);
                };
    };
};
void eval_all_torellon_ops(int t, double complex *numerical_tor_out)
{
    static double complex *numerical_op = NULL;
    if (numerical_op == NULL)
    {
        numerical_op = malloc(total_n_tor_op * sizeof(double complex));
    }
    request_space_tors_evaluation();
    eval_time_momentum_torellons(t, 0, 0, 0, numerical_op);
    eval_time_momentum_torellons(t, 0, 0, 1, numerical_op);
    eval_time_momentum_torellons(t, 0, 0, 1, numerical_op);
    eval_time_momentum_torellons(t, 0, 1, 0, numerical_op);
    eval_time_momentum_torellons(t, 1, 0, 0, numerical_op);
    eval_time_momentum_torellons(t, 1, 0, 1, numerical_op);
    eval_time_momentum_torellons(t, 1, 1, 0, numerical_op);
    for (int i = 0; i < total_n_tor_op; i++)
        numerical_tor_out[i] += numerical_op[i];
}

void collect_1pt_torellon_functions(cor_list *lcor, double complex *tor_storage)
{
    int n1, i;
    static double complex *tor1_bf;
    static int n_total_active_slices = 0;
    static int *listactive = NULL;
    if (listactive == NULL)
    {
        listactive = malloc(sizeof(int) * GLB_T);
        for (i = 0; i < GLB_T; i++)
            listactive[i] = -1;

        for (i = 0; i < lcor->n_entries; i++)
        {
            listactive[lcor->list[i].t1] = 1;
            listactive[lcor->list[i].t2] = 1;
        }

        for (i = 0; i < GLB_T; i++)
            if (listactive[i] == 1)
            {
                listactive[i] = n_total_active_slices;
                n_total_active_slices++;
            }
    }

#ifdef WITH_MPI
    static int *t_to_proc = NULL;
    static int *listsent = NULL;
    int t1, t2;

    if (t_to_proc == NULL)
    {
        listsent = malloc(sizeof(int) * GLB_T);

        tor1_bf = malloc(sizeof(double complex) * total_n_tor_op * n_total_active_slices);

        t_to_proc = malloc(sizeof(int) * GLB_T);
        for (i = 0; i < GLB_T; i++)
        {
            int x[4] = {i, 0, 0, 0}, p[4];
            glb_to_proc(x, p);
            t_to_proc[i] = p[0];
        }
    }

    for (i = 0; i < GLB_T; i++)
        listsent[i] = -1;

    static double complex *tor2;
    static double complex *tor1;
    MPI_Request req_1pt[GLB_T];

    for (int icor = 0; icor < lcor->n_entries; icor++)
    {

        t1 = glbT_to_active_slices[lcor->list[icor].t1];
        t2 = glbT_to_active_slices[lcor->list[icor].t2];
        tor1 = tor1_bf + total_n_tor_op * listactive[lcor->list[icor].t1];
        tor2 = tor1_bf + total_n_tor_op * listactive[lcor->list[icor].t2];

        if (listsent[lcor->list[icor].t1] == -1)
        {
            if (t1 != -1)
            {
                listsent[lcor->list[icor].t1] = 0;
                if (PID == 0)
                {
                    memcpy(tor1, tor_storage + t1 * total_n_tor_op, sizeof(double complex) * total_n_tor_op);
                }
                else
                {
                    MPI_Isend((double *)(tor_storage + t1 * total_n_tor_op), total_n_tor_op * 2, MPI_DOUBLE, 0, lcor->list[icor].t1, cart_comm, req_1pt + lcor->list[icor].t1);
                }
            }

            if (PID == 0 && t1 == -1)
            {
                MPI_Irecv((double *)(tor1), total_n_tor_op * 2, MPI_DOUBLE, t_to_proc[lcor->list[icor].t1], lcor->list[icor].t1, cart_comm, req_1pt + lcor->list[icor].t1);
                listsent[lcor->list[icor].t1] = 1;
            }
        }

        if (listsent[lcor->list[icor].t2] == -1)
        {
            listsent[lcor->list[icor].t2] = 0;
            if (t2 != -1)
            {
                if (PID == 0)
                {
                    memcpy(tor2, tor_storage + t2 * total_n_tor_op, sizeof(double complex) * total_n_tor_op);
                }
                else
                {
                    MPI_Isend((double *)(tor_storage + t2 * total_n_tor_op), total_n_tor_op * 2, MPI_DOUBLE, 0, lcor->list[icor].t2, cart_comm, req_1pt + lcor->list[icor].t2);
                }
            }

            if (PID == 0 && t2 == -1)
            {
                MPI_Irecv((double *)(tor2), total_n_tor_op * 2, MPI_DOUBLE, t_to_proc[lcor->list[icor].t2], lcor->list[icor].t2, cart_comm, req_1pt + lcor->list[icor].t2);
                listsent[lcor->list[icor].t2] = 1;
            }
        }
    }
    if (PID == 0)
        for (i = 0; i < GLB_T; i++)
            if (listsent[i] == 1)
            {
                MPI_Wait(req_1pt + i, MPI_STATUS_IGNORE);
            }

#else
    tor1_bf = tor_storage;
#endif

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=A1plusOhP Irrep ev=1/1 Charge=+ nop=%d\n", 1);
    lprintf("Measure ML", 0, "Tor id= 0 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 0; i < 1; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,1) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n", 1);
    lprintf("Measure ML", 0, "Tor id= 1 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 1; i < 2; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,1) Irrep=E2Dic4 Irrep ev=1/2 Charge=+ nop=%d\n", 1);
    lprintf("Measure ML", 0, "Tor id= 2 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 2; i < 3; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,1,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n", 1);
    lprintf("Measure ML", 0, "Tor id= 3 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 3; i < 4; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,0,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n", 2);
    lprintf("Measure ML", 0, "Tor id= 4 5 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 4; i < 6; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,0,1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 1);
    lprintf("Measure ML", 0, "Tor id= 6 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 6; i < 7; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,1,0) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 2);
    lprintf("Measure ML", 0, "Tor id= 7 8 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (i = 7; i < 9; i++)
                lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op * listactive[n1]]),
                        cimag(tor1_bf[i + total_n_tor_op * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }
}
void report_tor_group_setup()
{
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(0,0,0) Irrep=A1plusOhP Charge=+");
    lprintf("INIT Measure ML", 0, " |0=|");
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(0,0,1) Irrep=A1Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |1=-xyxzy-z|");
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(0,0,1) Irrep=E2Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |2=-xyyx,na|");
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(0,1,0) Irrep=A1Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |3=-yxy|");
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(1,0,0) Irrep=A1Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |4=-yxxy|5=-yxy|");
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(1,0,1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |6=-yxy|");
    lprintf("INIT Measure ML", 0, "\n1pt_tor Irrep multiplets Total P=(1,1,0) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |7=-xyx-zyz|8=-x-zyxz|");
}
