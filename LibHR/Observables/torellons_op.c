
#include <stdlib.h>
#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "utils.h"
#include "glueballs.h"
#include <string.h>
 #define ntors 19
 static double PI=3.141592653589793238462643383279502884197;
 static double complex *mom_def_Cp_poly_paths=NULL;
 static double complex *mom_def_Cm_poly_paths=NULL;
 static double complex *path_storage=NULL;
 static double complex poly1(int in)
{
return polyleg(in,2)->tr;
}

static double complex poly4(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly2(int in)
{
return polyleg(in,3)->tr;
}

static double complex poly5(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly6(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res, res1, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly7(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res, res1, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly0(int in)
{
return polyleg(in,1)->tr;
}

static double complex poly3(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly8(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly9(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res, res1, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly10(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res, res1, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly11(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly12(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly13(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res, res1, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static double complex poly14(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res, res1, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res1,res,wl->p[1]);
_suNg_trace(p,res1);
return p;
}

static void diPoly_p_0_0_0_Ir_1_C_1_n_1(double complex * tor_out)
{
*tor_out =+(4.)*mom_def_Cm_poly_paths[1]*mom_def_Cm_poly_paths[4]+(4.)*mom_def_Cm_poly_paths[2]*mom_def_Cm_poly_paths[5]+(4.)*mom_def_Cm_poly_paths[1]*mom_def_Cm_poly_paths[6]+(4.)*mom_def_Cm_poly_paths[2]*mom_def_Cm_poly_paths[7]+(4.)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[3]+(4.)*mom_def_Cm_poly_paths[2]*mom_def_Cm_poly_paths[8]+(4.)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[9]+(4.)*mom_def_Cm_poly_paths[2]*mom_def_Cm_poly_paths[10]+(4.)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[11]+(4.)*mom_def_Cm_poly_paths[1]*mom_def_Cm_poly_paths[12]+(4.)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[13]+(4.)*mom_def_Cm_poly_paths[1]*mom_def_Cm_poly_paths[14]+(4.)*mom_def_Cp_poly_paths[1]*mom_def_Cp_poly_paths[4]+(4.)*mom_def_Cp_poly_paths[2]*mom_def_Cp_poly_paths[5]+(4.)*mom_def_Cp_poly_paths[1]*mom_def_Cp_poly_paths[6]+(4.)*mom_def_Cp_poly_paths[2]*mom_def_Cp_poly_paths[7]+(4.)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[3]+(4.)*mom_def_Cp_poly_paths[2]*mom_def_Cp_poly_paths[8]+(4.)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[9]+(4.)*mom_def_Cp_poly_paths[2]*mom_def_Cp_poly_paths[10]+(4.)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[11]+(4.)*mom_def_Cp_poly_paths[1]*mom_def_Cp_poly_paths[12]+(4.)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[13]+(4.)*mom_def_Cp_poly_paths[1]*mom_def_Cp_poly_paths[14];
}

static double complex poly15(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[0]);
_suNg_trace(p,res1);
return p;
}

static double complex poly16(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[0]);
_suNg_trace(p,res1);
return p;
}

static double complex poly17(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[0]);
_suNg_trace(p,res1);
return p;
}

static double complex poly18(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
double complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res1,res,wl->p[0]);
_suNg_trace(p,res1);
return p;
}

static void diPoly_p_1_0_0_Ir_1_C_m1_n_1(double complex * tor_out)
{
*tor_out =+(2.)*mom_def_Cm_poly_paths[15]*mom_def_Cp_poly_paths[0]+(2.)*mom_def_Cm_poly_paths[16]*mom_def_Cp_poly_paths[0]+(2.)*mom_def_Cm_poly_paths[17]*mom_def_Cp_poly_paths[0]+(2.)*mom_def_Cm_poly_paths[18]*mom_def_Cp_poly_paths[0]+(-2.)*mom_def_Cm_poly_paths[0]*mom_def_Cp_poly_paths[15]+(-2.)*mom_def_Cm_poly_paths[0]*mom_def_Cp_poly_paths[16]+(-2.)*mom_def_Cm_poly_paths[0]*mom_def_Cp_poly_paths[17]+(-2.)*mom_def_Cm_poly_paths[0]*mom_def_Cp_poly_paths[18];
}

static void diPoly_p_1_0_0_Ir_3_C_1_n_1(double complex * tor_out)
{
*tor_out =+(-2.8284271247461901)*mom_def_Cm_poly_paths[11]*mom_def_Cp_poly_paths[0]+(2.8284271247461901)*mom_def_Cm_poly_paths[13]*mom_def_Cp_poly_paths[0]+(2.8284271247461901)*mom_def_Cm_poly_paths[0]*mom_def_Cp_poly_paths[11]+(-2.8284271247461901)*mom_def_Cm_poly_paths[0]*mom_def_Cp_poly_paths[13];
}

static void diPoly_p_1_0_0_Ir_3_C_1_n_2(double complex * tor_out)
{
*tor_out =+(2.8284271247461901)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[3]+(-2.8284271247461901)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[9]+(2.8284271247461901)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[3]+(-2.8284271247461901)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[9];
}

static void diPoly_p_1_1_0_Ir_1_C_1_n_1(double complex * tor_out)
{
*tor_out =+(2.)*mom_def_Cm_poly_paths[1]*mom_def_Cm_poly_paths[4]+(2.)*mom_def_Cm_poly_paths[0]*mom_def_Cm_poly_paths[3]+(2.)*mom_def_Cp_poly_paths[1]*mom_def_Cp_poly_paths[4]+(2.)*mom_def_Cp_poly_paths[0]*mom_def_Cp_poly_paths[3];
}

static int last_t = -10;
void request_space_tors_evaluation(){last_t=-10;}
  static void eval_time_momentum_torellons(int t, int px, int py, int pz)
  {
    int n_x, n_y, n_z, idx, in;
    double complex ce = 0.;
    if(path_storage==NULL)
      {        path_storage = malloc(ntors * X * Y * Z * sizeof(double complex));
        mom_def_Cp_poly_paths = malloc(ntors * sizeof(double complex));
        mom_def_Cm_poly_paths = malloc(ntors * sizeof(double complex));
        for (in = 0; in < ntors * X * Y * Z; in++)
            path_storage[in] = 0.;
    }
for(in = 0; in < ntors; in++)
{
mom_def_Cp_poly_paths[in]=0.;
mom_def_Cm_poly_paths[in]=0.;
}
if (t != last_t)
{
last_t=t;
for (n_y = 0; n_y < Y; n_y++)
for (n_z = 0; n_z < Z; n_z++)
for (n_x = 0; n_x < X; n_x++)
{
in = ipt(t, n_x, n_y, n_z);
ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
idx = ntors * (n_x + X * (n_y + Y * n_z));
path_storage[0+idx]= poly0(in);
mom_def_Cp_poly_paths[0]+=ce*creal(path_storage[0+idx]);
mom_def_Cm_poly_paths[0]+=I*ce*cimag(path_storage[0+idx]);
path_storage[3+idx]= poly3(in);
mom_def_Cp_poly_paths[3]+=ce*creal(path_storage[3+idx]);
mom_def_Cm_poly_paths[3]+=I*ce*cimag(path_storage[3+idx]);
path_storage[9+idx]= poly9(in);
mom_def_Cp_poly_paths[9]+=ce*creal(path_storage[9+idx]);
mom_def_Cm_poly_paths[9]+=I*ce*cimag(path_storage[9+idx]);
path_storage[11+idx]= poly11(in);
mom_def_Cp_poly_paths[11]+=ce*creal(path_storage[11+idx]);
mom_def_Cm_poly_paths[11]+=I*ce*cimag(path_storage[11+idx]);
path_storage[13+idx]= poly13(in);
mom_def_Cp_poly_paths[13]+=ce*creal(path_storage[13+idx]);
mom_def_Cm_poly_paths[13]+=I*ce*cimag(path_storage[13+idx]);
path_storage[15+idx]= poly15(in);
mom_def_Cp_poly_paths[15]+=ce*creal(path_storage[15+idx]);
mom_def_Cm_poly_paths[15]+=I*ce*cimag(path_storage[15+idx]);
path_storage[16+idx]= poly16(in);
mom_def_Cp_poly_paths[16]+=ce*creal(path_storage[16+idx]);
mom_def_Cm_poly_paths[16]+=I*ce*cimag(path_storage[16+idx]);
path_storage[17+idx]= poly17(in);
mom_def_Cp_poly_paths[17]+=ce*creal(path_storage[17+idx]);
mom_def_Cm_poly_paths[17]+=I*ce*cimag(path_storage[17+idx]);
path_storage[18+idx]= poly18(in);
mom_def_Cp_poly_paths[18]+=ce*creal(path_storage[18+idx]);
mom_def_Cm_poly_paths[18]+=I*ce*cimag(path_storage[18+idx]);
};
for (n_z = 0; n_z < Z; n_z++)
for (n_x = 0; n_x < X; n_x++)
for (n_y = 0; n_y < Y; n_y++)
{
in = ipt(t, n_x, n_y, n_z);
ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
idx = ntors * (n_x + X * (n_y + Y * n_z));
path_storage[1+idx]= poly1(in);
mom_def_Cp_poly_paths[1]+=ce*creal(path_storage[1+idx]);
mom_def_Cm_poly_paths[1]+=I*ce*cimag(path_storage[1+idx]);
path_storage[4+idx]= poly4(in);
mom_def_Cp_poly_paths[4]+=ce*creal(path_storage[4+idx]);
mom_def_Cm_poly_paths[4]+=I*ce*cimag(path_storage[4+idx]);
path_storage[6+idx]= poly6(in);
mom_def_Cp_poly_paths[6]+=ce*creal(path_storage[6+idx]);
mom_def_Cm_poly_paths[6]+=I*ce*cimag(path_storage[6+idx]);
path_storage[12+idx]= poly12(in);
mom_def_Cp_poly_paths[12]+=ce*creal(path_storage[12+idx]);
mom_def_Cm_poly_paths[12]+=I*ce*cimag(path_storage[12+idx]);
path_storage[14+idx]= poly14(in);
mom_def_Cp_poly_paths[14]+=ce*creal(path_storage[14+idx]);
mom_def_Cm_poly_paths[14]+=I*ce*cimag(path_storage[14+idx]);
};
for (n_x = 0; n_x < X; n_x++)
for (n_y = 0; n_y < Y; n_y++)
for (n_z = 0; n_z < Z; n_z++)
{
in = ipt(t, n_x, n_y, n_z);
ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
idx = ntors * (n_x + X * (n_y + Y * n_z));
path_storage[2+idx]= poly2(in);
mom_def_Cp_poly_paths[2]+=ce*creal(path_storage[2+idx]);
mom_def_Cm_poly_paths[2]+=I*ce*cimag(path_storage[2+idx]);
path_storage[5+idx]= poly5(in);
mom_def_Cp_poly_paths[5]+=ce*creal(path_storage[5+idx]);
mom_def_Cm_poly_paths[5]+=I*ce*cimag(path_storage[5+idx]);
path_storage[7+idx]= poly7(in);
mom_def_Cp_poly_paths[7]+=ce*creal(path_storage[7+idx]);
mom_def_Cm_poly_paths[7]+=I*ce*cimag(path_storage[7+idx]);
path_storage[8+idx]= poly8(in);
mom_def_Cp_poly_paths[8]+=ce*creal(path_storage[8+idx]);
mom_def_Cm_poly_paths[8]+=I*ce*cimag(path_storage[8+idx]);
path_storage[10+idx]= poly10(in);
mom_def_Cp_poly_paths[10]+=ce*creal(path_storage[10+idx]);
mom_def_Cm_poly_paths[10]+=I*ce*cimag(path_storage[10+idx]);
};
}
else{
for (n_y = 0; n_y < Y; n_y++)
for (n_z = 0; n_z < Z; n_z++)
for (n_x = 0; n_x < X; n_x++)
{
in = ipt(t, n_x, n_y, n_z);
ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
idx = ntors * (n_x + X * (n_y + Y * n_z));
mom_def_Cp_poly_paths[0]+=ce*creal(path_storage[0+idx]);
mom_def_Cm_poly_paths[0]+=I*ce*cimag(path_storage[0+idx]);
mom_def_Cp_poly_paths[3]+=ce*creal(path_storage[3+idx]);
mom_def_Cm_poly_paths[3]+=I*ce*cimag(path_storage[3+idx]);
mom_def_Cp_poly_paths[9]+=ce*creal(path_storage[9+idx]);
mom_def_Cm_poly_paths[9]+=I*ce*cimag(path_storage[9+idx]);
mom_def_Cp_poly_paths[11]+=ce*creal(path_storage[11+idx]);
mom_def_Cm_poly_paths[11]+=I*ce*cimag(path_storage[11+idx]);
mom_def_Cp_poly_paths[13]+=ce*creal(path_storage[13+idx]);
mom_def_Cm_poly_paths[13]+=I*ce*cimag(path_storage[13+idx]);
mom_def_Cp_poly_paths[15]+=ce*creal(path_storage[15+idx]);
mom_def_Cm_poly_paths[15]+=I*ce*cimag(path_storage[15+idx]);
mom_def_Cp_poly_paths[16]+=ce*creal(path_storage[16+idx]);
mom_def_Cm_poly_paths[16]+=I*ce*cimag(path_storage[16+idx]);
mom_def_Cp_poly_paths[17]+=ce*creal(path_storage[17+idx]);
mom_def_Cm_poly_paths[17]+=I*ce*cimag(path_storage[17+idx]);
mom_def_Cp_poly_paths[18]+=ce*creal(path_storage[18+idx]);
mom_def_Cm_poly_paths[18]+=I*ce*cimag(path_storage[18+idx]);
};
for (n_z = 0; n_z < Z; n_z++)
for (n_x = 0; n_x < X; n_x++)
for (n_y = 0; n_y < Y; n_y++)
{
in = ipt(t, n_x, n_y, n_z);
ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
idx = ntors * (n_x + X * (n_y + Y * n_z));
mom_def_Cp_poly_paths[1]+=ce*creal(path_storage[1+idx]);
mom_def_Cm_poly_paths[1]+=I*ce*cimag(path_storage[1+idx]);
mom_def_Cp_poly_paths[4]+=ce*creal(path_storage[4+idx]);
mom_def_Cm_poly_paths[4]+=I*ce*cimag(path_storage[4+idx]);
mom_def_Cp_poly_paths[6]+=ce*creal(path_storage[6+idx]);
mom_def_Cm_poly_paths[6]+=I*ce*cimag(path_storage[6+idx]);
mom_def_Cp_poly_paths[12]+=ce*creal(path_storage[12+idx]);
mom_def_Cm_poly_paths[12]+=I*ce*cimag(path_storage[12+idx]);
mom_def_Cp_poly_paths[14]+=ce*creal(path_storage[14+idx]);
mom_def_Cm_poly_paths[14]+=I*ce*cimag(path_storage[14+idx]);
};
for (n_x = 0; n_x < X; n_x++)
for (n_y = 0; n_y < Y; n_y++)
for (n_z = 0; n_z < Z; n_z++)
{
in = ipt(t, n_x, n_y, n_z);
ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
idx = ntors * (n_x + X * (n_y + Y * n_z));
mom_def_Cp_poly_paths[2]+=ce*creal(path_storage[2+idx]);
mom_def_Cm_poly_paths[2]+=I*ce*cimag(path_storage[2+idx]);
mom_def_Cp_poly_paths[5]+=ce*creal(path_storage[5+idx]);
mom_def_Cm_poly_paths[5]+=I*ce*cimag(path_storage[5+idx]);
mom_def_Cp_poly_paths[7]+=ce*creal(path_storage[7+idx]);
mom_def_Cm_poly_paths[7]+=I*ce*cimag(path_storage[7+idx]);
mom_def_Cp_poly_paths[8]+=ce*creal(path_storage[8+idx]);
mom_def_Cm_poly_paths[8]+=I*ce*cimag(path_storage[8+idx]);
mom_def_Cp_poly_paths[10]+=ce*creal(path_storage[10+idx]);
mom_def_Cm_poly_paths[10]+=I*ce*cimag(path_storage[10+idx]);
};

}
};
void eval_all_torellon_ops(int t, double complex *numerical_tor_out)
{
    static double complex *numerical_op = NULL;
    if (numerical_op == NULL)
    {
        numerical_op = malloc(total_n_tor_op * sizeof(double complex));
    }
    request_space_tors_evaluation();
    eval_time_momentum_torellons(t,0,0,0);
    diPoly_p_0_0_0_Ir_1_C_1_n_1(numerical_op+0);
    eval_time_momentum_torellons(t,1,0,0);
    diPoly_p_1_0_0_Ir_1_C_m1_n_1(numerical_op+1);
    diPoly_p_1_0_0_Ir_3_C_1_n_1(numerical_op+2);
    diPoly_p_1_0_0_Ir_3_C_1_n_2(numerical_op+3);
    eval_time_momentum_torellons(t,1,1,0);
    diPoly_p_1_1_0_Ir_1_C_1_n_1(numerical_op+4);
    for(int i=0;i<total_n_tor_op;i++)
        *(numerical_tor_out+i)+=*(numerical_op+i);
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
    int t1,t2;

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


    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=A1plusOhP Irrep ev=1/1 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 0 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 0; i < 1; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_glue_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_glue_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,0,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 1 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 1; i < 2; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_glue_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_glue_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,0,0) Irrep=E2Dic4 Irrep ev=1/2 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 2 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 2; i < 3; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_glue_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_glue_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,0,0) Irrep=E2Dic4 Irrep ev=2/2 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 3 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 3; i < 4; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_glue_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_glue_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,1,0) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 4 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 4; i < 5; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_glue_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_glue_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }
}
void report_tor_group_setup()
{
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(0,0,0) Irrep=A1plusOhP Charge=+");
lprintf("INIT Measure ML",0," |0=-xyyx|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(1,0,0) Irrep=A1Dic4 Charge=-");
lprintf("INIT Measure ML",0," |1=-yxy|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(1,0,0) Irrep=E2Dic4 Charge=+");
lprintf("INIT Measure ML",0," |2=-zxxz,na|na,3=-yxxy|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(1,1,0) Irrep=A1Dic2 Charge=+");
lprintf("INIT Measure ML",0," |4=-xyyx|");
}
