
#include <stdlib.h>
#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "utils.h"
#include "glueballs.h"
#include <string.h>
#define ntors 57
static hr_complex *tor_path_storage=NULL;
static hr_complex poly0(int in)
{
return polyleg(in,1)->tr;
}

static hr_complex poly1(int in)
{
return polyleg(in,2)->tr;
}

static hr_complex poly2(int in)
{
return polyleg(in,3)->tr;
}

static inline hr_complex diPoly_p_0_0_0_Ir_1_C_1_n_1(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(16.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[0+idx])+(16.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[1+idx])+(16.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[2+idx])+(16.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[0+idx])+(16.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[1+idx])+(16.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[2+idx]);
}

static hr_complex poly29(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger( res, *w2, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly37(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly12(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly13(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly45(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger( res, *w2, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly53(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly16(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly17(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly30(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger( res, *w2, *w1);

site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly38(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly20(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res, res1, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly21(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res, res1, *w1);

site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly46(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger( res, *w2, *w1);

site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly54(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly24(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res, res1, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly25(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res, res1, *w1);

site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly14(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger( res, *w2, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly22(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
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

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly28(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly9(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly41(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger( res, *w2, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly31(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly32(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly49(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly15(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger( res, *w2, *w1);

site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly23(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly35(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res, res1, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly36(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res, res1, *w1);

site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly42(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger( res, *w2, *w1);

site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly39(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res, res1, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly40(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res, res1, *w1);

site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly50(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,2);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly10(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger( res, *w2, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly43(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly44(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly18(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
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

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly26(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger( res, *w2, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly47(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly48(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly33(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly11(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger( res, *w2, *w1);

site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly51(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res, res1, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly52(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res, res1, *w1);

site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly19(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly27(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger( res, *w2, *w1);

site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly55(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res, res1, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly56(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res, res1, *w1);

site = idn_wrk(site, 1);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly34(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

w2 = pu_gauge_wrk(site,3);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
site = idn_wrk(site, 3);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg_dagger(res1, res, *w1);

w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res, res1, *w1);

site = iup_wrk(site, 1);
site = idn_wrk(site, 2);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg_dagger(res1, res, *w1);

wl=polyleg(in,1);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static inline hr_complex diPoly_p_0_0_0_Ir_3_C_1_n_1(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[29+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[37+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[12+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[13+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[45+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[53+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[16+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[17+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[30+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[38+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[20+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[21+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[46+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[54+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[24+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[25+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[14+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[22+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[28+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[9+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[41+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[31+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[32+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[49+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[15+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[23+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[35+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[36+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[42+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[39+idx])+(-2.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[40+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[50+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[10+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[43+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[44+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[18+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[26+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[47+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[48+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[33+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[11+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[51+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[52+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[19+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[27+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[55+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[56+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[34+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[29+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[37+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[12+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[13+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[45+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[53+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[16+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[17+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[30+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[38+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[20+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[21+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[46+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[54+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[24+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[25+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[14+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[22+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[28+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[9+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[41+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[31+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[32+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[49+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[15+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[23+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[35+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[36+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[42+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[39+idx])+(-2.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[40+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[50+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[10+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[43+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[44+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[18+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[26+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[47+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[48+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[33+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[11+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[51+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[52+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[19+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[27+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[55+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[56+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[34+idx]);
}

static inline hr_complex diPoly_p_0_0_0_Ir_3_C_1_n_2(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[12+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[13+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[45+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[53+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[20+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[21+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[46+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[54+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[28+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[9+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[41+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[49+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[35+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[36+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[42+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[50+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[10+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[43+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[44+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[18+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[26+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[47+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[48+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[33+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[11+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[51+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[52+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[19+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[27+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[55+idx])+(1.73205080756887729)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[56+idx])+(-1.73205080756887729)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[34+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[12+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[13+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[45+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[53+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[20+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[21+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[46+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[54+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[28+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[9+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[41+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[49+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[35+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[36+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[42+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[50+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[10+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[43+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[44+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[18+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[26+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[47+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[48+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[33+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[11+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[51+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[52+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[19+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[27+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[55+idx])+(1.73205080756887729)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[56+idx])+(-1.73205080756887729)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[34+idx]);
}

static inline hr_complex diPoly_p_0_0_0_Ir_4_C_1_n_1(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[29+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[37+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[12+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[13+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[45+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[53+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[16+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[17+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[30+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[38+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[20+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[21+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[46+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[54+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[24+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[25+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[14+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[22+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[28+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[9+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[41+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[31+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[32+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[49+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[15+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[23+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[35+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[36+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[42+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[39+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[40+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[50+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[10+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[43+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[44+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[18+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[26+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[47+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[48+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[33+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[11+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[51+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[52+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[19+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[27+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[55+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[56+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[34+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[29+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[37+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[12+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[13+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[45+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[53+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[16+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[17+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[30+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[38+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[20+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[21+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[46+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[54+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[24+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[25+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[14+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[22+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[28+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[9+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[41+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[31+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[32+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[49+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[15+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[23+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[35+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[36+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[42+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[39+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[40+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[50+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[10+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[43+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[44+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[18+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[26+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[47+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[48+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[33+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[11+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[51+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[52+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[19+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[27+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[55+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[56+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[34+idx]);
}

static inline hr_complex diPoly_p_0_0_0_Ir_4_C_1_n_2(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[29+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[37+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[12+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[13+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[45+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[53+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[16+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[17+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[30+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[38+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[20+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[21+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[46+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[54+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[24+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[25+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[14+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[22+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[28+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[9+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[41+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[31+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[32+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[49+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[15+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[23+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[35+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[36+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[42+idx])+(-1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[39+idx])+(1.41421356237309505)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[40+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[50+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[10+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[43+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[44+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[18+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[26+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[47+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[48+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[33+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[11+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[51+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[52+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[19+idx])+(-1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[27+idx])+(-1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[55+idx])+(1.41421356237309505)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[56+idx])+(1.41421356237309505)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[34+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[29+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[37+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[12+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[13+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[45+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[53+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[16+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[17+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[30+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[38+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[20+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[21+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[46+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[54+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[24+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[25+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[14+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[22+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[28+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[9+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[41+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[31+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[32+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[49+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[15+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[23+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[35+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[36+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[42+idx])+(-1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[39+idx])+(1.41421356237309505)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[40+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[50+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[10+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[43+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[44+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[18+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[26+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[47+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[48+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[33+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[11+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[51+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[52+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[19+idx])+(-1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[27+idx])+(-1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[55+idx])+(1.41421356237309505)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[56+idx])+(1.41421356237309505)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[34+idx]);
}

static inline hr_complex diPoly_p_0_0_0_Ir_4_C_1_n_3(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[29+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[37+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[12+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[13+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[45+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[53+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[16+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[17+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[30+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[38+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[20+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[21+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[46+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[54+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[24+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[25+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[14+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[22+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[28+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[9+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[41+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[31+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[32+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[49+idx])+(-1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[15+idx])+(-1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[23+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[35+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[36+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[42+idx])+(1.+I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[39+idx])+(1.-I*1.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[40+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[50+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[10+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[43+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[44+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[18+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[26+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[47+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[48+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[33+idx])+(1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[11+idx])+(-1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[51+idx])+(1.+I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[52+idx])+(1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[19+idx])+(1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[27+idx])+(-1.+I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[55+idx])+(-1.-I*1.)*cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[56+idx])+(-1.-I*1.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[34+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[29+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[37+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[12+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[13+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[45+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[53+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[16+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[17+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[30+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[38+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[20+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[21+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[46+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[54+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[24+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[25+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[14+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[22+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[28+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[9+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[41+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[31+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[32+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[49+idx])+(-1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[15+idx])+(-1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[23+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[35+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[36+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[42+idx])+(1.+I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[39+idx])+(1.-I*1.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[40+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[50+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[10+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[43+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[44+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[18+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[26+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[47+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[48+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[33+idx])+(1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[11+idx])+(-1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[51+idx])+(1.+I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[52+idx])+(1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[19+idx])+(1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[27+idx])+(-1.+I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[55+idx])+(-1.-I*1.)*creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[56+idx])+(-1.-I*1.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[34+idx]);
}

static hr_complex poly3(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
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

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly4(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 2);
w2 = pu_gauge_wrk(site,2);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static inline hr_complex diPoly_p_0_1_0_Ir_1_C_1_n_1(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +(4.)*cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[3+idx])+(4.)*cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[4+idx])+(4.)*creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[3+idx])+(4.)*creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[4+idx]);
}

static hr_complex poly5(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly6(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 1);
w2 = pu_gauge_wrk(site,1);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 3);
w1 = pu_gauge_wrk(site,1);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,3);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly7(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
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

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static hr_complex poly8(int in)
{
suNg *w1, *w2;
suNg res, res1;
int site=in;
hr_complex p;
wilson_lines *wl;

site = idn_wrk(site, 3);
w2 = pu_gauge_wrk(site,3);

_suNg_dagger(res1,*w2);
w2 = &res1;
w1 = pu_gauge_wrk(site,2);
_suNg_times_suNg( res, *w2, *w1);

site = iup_wrk(site, 2);
w1 = pu_gauge_wrk(site,3);
_suNg_times_suNg(res1, res, *w1);

wl=polyleg(in,2);

_suNg_times_suNg(res,res1,wl->p[0]);
_suNg_trace(p,res);
return p;
}

static inline hr_complex diPoly_p_1_1_1_Ir_1_C_1_n_1(int x1, int y1, int z1)
{
    int idx = ntors * (x1 + X * (y1 + Y * z1));
    
return +cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[5+idx])+cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[6+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[3+idx])+cimag(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*cimag(tor_path_storage[4+idx])+cimag(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[7+idx])+cimag(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*cimag(tor_path_storage[8+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[5+idx])+creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[6+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[3+idx])+creal(tor_path_storage[2+ntors * ((x1+X/2)%X + X * ((y1+Y/2)%Y + Y * z1))])*creal(tor_path_storage[4+idx])+creal(tor_path_storage[0+ntors * (x1 + X * ((y1+Y/2)%Y + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[7+idx])+creal(tor_path_storage[1+ntors * ((x1+X/2)%X + X * (y1 + Y *((z1+Z/2)%Z)))])*creal(tor_path_storage[8+idx]);
}

static void eval_time_momentum_torellons(int t, hr_complex * np, hr_complex **pf)
  {
    int nnx, nny, nnz, idx=0, in;
    if(tor_path_storage==NULL)
      {        tor_path_storage = malloc(ntors * X * Y * Z * sizeof(hr_complex));
        for (in = 0; in < ntors * X * Y * Z; in++)
            tor_path_storage[in] = 0.;
    }
for (nny = 0; nny < Y; nny++)
for (nnz = 0; nnz < Z; nnz++)
{
for (nnx = 0; nnx < X; nnx++)
{
in = ipt(t, nnx, nny, nnz);
idx = ntors *(nnx + X * (nny + Y * nnz));
tor_path_storage[0+idx]= poly0(in);
tor_path_storage[3+idx]= poly3(in);
tor_path_storage[7+idx]= poly7(in);
tor_path_storage[9+idx]= poly9(in);
tor_path_storage[26+idx]= poly26(in);
tor_path_storage[27+idx]= poly27(in);
tor_path_storage[28+idx]= poly28(in);
tor_path_storage[33+idx]= poly33(in);
tor_path_storage[34+idx]= poly34(in);
tor_path_storage[35+idx]= poly35(in);
tor_path_storage[36+idx]= poly36(in);
tor_path_storage[41+idx]= poly41(in);
tor_path_storage[42+idx]= poly42(in);
tor_path_storage[43+idx]= poly43(in);
tor_path_storage[44+idx]= poly44(in);
tor_path_storage[49+idx]= poly49(in);
tor_path_storage[50+idx]= poly50(in);
tor_path_storage[51+idx]= poly51(in);
tor_path_storage[52+idx]= poly52(in);
}
pf[0][nny + Y * (nnz + Z * t)] += tor_path_storage[0 + ntors * (X * (nny + Y * nnz))];
}
for (nnz = 0; nnz < Z; nnz++)
for (nnx = 0; nnx < X; nnx++)
{
for (nny = 0; nny < Y; nny++)
{
in = ipt(t, nnx, nny, nnz);
idx = ntors *(nnx + X * (nny + Y * nnz));
tor_path_storage[1+idx]= poly1(in);
tor_path_storage[5+idx]= poly5(in);
tor_path_storage[8+idx]= poly8(in);
tor_path_storage[10+idx]= poly10(in);
tor_path_storage[11+idx]= poly11(in);
tor_path_storage[12+idx]= poly12(in);
tor_path_storage[13+idx]= poly13(in);
tor_path_storage[18+idx]= poly18(in);
tor_path_storage[19+idx]= poly19(in);
tor_path_storage[20+idx]= poly20(in);
tor_path_storage[21+idx]= poly21(in);
tor_path_storage[45+idx]= poly45(in);
tor_path_storage[46+idx]= poly46(in);
tor_path_storage[47+idx]= poly47(in);
tor_path_storage[48+idx]= poly48(in);
tor_path_storage[53+idx]= poly53(in);
tor_path_storage[54+idx]= poly54(in);
tor_path_storage[55+idx]= poly55(in);
tor_path_storage[56+idx]= poly56(in);
}
pf[1][nnx + X * (nnz + Z * t)] += tor_path_storage[1 + ntors *(nnx + X * Y * nnz)];
}
for (nnx = 0; nnx < X; nnx++)
for (nny = 0; nny < Y; nny++)
{
for (nnz = 0; nnz < Z; nnz++)
{
in = ipt(t, nnx, nny, nnz);
idx = ntors * (nnx + X * (nny + Y * nnz));
tor_path_storage[2+idx]= poly2(in);
tor_path_storage[4+idx]= poly4(in);
tor_path_storage[6+idx]= poly6(in);
tor_path_storage[14+idx]= poly14(in);
tor_path_storage[15+idx]= poly15(in);
tor_path_storage[16+idx]= poly16(in);
tor_path_storage[17+idx]= poly17(in);
tor_path_storage[22+idx]= poly22(in);
tor_path_storage[23+idx]= poly23(in);
tor_path_storage[24+idx]= poly24(in);
tor_path_storage[25+idx]= poly25(in);
tor_path_storage[29+idx]= poly29(in);
tor_path_storage[30+idx]= poly30(in);
tor_path_storage[31+idx]= poly31(in);
tor_path_storage[32+idx]= poly32(in);
tor_path_storage[37+idx]= poly37(in);
tor_path_storage[38+idx]= poly38(in);
tor_path_storage[39+idx]= poly39(in);
tor_path_storage[40+idx]= poly40(in);
}
pf[2][nnx + X * (nny  + Y * t)] += tor_path_storage[2 + ntors * (nnx + X * nny )];
}
hr_complex ce = I * 2.0 * 3.141592653589793238462643383279502884197 / GLB_X;for (nnx = 0; nnx < X; nnx++)
for (nny = 0; nny < Y; nny++)
for (nnz = 0; nnz < Z; nnz++)
{
np[0] +=diPoly_p_0_0_0_Ir_1_C_1_n_1(nnx,nny,nnz);
np[1] +=diPoly_p_0_0_0_Ir_3_C_1_n_1(nnx,nny,nnz);
np[2] +=diPoly_p_0_0_0_Ir_3_C_1_n_2(nnx,nny,nnz);
np[3] +=diPoly_p_0_0_0_Ir_4_C_1_n_1(nnx,nny,nnz);
np[4] +=diPoly_p_0_0_0_Ir_4_C_1_n_2(nnx,nny,nnz);
np[5] +=diPoly_p_0_0_0_Ir_4_C_1_n_3(nnx,nny,nnz);
np[6] +=cexp(-ce * (double)(nny))*diPoly_p_0_1_0_Ir_1_C_1_n_1(nnx,nny,nnz);
np[7] +=cexp(-ce * (double)(nnx + nny + nnz))*diPoly_p_1_1_1_Ir_1_C_1_n_1(nnx,nny,nnz);
}
}
void eval_all_torellon_ops(int t, hr_complex *numerical_tor_out, hr_complex ** polyf)
{
    static hr_complex *numerical_op = NULL;
    if (numerical_op == NULL)
    {
        numerical_op = malloc(total_n_tor_op * sizeof(hr_complex));
    } 
   for (int i = 0; i < total_n_tor_op; i++)
       numerical_op[i] = 0;
eval_time_momentum_torellons(t,numerical_op,polyf);
    for(int i=0;i<total_n_tor_op;i++)
       numerical_tor_out[i]+= numerical_op[i];
}

void collect_1pt_torellon_functions(cor_list *lcor, hr_complex *tor_storage, hr_complex ** polyf)
{
    int n1, n2, n3, i;
    static hr_complex *tor1_bf;
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

        tor1_bf = malloc(sizeof(hr_complex) * total_n_tor_op * n_total_active_slices);

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

    static hr_complex *tor2;
    static hr_complex *tor1;
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
                    memcpy(tor1, tor_storage + t1 * total_n_tor_op, sizeof(hr_complex) * total_n_tor_op);
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
                    memcpy(tor2, tor_storage + t2 * total_n_tor_op, sizeof(hr_complex) * total_n_tor_op);
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
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=EplusOhP Irrep ev=1/2 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 1 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 1; i < 2; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=EplusOhP Irrep ev=2/2 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 2 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 2; i < 3; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=T1plusOhP Irrep ev=1/3 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 3 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 3; i < 4; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=T1plusOhP Irrep ev=2/3 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 4 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 4; i < 5; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,0,0) Irrep=T1plusOhP Irrep ev=3/3 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 5 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 5; i < 6; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(0,1,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 6 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 6; i < 7; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1ptTor function P=(1,1,1) Irrep=A1Dic3 Irrep ev=1/1 Charge=+ nop=%d\n",1 );
    lprintf("Measure ML", 0, "Tor id= 7 \n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
                for (i = 7; i < 8; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(tor1_bf[i + total_n_tor_op  * listactive[n1]]),
                            cimag(tor1_bf[i + total_n_tor_op  * listactive[n1]]));
            lprintf("Measure ML", 0, "\n");
        }

    if (polyf == NULL)
          return;

    hr_complex *lpoly = NULL;
    hr_complex *gpoly = NULL;
    hr_complex *pcor = NULL;
    if (lpoly == NULL)
    {
        lpoly = malloc(T * sizeof(hr_complex));
#ifdef WITH_MPI
        gpoly = malloc(GLB_T * sizeof(hr_complex));
#else
        gpoly = lpoly;
#endif
        pcor = malloc(GLB_T * sizeof(hr_complex));
    }
    for (n1 = 0; n1 < GLB_T; n1++)
        pcor[n1] = 0.;

    for (n1 = 0; n1 < Y; n1++)
        for (n2 = 0; n2 < Z; n2++)
        {
            for (n3 = 0; n3 < T; n3++)
                lpoly[n3] = polyf[0][(n1 + Y * (n2 + Z * n3))];
#ifdef WITH_MPI
            MPI_Gather((double *)lpoly, 2 * T, MPI_DOUBLE, (double *)gpoly, 2 * T, MPI_DOUBLE, 0, GLB_COMM);
#endif
            for (i = 0; i < lcor->n_entries; i++)
                pcor[abs(lcor->list[i].t2 - lcor->list[i].t1)] += conj(gpoly[lcor->list[i].t1]) * gpoly[lcor->list[i].t2] / (3.0 * lcor->list[i].n_pairs * Y * Z);
        }
    for (n1 = 0; n1 < X; n1++)
        for (n2 = 0; n2 < Z; n2++)
        {
            for (n3 = 0; n3 < T; n3++)
                lpoly[n3] = polyf[1][(n1 + X * (n2 + Z * n3))];
#ifdef WITH_MPI
            MPI_Gather((double *)lpoly, 2 * T, MPI_DOUBLE, (double *)gpoly, 2 * T, MPI_DOUBLE, 0, GLB_COMM);
#endif

            for (i = 0; i < lcor->n_entries; i++)
                pcor[abs(lcor->list[i].t2 - lcor->list[i].t1)] += conj(gpoly[lcor->list[i].t1]) * gpoly[lcor->list[i].t2] / (3.0 * lcor->list[i].n_pairs * X * Z);
        }
    for (n1 = 0; n1 < X; n1++)
        for (n2 = 0; n2 < Y; n2++)
        {
            for (n3 = 0; n3 < T; n3++)
                lpoly[n3] = polyf[2][(n1 + X * (n2 + Y * n3))];
#ifdef WITH_MPI
            MPI_Gather((double *)lpoly, 2 * T, MPI_DOUBLE, (double *)gpoly, 2 * T, MPI_DOUBLE, 0, GLB_COMM);
#endif
            for (i = 0; i < lcor->n_entries; i++)
                pcor[abs(lcor->list[i].t2 - lcor->list[i].t1)] += conj(gpoly[lcor->list[i].t1]) * gpoly[lcor->list[i].t2] / (3.0 * lcor->list[i].n_pairs * X * Y);
        }
   for (n1 = 0; n1 < GLB_T; n1++)
        lprintf("Measure ML", 0, " Polyakov Cor dt=%d ( %.10e %.10e )\n", n1, creal(pcor[n1]), cimag(pcor[n1]));
}
void report_tor_group_setup()
{
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(0,0,0) Irrep=A1plusOhP Charge=+");
lprintf("INIT Measure ML",0," |0=|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(0,0,0) Irrep=EplusOhP Charge=+");
lprintf("INIT Measure ML",0," |1=-x-yxzy,2=-xy-zxz|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(0,0,0) Irrep=T1plusOhP Charge=+");
lprintf("INIT Measure ML",0," |3=-x-yxzy,4=-x-yxzy,5=-x-yxzy|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(0,1,0) Irrep=A1Dic4 Charge=+");
lprintf("INIT Measure ML",0," |6=-yxy|");
lprintf("INIT Measure ML",0,"\n1pt_tor Irrep multiplets Total P=(1,1,1) Irrep=A1Dic3 Charge=+");
lprintf("INIT Measure ML",0," |7=-xyx|");
}
