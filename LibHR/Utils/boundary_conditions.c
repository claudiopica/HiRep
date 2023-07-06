/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *
* All rights reserved.                                                      *
\***************************************************************************/

#include "utils.h"
#include "libhr_core.h"
#include "memory.h"
#include "io.h"
#include <math.h>

static int init = 0;
static BCs_pars_t BCs_pars;

#ifdef PLAQ_WEIGHTS
#ifdef GAUGE_SPATIAL_TWIST
static void init_plaq_twisted_BCs();
#endif

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
static void init_plaq_SF_BCs(double ct);
static void init_gf_SF_BCs(suNg *dn, suNg *up);
#endif
#endif

void init_BCs(BCs_pars_t *pars) {
    error(init == 1, 1, "init_BCs [boundary_conditions.c]", "BCs already initialized");
    init = 1;

#ifdef PLAQ_WEIGHTS
    plaq_weight = malloc(sizeof(double) * glattice.gsize_gauge * 16);
    rect_weight = malloc(sizeof(double) * glattice.gsize_gauge * 16);
    for (int i = 0; i < 16 * glattice.gsize_gauge; i++) {
        rect_weight[i] = 1.0;
        plaq_weight[i] = 1.0;
    }
#endif

    BCs_pars.fermion_twisting_theta[0] = 0.;
    BCs_pars.fermion_twisting_theta[1] = 0.;
    BCs_pars.fermion_twisting_theta[2] = 0.;
    BCs_pars.fermion_twisting_theta[3] = 0.;
    BCs_pars.gauge_boundary_improvement_cs = 1.;
    BCs_pars.gauge_boundary_improvement_ct = 1.;
    BCs_pars.chiSF_boundary_improvement_ds = 1.;
    BCs_pars.SF_BCs = 0;

    if (pars != NULL) { BCs_pars = *pars; }

    lprintf("BCS", 0,
            "Gauge field: "
#if defined(BC_T_OPEN)
            "OPEN"
#elif defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
            "DIRICHLET"
#else
            "PERIODIC"
#endif
            " x "
#if defined(GAUGE_SPATIAL_TWIST)
            "TWISTED TWISTED TWISTED"
#else
            "PERIODIC PERIODIC PERIODIC"
#endif
            "\n");

    lprintf("BCS", 0,
            "Fermion fields: "
#if defined(BC_T_SF_ROTATED)
            "OPEN"
#elif defined(BC_T_OPEN) || defined(BC_T_SF)
            "DIRICHLET"
#elif defined(BC_T_ANTIPERIODIC)
            "ANTIPERIODIC"
#elif defined(BC_T_THETA)
            "THETA"
#else
            "PERIODIC"
#endif
            " x "
#if defined(BC_X_ANTIPERIODIC)
            "ANTIPERIODIC "
#elif defined(BC_X_THETA)
            "THETA "
#elif defined(GAUGE_SPATIAL_TWIST)
            "TWISTED "
#else
            "PERIODIC "
#endif
#if defined(BC_Y_ANTIPERIODIC)
            "ANTIPERIODIC "
#elif defined(BC_Y_THETA)
            "THETA "
#elif defined(GAUGE_SPATIAL_TWIST)
            "TWISTED "
#else
            "PERIODIC "
#endif
#if defined(BC_Z_ANTIPERIODIC)
            "ANTIPERIODIC"
#elif defined(BC_Z_THETA)
            "THETA"
#elif defined(GAUGE_SPATIAL_TWIST)
            "TWISTED"
#else
            "PERIODIC"
#endif
            "\n");

#if defined(BC_T_SF_ROTATED)
    lprintf("BCS", 0, "Chirally rotated Schroedinger Functional ds=%e BCs=%d(1=Background, 0=no Background)\n",
            BCs_pars.chiSF_boundary_improvement_ds, BCs_pars.SF_BCs);
#elif defined(BC_T_SF)
    lprintf("BCS", 0, "Basic Schroedinger Functional BCs=%d (1=Background, 0=no Background)\n", BCs_pars.SF_BCs);
#endif

#ifdef FERMION_THETA

#ifndef BC_T_THETA
    BCs_pars.fermion_twisting_theta[0] = 0.;
#endif
#ifndef BC_X_THETA
    BCs_pars.fermion_twisting_theta[1] = 0.;
#endif
#ifndef BC_Y_THETA
    BCs_pars.fermion_twisting_theta[2] = 0.;
#endif
#ifndef BC_Z_THETA
    BCs_pars.fermion_twisting_theta[3] = 0.;
#endif

    lprintf("BCS", 0, "Fermion twisting theta angles = %e %e %e %e\n", BCs_pars.fermion_twisting_theta[0],
            BCs_pars.fermion_twisting_theta[1], BCs_pars.fermion_twisting_theta[2], BCs_pars.fermion_twisting_theta[3]);

    eitheta[0] =
        cos(BCs_pars.fermion_twisting_theta[0] / (double)GLB_T) + I * sin(BCs_pars.fermion_twisting_theta[0] / (double)GLB_T);
    eitheta[1] =
        cos(BCs_pars.fermion_twisting_theta[1] / (double)GLB_X) + I * sin(BCs_pars.fermion_twisting_theta[1] / (double)GLB_X);
    eitheta[2] =
        cos(BCs_pars.fermion_twisting_theta[2] / (double)GLB_Y) + I * sin(BCs_pars.fermion_twisting_theta[2] / (double)GLB_Y);
    eitheta[3] =
        cos(BCs_pars.fermion_twisting_theta[3] / (double)GLB_Z) + I * sin(BCs_pars.fermion_twisting_theta[3] / (double)GLB_Z);

#endif /* FERMION_THETA */

#ifdef BC_T_OPEN
    lprintf("BCS", 0, "Open BC gauge boundary term ct=%e cs=%e\n", BCs_pars.gauge_boundary_improvement_ct,
            BCs_pars.gauge_boundary_improvement_cs);
    init_plaq_open_BCs(plaq_weight, rect_weight, BCs_pars.gauge_boundary_improvement_ct,
                       BCs_pars.gauge_boundary_improvement_cs);
#endif

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
    lprintf("BCS", 0, "BC gauge boundary term ct=%e\n", BCs_pars.gauge_boundary_improvement_ct);
    if (BCs_pars.SF_BCs == 0) {
        _suNg_unit(BCs_pars.gauge_boundary_dn);
        _suNg_unit(BCs_pars.gauge_boundary_up);
    } else {
        init_gf_SF_BCs(&(BCs_pars.gauge_boundary_dn), &(BCs_pars.gauge_boundary_up));
    }
    init_plaq_SF_BCs(BCs_pars.gauge_boundary_improvement_ct);
#endif

#ifdef GAUGE_SPATIAL_TWIST
    init_plaq_twisted_BCs();
#endif
}

void free_BCs(void) {
    if (init == 0) { return; }

    init = 0;

#ifdef PLAQ_WEIGHTS
    if (plaq_weight != NULL) { free(plaq_weight); }
    if (rect_weight != NULL) { free(rect_weight); }
#endif
}

#ifdef BC_T_ANTIPERIODIC
static void sp_T_antiperiodic_BCs();
#endif
#ifdef BC_X_ANTIPERIODIC
static void sp_X_antiperiodic_BCs();
#endif
#ifdef BC_Y_ANTIPERIODIC
static void sp_Y_antiperiodic_BCs();
#endif
#ifdef BC_Z_ANTIPERIODIC
static void sp_Z_antiperiodic_BCs();
#endif
/*static void sp_spatial_theta_BCs(double theta);*/
#ifdef BC_T_SF_ROTATED
static void chiSF_ds_BT(double ds);
#endif

void apply_BCs_on_represented_gauge_field() {
#ifdef BC_T_ANTIPERIODIC
    sp_T_antiperiodic_BCs();
#endif
#ifdef BC_X_ANTIPERIODIC
    sp_X_antiperiodic_BCs();
#endif
#ifdef BC_Y_ANTIPERIODIC
    sp_Y_antiperiodic_BCs();
#endif
#ifdef BC_Z_ANTIPERIODIC
    sp_Z_antiperiodic_BCs();
#endif
#ifdef BC_T_SF_ROTATED
#ifndef ALLOCATE_REPR_GAUGE_FIELD
#error The represented gauge field must be allocated!!!
#endif
    chiSF_ds_BT(BCs_pars.chiSF_boundary_improvement_ds);
#endif
}

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
static void gf_SF_BCs(suNg *dn, suNg *up);
#endif
#ifdef BC_T_OPEN
static void gf_open_BCs();
#endif

void apply_BCs_on_fundamental_gauge_field(void) {
    complete_sendrecv_gfield(u_gauge);
#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
    gf_SF_BCs(&BCs_pars.gauge_boundary_dn, &BCs_pars.gauge_boundary_up);
#endif
#ifdef BC_T_OPEN
    gf_open_BCs();
#endif
}

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
static void mf_Dirichlet_BCs(suNg_av_field *force);
#endif
#ifdef BC_T_OPEN
static void mf_open_BCs(suNg_av_field *force);
#endif

void apply_BCs_on_momentum_field(suNg_av_field *force) {
#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
    mf_Dirichlet_BCs(force);
#endif
#ifdef BC_T_OPEN
    mf_open_BCs(force);
#endif
}

#if defined(BC_T_SF)
static void sf_Dirichlet_BCs(spinor_field *sp);
static void sf_Dirichlet_BCs_flt(spinor_field_flt *sp);
#endif
#if defined(BC_T_SF_ROTATED)
static void sf_open_BCs(spinor_field *sp);
static void sf_open_BCs_flt(spinor_field_flt *sp);
#endif
#ifdef BC_T_OPEN
static void sf_open_v2_BCs(spinor_field *sf);
static void sf_open_v2_BCs_flt(spinor_field_flt *sf);
#endif

void apply_BCs_on_spinor_field(spinor_field *sp) {
#ifdef BC_T_SF
    sf_Dirichlet_BCs(sp);
#endif
#ifdef BC_T_SF_ROTATED
    sf_open_BCs(sp);
#endif
#ifdef BC_T_OPEN
    sf_open_v2_BCs(sp);
#endif
}

void apply_BCs_on_spinor_field_flt(spinor_field_flt *sp) {
#if defined(BC_T_SF)
    sf_Dirichlet_BCs_flt(sp);
#endif
#if defined(BC_T_SF_ROTATED)
    sf_open_BCs_flt(sp);
#endif
#ifdef BC_T_OPEN
    sf_open_v2_BCs_flt(sp);
#endif
}

#if (defined(WITH_EXPCLOVER) || defined(WITH_CLOVER)) && defined(BC_T_OPEN)
static void cl_open_BCs(suNfc_field *);
#endif
#if (defined(WITH_EXPCLOVER) || defined(WITH_CLOVER)) && (defined(BC_T_SF) || defined(BC_T_SF_ROTATED))
static void cl_SF_BCs(suNfc_field *);
#endif

void apply_BCs_on_clover_term(suNfc_field *cl) {
#if (defined(WITH_EXPCLOVER) || defined(WITH_CLOVER)) && defined(BC_T_OPEN)
    cl_open_BCs(cl);
#endif

#if (defined(WITH_EXPCLOVER) || defined(WITH_CLOVER)) && (defined(BC_T_SF) || defined(BC_T_SF_ROTATED))
    cl_SF_BCs(cl);
#endif
}

/***************************************************************************/
/* BOUNDARY CONDITIONS TO BE APPLIED ON THE REPRESENTED GAUGE FIELD        */
/***************************************************************************/

#ifdef BC_T_ANTIPERIODIC
static void sp_T_antiperiodic_BCs() {
    if (COORD[0] == 0) {
        int index;
        int ix, iy, iz;
        suNf *u;
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(2 * T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge_f(index, 0);
                        _suNf_minus(*u, *u);
                    }
                }
            }
        }
    }
}
#endif
#ifdef BC_X_ANTIPERIODIC
static void sp_X_antiperiodic_BCs() {
    if (COORD[1] == 0) {
        int index;
        int it, iy, iz;
        suNf *u;
        for (it = 0; it < T_EXT; ++it) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(it, 2 * X_BORDER, iy, iz);
                    if (index != -1) {
                        u = pu_gauge_f(index, 1);
                        _suNf_minus(*u, *u);
                    }
                }
            }
        }
    }
}
#endif
#ifdef BC_Y_ANTIPERIODIC
static void sp_Y_antiperiodic_BCs() {
    if (COORD[2] == 0) {
        int index;
        int ix, it, iz;
        suNf *u;
        for (it = 0; it < T_EXT; ++it) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(it, ix, 2 * Y_BORDER, iz);
                    if (index != -1) {
                        u = pu_gauge_f(index, 2);
                        _suNf_minus(*u, *u);
                    }
                }
            }
        }
    }
}
#endif
#ifdef BC_Z_ANTIPERIODIC
static void sp_Z_antiperiodic_BCs() {
    if (COORD[3] == 0) {
        int index;
        int ix, iy, it;
        suNf *u;
        for (it = 0; it < T_EXT; ++it) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    index = ipt_ext(it, ix, iy, 2 * Z_BORDER);
                    if (index != -1) {
                        u = pu_gauge_f(index, 3);
                        _suNf_minus(*u, *u);
                    }
                }
            }
        }
    }
}
#endif

/*
#ifndef REPR_ADJOINT
static void sp_spatial_theta_BCs(double theta) {
  complex phase[3];
  phase[0].re=cos(theta/GLB_X);
  phase[0].im=sin(theta/GLB_X);
  phase[1].re=cos(theta/GLB_Y);
  phase[1].im=sin(theta/GLB_Y);
  phase[2].re=cos(theta/GLB_Z);
  phase[2].im=sin(theta/GLB_Z);

  int index,it,ix,iy,iz,mu;
  suNf Rtmp;
  suNf *Ru;
  for (it=0;it<T_EXT;++it) for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
    index=ipt_ext(it,ix,iy,iz);
    if(index!=-1) {
      for (mu=1;mu<4;mu++) {
        Ru = pu_gauge_f(index,mu);
        Rtmp = *Ru;
        _suNf_mulc(*Ru,phase[mu-1],Rtmp);
      }
    }
  }
}
#endif
*/

#ifdef BC_T_SF_ROTATED
static void chiSF_ds_BT(double ds) {
    if (COORD[0] == 0) {
        int index;
        int ix, iy, iz;
        suNf *u;
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge_f(index, 1);
                        _suNf_mul(*u, ds, *u);
                        u = pu_gauge_f(index, 2);
                        _suNf_mul(*u, ds, *u);
                        u = pu_gauge_f(index, 3);
                        _suNf_mul(*u, ds, *u);
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        int index;
        int ix, iy, iz;
        suNf *u;
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge_f(index, 1);
                        _suNf_mul(*u, ds, *u);
                        u = pu_gauge_f(index, 2);
                        _suNf_mul(*u, ds, *u);
                        u = pu_gauge_f(index, 3);
                        _suNf_mul(*u, ds, *u);
                    }
                }
            }
        }
    }
}
#endif

/***************************************************************************/
/* BOUNDARY CONDITIONS TO BE APPLIED ON THE FUNDAMENTAL GAUGE FIELD        */
/***************************************************************************/

#define PI 3.141592653589793238462643383279502884197
#define ST 1.414213562373095048801688724209698078570

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
#ifndef GAUGE_SON

#if NG == 2

#ifndef HALFBG_SF
static double SF_eta = PI / 4.0;
static double SF_phi0_up[NG] = { -PI, PI };
#else
static double SF_eta = PI / 8.0;
static double SF_phi0_up[NG] = { -PI / 2, PI / 2 };
#endif

static double SF_phi0_dn[NG] = { 0., 0. };
static double SF_phi1_dn[NG] = { -1., 1. };
static double SF_phi1_up[NG] = { 1., -1. };

#elif NG == 3

static double SF_eta = 0.;
static double SF_phi0_dn[NG] = { -PI / 3., 0., PI / 3. };
static double SF_phi1_dn[NG] = { 1., -.5, -.5 };
//  static double SF_phi0_up[NG] = {-PI, PI / 3., 2. * PI / 3.};
static double SF_phi0_up[NG] = { -5. * PI / 3., 2. * PI / 3., PI };
static double SF_phi1_up[NG] = { -1., .5, .5 };

#elif NG == 4

static double SF_eta = 0.;
static double SF_phi0_dn[NG] = { -ST * PI / 4., ST *PI / 4. - PI / 2., PI / 2. - ST *PI / 4., ST *PI / 4. };
static double SF_phi1_dn[NG] = { -.5, -.5, .5, .5 };
static double SF_phi0_up[NG] = { -ST * PI / 4. - PI / 2., -PI + ST *PI / 4., PI - ST *PI / 4., PI / 2. + ST *PI / 4. };
static double SF_phi1_up[NG] = { .5, .5, -.5, -.5 };

#endif

static void init_gf_SF_BCs(suNg *dn, suNg *up) {
#if (NG > 4)

    error(0, 1, "init_gf_SF_BCs " __FILE__, "SF boundary conditions not defined at this NG");

#else
    int k;

    _suNg_zero(*dn);
    for (k = 0; k < NG; k++) {
        dn->c[(1 + NG) * k] = cos((SF_phi0_dn[k] + SF_phi1_dn[k] * SF_eta) / (GLB_T - 2)) +
                              I * sin((SF_phi0_dn[k] + SF_phi1_dn[k] * SF_eta) / (GLB_T - 2));
    }

    _suNg_zero(*up);
    for (k = 0; k < NG; k++) {
        up->c[(1 + NG) * k] = cos((SF_phi0_up[k] + SF_phi1_up[k] * SF_eta) / (GLB_T - 2)) +
                              I * sin((SF_phi0_up[k] + SF_phi1_up[k] * SF_eta) / (GLB_T - 2));
    }

#if defined(BC_T_SF)
    lprintf("BCS", 0, "SF boundary phases phi0  ( ");
    for (k = 0; k < NG; k++) {
        lprintf("BCS", 0, "%lf ", SF_phi0_dn[k] / 2.);
    }
    lprintf("BCS", 0, ")\n");
    lprintf("BCS", 0, "SF boundary phases phi0'  ( ");
    for (k = 0; k < NG; k++) {
        lprintf("BCS", 0, "%lf ", SF_phi0_up[k] / 2.);
    }
    lprintf("BCS", 0, ")\n");
#else
    lprintf("BCS", 0, "SF boundary phases phi0'  ( ");
    for (k = 0; k < NG; k++) {
        lprintf("BCS", 0, "%lf ", SF_phi0_up[k] / 2.);
    }
    lprintf("BCS", 0, ")\n");
    lprintf("BCS", 0, "SF boundary phases phi1'  ( ");
    for (k = 0; k < NG; k++) {
        lprintf("BCS", 0, "%lf ", SF_phi1_up[k]);
    }
    lprintf("BCS", 0, ")\n");

    lprintf("BCS", 0, "SF boundary phases phi0  ( ");
    for (k = 0; k < NG; k++) {
        lprintf("BCS", 0, "%lf ", SF_phi0_dn[k] / 2.);
    }
    lprintf("BCS", 0, ")\n");
    lprintf("BCS", 0, "SF boundary phases phi1  ( ");
    for (k = 0; k < NG; k++) {
        lprintf("BCS", 0, "%lf ", SF_phi1_dn[k]);
    }
    lprintf("BCS", 0, ")\n");
#endif
#endif
}
#else
static void init_gf_SF_BCs(suNg *dn, suNg *up) {
    error(1, 1, "init_gf_SF_BCs", "SF not implemented for real gauge group");
}
#endif
#endif

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
static void gf_SF_BCs(suNg *dn, suNg *up) {
    int index;
    int ix, iy, iz;
    suNg *u;

    if (COORD[0] == 0) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            u = pu_gauge(index, 0);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 1);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 2);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 3);
                            _suNg_unit(*u);
                        }
                    }
                }
            }
        }
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge(index, 0);
                        _suNg_unit(*u);
                        u = pu_gauge(index, 1);
                        _suNg_unit(*u);
                        u = pu_gauge(index, 2);
                        _suNg_unit(*u);
                        u = pu_gauge(index, 3);
                        _suNg_unit(*u);
                    }
                }
            }
        }
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge(index, 1);
                        *u = *dn;
                        u = pu_gauge(index, 2);
                        *u = *dn;
                        u = pu_gauge(index, 3);
                        *u = *dn;
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge(index, 0);
                        _suNg_unit(*u);
                        u = pu_gauge(index, 1);
                        *u = *up;
                        u = pu_gauge(index, 2);
                        *u = *up;
                        u = pu_gauge(index, 3);
                        *u = *up;
                    }
                }
            }
        }
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            u = pu_gauge(index, 0);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 1);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 2);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 3);
                            _suNg_unit(*u);
                        }
                    }
                }
            }
        }
    }
}

void SF_classical_solution(void) {
    double x0;
    int ix, iy, iz, it, k, index;
    suNg U, *u;

    for (it = 0; it < T_EXT; ++it) {
        x0 = (double)(zerocoord[0] + it - T_BORDER);

        if (x0 >= 1 && x0 < GLB_T - 1) {
            _suNg_zero(U);

            for (k = 0; k < NG; k++) {
                U.c[(1 + NG) * k] = cos(((SF_phi0_dn[k] + SF_phi1_dn[k] * SF_eta) * (GLB_T - 1.0 - x0) +
                                         (SF_phi0_up[k] + SF_phi1_up[k] * SF_eta) * (x0 - 1.0)) /
                                        (GLB_T - 2.0) / (GLB_T - 2.0)) +
                                    I * sin(((SF_phi0_dn[k] + SF_phi1_dn[k] * SF_eta) * (GLB_T - 1.0 - x0) +
                                             (SF_phi0_up[k] + SF_phi1_up[k] * SF_eta) * (x0 - 1.0)) /
                                            (GLB_T - 2.0) / (GLB_T - 2));
            }

            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(it, ix, iy, iz);
                        if (index != -1) {
                            u = pu_gauge(index, 0);
                            _suNg_unit(*u);
                            u = pu_gauge(index, 1);
                            *u = U;
                            u = pu_gauge(index, 2);
                            *u = U;
                            u = pu_gauge(index, 3);
                            *u = U;
                        }
                    }
                }
            }
        }
    }
    apply_BCs_on_fundamental_gauge_field();
}
#endif

#ifdef BC_T_OPEN
static void gf_open_BCs() {
    int index;
    int ix, iy, iz;
    suNg *u;

    if (COORD[0] == 0) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            u = pu_gauge(index, 0);
                            _suNg_zero(*u);
                            u = pu_gauge(index, 1);
                            _suNg_zero(*u);
                            u = pu_gauge(index, 2);
                            _suNg_zero(*u);
                            u = pu_gauge(index, 3);
                            _suNg_zero(*u);
                        }
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        u = pu_gauge(index, 0);
                        _suNg_zero(*u);
                    }
                }
            }
        }
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            u = pu_gauge(index, 0);
                            _suNg_zero(*u);
                            u = pu_gauge(index, 1);
                            _suNg_zero(*u);
                            u = pu_gauge(index, 2);
                            _suNg_zero(*u);
                            u = pu_gauge(index, 3);
                            _suNg_zero(*u);
                        }
                    }
                }
            }
        }
    }
}
#endif

/***************************************************************************/
/* BOUNDARY CONDITIONS TO BE APPLIED ON THE MOMENTUM FIELDS                */
/***************************************************************************/

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
static void mf_Dirichlet_BCs(suNg_av_field *force) {
    int ix, iy, iz, index;

    if (COORD[0] == 0) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                        }
                    }
                }
            }
        }
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                    }

                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1) {
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                    }
                }
            }
        }
    }

    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                        _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                    }
                }
            }
        }
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                        }
                    }
                }
            }
        }
    }
}
#endif
#ifdef BC_T_OPEN
static void mf_open_BCs(suNg_av_field *force) {
    int ix, iy, iz, index;

    if (COORD[0] == 0) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                        }
                    }
                }
            }
        }
    }

    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) { _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0)); }
                }
            }
        }
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 0));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 1));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 2));
                            _algebra_vector_zero_g(*_4FIELD_AT(force, index, 3));
                        }
                    }
                }
            }
        }
    }
}
#endif

/***************************************************************************/
/* BOUNDARY CONDITIONS TO BE APPLIED ON THE CLOVER TERM                    */
/***************************************************************************/
#if (defined(WITH_EXPCLOVER) || defined(WITH_CLOVER)) && (defined(BC_T_SF) || defined(BC_T_SF_ROTATED))
static void cl_SF_BCs(suNfc_field *cl) {
    int index;
    suNfc u;
    _suNfc_zero(u);

    // These should reflect the boundary conditions imposed on the spinor fields
    if (COORD[0] == 0) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    if (T_BORDER > 0) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            for (int mu = 0; mu < 4; mu++) {
                                *_4FIELD_AT(cl, index, mu) = u;
                            }
                        }
                    }
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        for (int mu = 0; mu < 4; mu++) {
                            *_4FIELD_AT(cl, index, mu) = u;
                        }
                    }
                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1) {
                        for (int mu = 0; mu < 4; mu++) {
                            *_4FIELD_AT(cl, index, mu) = u;
                        }
                    }
                }
            }
        }
    }
}
#endif
#if (defined(WITH_EXPCLOVER) || defined(WITH_CLOVER)) && defined(BC_T_OPEN)
static void cl_open_BCs(suNfc_field *cl) {
    int index;
    suNfc u;
    _suNfc_zero(u);

    // These should reflect the boundary conditions imposed on the spinor fields
    if (COORD[0] == 0) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    if (T_BORDER > 0) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            for (int mu = 0; mu < 4; mu++) {
                                *_4FIELD_AT(cl, index, mu) = u;
                            }
                        }
                    }
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        for (int mu = 0; mu < 4; mu++) {
                            *_4FIELD_AT(cl, index, mu) = u;
                        }
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        for (int mu = 0; mu < 4; mu++) {
                            *_4FIELD_AT(cl, index, mu) = u;
                        }
                    }
                    if (T_BORDER > 0) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            for (int mu = 0; mu < 4; mu++) {
                                *_4FIELD_AT(cl, index, mu) = u;
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif
/***************************************************************************/
/* BOUNDARY CONDITIONS TO BE APPLIED ON THE SPINOR FIELDS                  */
/***************************************************************************/
#if defined(BC_T_SF)
static void sf_Dirichlet_BCs(spinor_field *sp) {
    int ix, iy, iz, index;
    if (COORD[0] == 0) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                }
            }
        }
    }
}

static void sf_Dirichlet_BCs_flt(spinor_field_flt *sp) {
    int ix, iy, iz, index;
    if (COORD[0] == 0) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                }
            }
        }
    }
}
#endif
#if defined(BC_T_SF_ROTATED)
static void sf_open_BCs(spinor_field *sp) {
    int ix, iy, iz, index;
    if (COORD[0] == 0) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                }
            }
        }
    }
}

static void sf_open_BCs_flt(spinor_field_flt *sp) {
    int ix, iy, iz, index;
    if (COORD[0] == 0) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1 && sp->type->master_shift <= index &&
                        sp->type->master_shift + sp->type->gsize_spinor > index) {
                        _spinor_zero_f(*_FIELD_AT(sp, index));
                    }
                }
            }
        }
    }
}
#endif

#ifdef BC_T_OPEN
static void sf_open_v2_BCs(spinor_field *sf) {
    int index;
    suNf_spinor u;
    _spinor_zero_f(u);

    // These should reflect the boundary conditions imposed on the clover field
    if (COORD[0] == 0) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    if (T_BORDER > 0) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1 && sf->type->master_shift <= index &&
                            sf->type->master_shift + sf->type->gsize_spinor > index) {
                            *_FIELD_AT(sf, index) = u;
                        }
                    }
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1 && sf->type->master_shift <= index &&
                        sf->type->master_shift + sf->type->gsize_spinor > index) {
                        *_FIELD_AT(sf, index) = u;
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1 && sf->type->master_shift <= index &&
                        sf->type->master_shift + sf->type->gsize_spinor > index) {
                        *_FIELD_AT(sf, index) = u;
                    }
                    if (T_BORDER > 0) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1 && sf->type->master_shift <= index &&
                            sf->type->master_shift + sf->type->gsize_spinor > index) {
                            *_FIELD_AT(sf, index) = u;
                        }
                    }
                }
            }
        }
    }
}

static void sf_open_v2_BCs_flt(spinor_field_flt *sf) {
    int index;
    suNf_spinor_flt u;
    _spinor_zero_f(u);

    // These should reflect the boundary conditions imposed on the clover field
    if (COORD[0] == 0) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    if (T_BORDER > 0) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1 && sf->type->master_shift <= index &&
                            sf->type->master_shift + sf->type->gsize_spinor > index) {
                            *_FIELD_AT(sf, index) = u;
                        }
                    }
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1 && sf->type->master_shift <= index &&
                        sf->type->master_shift + sf->type->gsize_spinor > index) {
                        *_FIELD_AT(sf, index) = u;
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (int ix = 0; ix < X_EXT; ix++) {
            for (int iy = 0; iy < Y_EXT; iy++) {
                for (int iz = 0; iz < Z_EXT; iz++) {
                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1 && sf->type->master_shift <= index &&
                        sf->type->master_shift + sf->type->gsize_spinor > index) {
                        *_FIELD_AT(sf, index) = u;
                    }
                    if (T_BORDER > 0) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1 && sf->type->master_shift <= index &&
                            sf->type->master_shift + sf->type->gsize_spinor > index) {
                            *_FIELD_AT(sf, index) = u;
                        }
                    }
                }
            }
        }
    }
}
#endif

/***************************************************************************/
/* BOUNDARY CONDITIONS TO BE APPLIED IN THE WILSON ACTION                  */
/***************************************************************************/

#ifdef PLAQ_WEIGHTS

#ifdef GAUGE_SPATIAL_TWIST
static void init_plaq_twisted_BCs() {
    error(plaq_weight == NULL, 1, "init_plaq_twisted_BCs [boundary_conditions.c]", "Structure plaq_weight not initialized yet");

    int loc[4] = { T, X, Y, Z };
    int mu, nu, rho, x[4], index;

    for (mu = 1; mu < 3; mu++) {
        for (nu = mu + 1; nu < 4; nu++) {
            rho = 6 - mu - nu;
            x[mu] = 1;
            x[nu] = 1;
            if (COORD[mu] == 0 && COORD[nu] == 0) {
                for (x[0] = 0; x[0] < T; x[0]++) {
                    for (x[rho] = 0; x[rho] < loc[rho]; x[rho]++) {
                        index = ipt(x[0], x[1], x[2], x[3]);
                        plaq_weight[index * 16 + mu * 4 + nu] *= -1.;
                        plaq_weight[index * 16 + nu * 4 + mu] *= -1.; /*IF COMPLEX, THE WEIGHT SHOULD BE C.C.*/
                    }
                }
            }
        }
    }

    lprintf("BCS", 0, "Twisted BCs. Dirac strings intersecting at ( X , Y , Z ) = ( 1 , 1 , 1 )\n");
}
#endif //GAUGE_SPATIAL_TWIST

void init_plaq_open_BCs(double *lplaq_weight, double *lrect_weight, double ct, double cs) {
    error(lplaq_weight == NULL, 1, "init_plaq_open_BCs [boundary_conditions.c]", "Structure plaq_weight not initialized yet");

    int mu, nu, ix, iy, iz, index;

    if (COORD[0] == 0) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            mu = 0;
                            for (nu = mu + 1; nu < 4; nu++) {
                                lplaq_weight[index * 16 + mu * 4 + nu] = 0;
                                lplaq_weight[index * 16 + nu * 4 + mu] = 0;
                                if (lrect_weight != NULL) { lrect_weight[index * 16 + mu * 4 + nu] = 0; }
                                if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0; }
                            }
                            for (mu = 1; mu < 3; mu++) {
                                for (nu = mu + 1; nu < 4; nu++) {
                                    lplaq_weight[index * 16 + mu * 4 + nu] = 0.5 * cs;
                                    lplaq_weight[index * 16 + nu * 4 + mu] = 0.5 * cs;
                                    if (lrect_weight != NULL) { lrect_weight[index * 16 + mu * 4 + nu] = 0.5 * cs; }
                                    if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0.5 * cs; }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        mu = 0;
                        for (nu = mu + 1; nu < 4; nu++) {
                            lplaq_weight[index * 16 + mu * 4 + nu] = ct;
                            lplaq_weight[index * 16 + nu * 4 + mu] = ct;
                        }
                        for (mu = 1; mu < 3; mu++) {
                            for (nu = mu + 1; nu < 4; nu++) {
                                lplaq_weight[index * 16 + mu * 4 + nu] = 0.5 * cs;
                                lplaq_weight[index * 16 + nu * 4 + mu] = 0.5 * cs;
                                if (lrect_weight != NULL) { lrect_weight[index * 16 + mu * 4 + nu] = 0.5 * cs; }
                                if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0.5 * cs; }
                            }
                        }
                    }
                }
            }
        }
    }

    if (COORD[0] == NP_T - 1) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            mu = 0;
                            for (nu = mu + 1; nu < 4; nu++) {
                                lplaq_weight[index * 16 + mu * 4 + nu] = ct;
                                lplaq_weight[index * 16 + nu * 4 + mu] = ct;
                            }
                            for (mu = 1; mu < 3; mu++) {
                                for (nu = mu + 1; nu < 4; nu++) {
                                    lplaq_weight[index * 16 + mu * 4 + nu] = 0.5 * cs;
                                    lplaq_weight[index * 16 + nu * 4 + mu] = 0.5 * cs;
                                    if (lrect_weight != NULL) { lrect_weight[index * 16 + mu * 4 + nu] = 0.5 * cs; }
                                    if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0.5 * cs; }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 2, ix, iy, iz);
                    if (index != -1) {
                        mu = 0;
                        for (nu = mu + 1; nu < 4; nu++) {
                            lplaq_weight[index * 16 + mu * 4 + nu] = ct;
                            lplaq_weight[index * 16 + nu * 4 + mu] = ct;
                            if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0; }
                        }
                    }

                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        mu = 0;
                        for (nu = mu + 1; nu < 4; nu++) {
                            lplaq_weight[index * 16 + mu * 4 + nu] = 0;
                            lplaq_weight[index * 16 + nu * 4 + mu] = 0;
                            if (lrect_weight != NULL) { lrect_weight[index * 16 + mu * 4 + nu] = 0; }
                            if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0; }
                        }
                        for (mu = 1; mu < 3; mu++) {
                            for (nu = mu + 1; nu < 4; nu++) {
                                lplaq_weight[index * 16 + mu * 4 + nu] = 0.5 * cs;
                                lplaq_weight[index * 16 + nu * 4 + mu] = 0.5 * cs;
                                if (lrect_weight != NULL) { lrect_weight[index * 16 + mu * 4 + nu] = 0.5 * cs; }
                                if (lrect_weight != NULL) { lrect_weight[index * 16 + nu * 4 + mu] = 0.5 * cs; }
                            }
                        }
                    }
                }
            }
        }
    }
}

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)

static void init_plaq_SF_BCs(double ct) {
    error(plaq_weight == NULL, 1, "init_plaq_open_BCs [boundary_conditions.c]", "Structure plaq_weight not initialized yet");
    int mu, nu, ix, iy, iz, index;

    if (COORD[0] == 0) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T_BORDER - 1, ix, iy, iz);
                        if (index != -1) {
                            for (mu = 0; mu < 3; mu++) {
                                for (nu = mu + 1; nu < 4; nu++) {
                                    plaq_weight[index * 16 + mu * 4 + nu] = 0.0; //0.5 * cs;
                                    plaq_weight[index * 16 + nu * 4 + mu] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T_BORDER, ix, iy, iz);
                    if (index != -1) {
                        for (mu = 0; mu < 3; mu++) {
                            for (nu = mu + 1; nu < 4; nu++) {
                                plaq_weight[index * 16 + mu * 4 + nu] = 0.0;
                                plaq_weight[index * 16 + nu * 4 + mu] = 0.0;
                            }
                        }
                    }
                    index = ipt_ext(T_BORDER + 1, ix, iy, iz);
                    if (index != -1) {
                        mu = 0;
                        for (nu = mu + 1; nu < 4; nu++) {
                            plaq_weight[index * 16 + mu * 4 + nu] = ct;
                            plaq_weight[index * 16 + nu * 4 + mu] = ct;
                        }

                        for (mu = 1; mu < 3; mu++) {
                            for (nu = mu + 1; nu < 4; nu++) {
                                plaq_weight[index * 16 + mu * 4 + nu] = 0.0;
                                plaq_weight[index * 16 + nu * 4 + mu] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }

    if (COORD[0] == NP_T - 1) {
        if (T_BORDER > 0) {
            for (ix = 0; ix < X_EXT; ++ix) {
                for (iy = 0; iy < Y_EXT; ++iy) {
                    for (iz = 0; iz < Z_EXT; ++iz) {
                        index = ipt_ext(T + T_BORDER, ix, iy, iz);
                        if (index != -1) {
                            for (mu = 0; mu < 3; mu++) {
                                for (nu = mu + 1; nu < 4; nu++) {
                                    plaq_weight[index * 16 + mu * 4 + nu] = 0.0;
                                    plaq_weight[index * 16 + nu * 4 + mu] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(T + T_BORDER - 2, ix, iy, iz);
                    if (index != -1) {
                        mu = 0;
                        for (nu = mu + 1; nu < 4; nu++) {
                            plaq_weight[index * 16 + mu * 4 + nu] = ct;
                            plaq_weight[index * 16 + nu * 4 + mu] = ct;
                        }
                    }

                    index = ipt_ext(T + T_BORDER - 1, ix, iy, iz);
                    if (index != -1) {
                        for (mu = 0; mu < 3; mu++) {
                            for (nu = mu + 1; nu < 4; nu++) {
                                plaq_weight[index * 16 + mu * 4 + nu] = 0.0;
                                plaq_weight[index * 16 + nu * 4 + mu] = 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif
#endif
