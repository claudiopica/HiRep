/* arXiv:1006.4518 [hep-lat] */

#include "utils.h"
#include "libhr_core.h"
#include "memory.h"
#include "IO/logger.h"
#include "inverters.h"
#include "update.h"

static suNg_field *ws_gf = NULL;
static suNg_field *ws_gf_tmp = NULL;
static suNg_field *Vprime = NULL;
static suNg_field *Vtmp = NULL;
static double *wf_plaq_weight = NULL;

void WF_initialize() {
#ifdef GAUGE_SPATIAL_TWIST
    error(0 == 0, 0, "WF_initialize", "WF has not yet been setup to work with GAUGE_SPATIAL_TWIST enabled");
#endif

    if (ws_gf == NULL) {
        ws_gf = alloc_suNg_field(&glattice);
        ws_gf_tmp = alloc_suNg_field(&glattice);
        Vprime = alloc_suNg_field(&glattice);
        Vtmp = alloc_suNg_field(&glattice);

#if defined(PLAQ_WEIGHTS)
#ifdef PURE_GAUGE_ANISOTROPY
        wf_plaq_weight = malloc(sizeof(double) * glattice.gsize_gauge * 16);
        for (int i = 0; i < 16 * glattice.gsize_gauge; i++) {
            wf_plaq_weight[i] = 1.0;
        }
#ifdef BC_T_OPEN
        init_plaq_open_BCs(wf_plaq_weight, NULL, 1.0, 1.0);
#endif

#else
        wf_plaq_weight = plaq_weight;
#endif
#endif
    }
}

void WF_set_bare_anisotropy(double *wf_chi) {
#ifdef PURE_GAUGE_ANISOTROPY
    WF_initialize();
    int ix, iy, iz, it, index, mu, nu;

    for (it = 0; it < T_EXT; ++it) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(it, ix, iy, iz);
                    if (index != -1) {
                        nu = 0;
                        for (mu = nu + 1; mu < 4; mu++) {
                            wf_plaq_weight[index * 16 + nu * 4 + mu] *= *wf_chi * *wf_chi;
                        }
                    }
                }
            }
        }
    }

#else
    error(0 == 0, 0, "WF_set_bare_anisotropy",
          "In order to use anisotropic lattice you must compile with PURE_GAUGE_ANISOTROPY enabled");
#endif
}

void WF_free() {
    if (ws_gf != NULL) {
        free_suNg_field(ws_gf);
        free_suNg_field(ws_gf_tmp);
        free_suNg_field(Vprime);
        free_suNg_field(Vtmp);
        if (wf_plaq_weight != NULL) {
#ifdef PLAQ_WEIGHTS
            if (wf_plaq_weight != plaq_weight) { free(wf_plaq_weight); }
#else
            free(wf_plaq_weight);
#endif
        }
    }
}

/*
 * d/dt V = Z(V) V
 * S_W = 1/g^2 \sum_{p oriented} Re tr ( 1 - V(p) )
 * d_L f(V) = T^a d/ds f(e^{sT^a} V)
 * T_a^dag = -T_a      tr T_a T_b = -delta_{ab}/2
 * Z(V) = -g^2 d_L S_W = \sum_{p oriented} d_L Re tr ( V(p) )
 * Z(V) = d_L Re ( tr V s + tr V^dag s^dag ) =
 *      = d_L ( tr V s + tr V^dag s^dag ) =
 *      = T^a tr T^a ( V s - s^dag V^dag )
 *      = - 1/2 ( V s - s^dag V^dag - 1/N tr (V s - s^dag V^dag) )
 */
static void Zeta(suNg_field *Zu, const suNg_field *U, const double alpha) {
    error(Zu->type != &glattice, 1, "wilson_flow.c", "'Z' in Zeta must be defined on the whole lattice");
    error(U->type != &glattice, 1, "wilson_flow.c", "'U' in Zeta must be defined on the whole lattice");
    suNg staple, tmp1, tmp2;
    int ix, iy, iz, it, i;

    for (it = 0; it < T; ++it) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    i = ipt(it, ix, iy, iz);

                    for (int mu = 0; mu < 4; ++mu) {
                        _suNg_zero(staple);
                        for (int nu = (mu + 1) % 4; nu != mu; nu = (nu + 1) % 4) {
                            suNg *u1 = _4FIELD_AT(U, iup(i, mu), nu);
                            suNg *u2 = _4FIELD_AT(U, iup(i, nu), mu);
                            suNg *u3 = _4FIELD_AT(U, i, nu);
                            _suNg_times_suNg_dagger(tmp1, *u1, *u2);
                            _suNg_times_suNg_dagger(tmp2, tmp1, *u3);
#ifdef PLAQ_WEIGHTS
                            _suNg_mul(tmp2, wf_plaq_weight[i * 16 + nu * 4 + mu], tmp2);
#endif
                            _suNg_add_assign(staple, tmp2);

                            int j = idn(i, nu);
                            u1 = _4FIELD_AT(U, iup(j, mu), nu);
                            u2 = _4FIELD_AT(U, j, mu);
                            u3 = _4FIELD_AT(U, j, nu);
                            _suNg_times_suNg(tmp1, *u2, *u1);
                            _suNg_dagger_times_suNg(tmp2, tmp1, *u3);

#ifdef PLAQ_WEIGHTS
                            _suNg_mul(tmp2, wf_plaq_weight[j * 16 + nu * 4 + mu], tmp2);
#endif
                            _suNg_add_assign(staple, tmp2);
                        }

                        _suNg_times_suNg(tmp1, *_4FIELD_AT(U, i, mu), staple);

                        _suNg_dagger(tmp2, tmp1);
                        _suNg_sub_assign(tmp1, tmp2);
#if !defined(GAUGE_SON) && !defined(WITH_QUATERNIONS)
                        double imtr;
                        _suNg_trace_im(imtr, tmp1);
                        imtr = imtr / NG;
                        for (int k = 0; k < NG * NG; k += NG + 1) {
                            tmp1.c[k] -= I * imtr;
                        }
#endif
#ifdef BC_T_OPEN
                        if (mu != 0 && ((it + zerocoord[0] == 0) || (it + zerocoord[0] == GLB_T - 1))) {
                            _suNg_mul(tmp1, -alpha, tmp1);
                        } else {
                            _suNg_mul(tmp1, -alpha / 2., tmp1);
                        }
#else
                        _suNg_mul(tmp1, -alpha / 2., tmp1);
#endif

                        _suNg_add_assign(*_4FIELD_AT(Zu, i, mu), tmp1);
                    }
                }
            }
        }
    }
}

void WilsonFlow1(suNg_field *V, const double epsilon) {
    zero(ws_gf);
    Zeta(ws_gf, V, epsilon);

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[2];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
            *_4FIELD_AT(V, ix, mu) = utmp[1];
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif
}

//define distance between complex matrices
double max_distance(suNg_field *V, suNg_field *Vprimel) {
    double d, tmp;
    suNg diff;
    _suNg_zero(diff);

    tmp = 0;
    _MASTER_FOR(&glattice, ix) {
        d = 0.;

        for (int mu = 0; mu < 4; ++mu) {
            _suNg_mul(diff, 1., *_4FIELD_AT(V, ix, mu));
            _suNg_sub_assign(diff, *_4FIELD_AT(Vprimel, ix, mu));
            _suNg_sqnorm(d, diff);
            if (d > tmp) { tmp = d; }
        }
    }
    global_max(&tmp, 1);

    return tmp / (double)NG;
}

// following 1301.4388
int WilsonFlow3_adaptative(suNg_field *V, double *epsilon, double *epsilon_new, double *delta) {
    double varepsilon, d;
    copy_suNg_field(Vtmp, V);
    zero(ws_gf);
    zero(ws_gf_tmp);
    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    Zeta(ws_gf, V, *epsilon / 4.); //ws_gf = Z0/4

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[2];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
            *_4FIELD_AT(V, ix, mu) = utmp[1]; // V = exp(Z0/4) W0
            _suNg_mul(*_4FIELD_AT(ws_gf_tmp, ix, mu), -4., *_4FIELD_AT(ws_gf, ix, mu)); //ws_gf_tmp = -Z0
            _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -17. / 9., *_4FIELD_AT(ws_gf, ix, mu)); //ws_gf =  -17*Z0/36
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    Zeta(ws_gf, V, 8. * *epsilon / 9.); // ws_gf = 8 Z1 /9 - 17 Z0/36
    Zeta(ws_gf_tmp, V, 2. * *epsilon); // ws_gf_tmp = 2 Z1 - Z0

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[4];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            suNg_Exp(&utmp[2], _4FIELD_AT(ws_gf_tmp, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu)); // utmp[1] = exp(8 Z1/9 - 17 Z0/36) W1
            _suNg_times_suNg(utmp[3], utmp[2], *_4FIELD_AT(V, ix, mu)); // utmp[4] = exp( Z1 -  Z0) W1
            *_4FIELD_AT(V, ix, mu) = utmp[1];
            *_4FIELD_AT(Vprime, ix, mu) = utmp[3];
            _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -1., *_4FIELD_AT(ws_gf, ix, mu));
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);

    start_sendrecv_suNg_field(Vprime);
    complete_sendrecv_suNg_field(Vprime);

#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    Zeta(ws_gf, V, 3. * *epsilon / 4.);

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[2];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
            *_4FIELD_AT(V, ix, mu) = utmp[1];
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    // now need to get the maximum of the distance
    d = max_distance(V, Vprime);
    varepsilon = *epsilon * pow(*delta / d, 1. / 3.);

    if (d > *delta) {
        copy_suNg_field(V, Vtmp);
        *epsilon_new = 0.5 * varepsilon;
        lprintf("WILSONFLOW", 20, "d > delta : must repeat the calculation with epsilon=%lf\n", *epsilon_new);
        return 1 == 0;
    } else {
        if (varepsilon < 1.0) {
            *epsilon_new = 0.95 * varepsilon;
        } else {
            *epsilon_new = 1.0;
        }

        return 1 == 1;
    }
}

void WilsonFlow3(suNg_field *V, const double epsilon) {
    zero(ws_gf);
    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    Zeta(ws_gf, V, epsilon / 4.);

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[2];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
            *_4FIELD_AT(V, ix, mu) = utmp[1];
            _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -17. / 9., *_4FIELD_AT(ws_gf, ix, mu));
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    Zeta(ws_gf, V, 8. * epsilon / 9.);

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[2];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
            *_4FIELD_AT(V, ix, mu) = utmp[1];
            _suNg_mul(*_4FIELD_AT(ws_gf, ix, mu), -1., *_4FIELD_AT(ws_gf, ix, mu));
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif

    Zeta(ws_gf, V, 3. * epsilon / 4.);

    _MASTER_FOR(&glattice, ix) {
        suNg utmp[2];
        for (int mu = 0; mu < 4; ++mu) {
            suNg_Exp(&utmp[0], _4FIELD_AT(ws_gf, ix, mu));
            _suNg_times_suNg(utmp[1], utmp[0], *_4FIELD_AT(V, ix, mu));
            *_4FIELD_AT(V, ix, mu) = utmp[1];
        }
    }

    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
#if defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    apply_BCs_on_fundamental_gauge_field();
#endif
}

static void WF_measure_and_store(suNg_field *V, storage_switch swc, data_storage_array **ret, int nmeas, int idmeas,
                                 double *t) {
    double TC;
#if defined(BC_T_ANTIPERIODIC) || defined(BC_T_PERIODIC) && !defined(PURE_GAUGE_ANISOTROPY)
    double E_density, Esym_density;
    int idx[2] = { nmeas + 1, 4 };
    if (swc == STORE && *ret == NULL) {
        *ret = allocate_data_storage_array(1);
        allocate_data_storage_element(*ret, 0, 2, idx); // ( nmeas+1 ) * ( 4 <=> t, E, Esym, TC)
    }
#else
    int j;
    double E_density[2 * GLB_T];
    double Esym_density[2 * GLB_T];
    double E_densityavg[2];
    double Esym_densityavg[2];
    int idx[3] = { nmeas + 1, GLB_T, 5 };
    if (swc == STORE && *ret == NULL) {
        *ret = allocate_data_storage_array(2);
        allocate_data_storage_element(
            *ret, 0, 3,
            idx); // ( nmeas+1 ) * ( GLB_T ) * ( 5 <=> t, Etime, Espace, Esym_densitytime, Esym_densityspace)
        idx[1] = 6;
        allocate_data_storage_element(
            *ret, 1, 2,
            idx); // ( nmeas+1 ) * ( 6 <=> t, Etime_tsum, Espace_tsum, Esym_densitytime_tsum, Esym_densityspace_tsum, TC)
    }
#endif

    TC = topo(V);

#if defined(BC_T_ANTIPERIODIC) || defined(BC_T_PERIODIC) && !defined(PURE_GAUGE_ANISOTROPY)

    E_density = E(V);
    Esym_density = Esym(V);
    lprintf("WILSONFLOW", 0, "WF (t,E,t2*E,Esym,t2*Esym,TC) = %1.8e %1.16e %1.16e %1.16e %1.16e %1.16e\n", *t, E_density,
            *t * *t * E_density, Esym_density, *t * *t * Esym_density, TC);
    if (swc == STORE) {
        idx[0] = idmeas - 1;
        idx[1] = 0;
        *data_storage_element(*ret, 0, idx) = *t;
        idx[1] = 1;
        *data_storage_element(*ret, 0, idx) = E_density;
        idx[1] = 2;
        *data_storage_element(*ret, 0, idx) = Esym_density;
        idx[1] = 3;
        *data_storage_element(*ret, 0, idx) = TC;
    }
#else

    E_T(E_density, V);
    Esym_T(Esym_density, V);

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
    E_density[0] = E_density[1] = Esym_density[0] = Esym_density[1] = Esym_density[2] = Esym_density[3] = 0.0;
    E_density[2 * GLB_T - 2] = E_density[2 * GLB_T - 1] = Esym_density[2 * GLB_T - 2] = Esym_density[2 * GLB_T - 1] = 0.0;
#elif defined(BC_T_OPEN)
    E_density[2 * (GLB_T - 1)] = 0.0;
    Esym_density[0] = Esym_density[1] = Esym_density[2 * GLB_T - 2] = Esym_density[2 * GLB_T - 1] = 0.0;
#endif

    if (swc == STORE) {
        for (j = 0; j < GLB_T; j++) {
            idx[0] = idmeas - 1;
            idx[1] = j;
            idx[2] = 0;
            *data_storage_element(*ret, 0, idx) = *t;
            idx[2] = 1;
            *data_storage_element(*ret, 0, idx) = E_density[2 * j];
            idx[2] = 2;
            *data_storage_element(*ret, 0, idx) = E_density[2 * j + 1];
            idx[2] = 3;
            *data_storage_element(*ret, 0, idx) = Esym_density[2 * j];
            idx[2] = 4;
            *data_storage_element(*ret, 0, idx) = Esym_density[2 * j + 1];
        }
    }

    E_densityavg[0] = E_densityavg[1] = Esym_densityavg[0] = Esym_densityavg[1] = 0.0;
    for (j = 0; j < GLB_T; j++) {
        lprintf("WILSONFLOW", 0, "WF (T,t,Etime,Espace,Esymtime,Esymspace) = %d %1.8e %1.16e %1.16e %1.16e %1.16e\n", j, *t,
                E_density[2 * j], E_density[2 * j + 1], Esym_density[2 * j], Esym_density[2 * j + 1]);
        E_densityavg[0] += E_density[2 * j];
        E_densityavg[1] += E_density[2 * j + 1];
        Esym_densityavg[0] += Esym_density[2 * j];
        Esym_densityavg[1] += Esym_density[2 * j + 1];
    }

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
    E_densityavg[0] /= GLB_T - 2;
    E_densityavg[1] /= GLB_T - 3;
    Esym_densityavg[0] /= GLB_T - 3;
    Esym_densityavg[1] /= GLB_T - 3;
#else
    E_densityavg[0] /= GLB_T;
    E_densityavg[1] /= GLB_T;
    Esym_densityavg[0] /= GLB_T;
    Esym_densityavg[1] /= GLB_T;
#endif

    lprintf(
        "WILSONFLOW", 0,
        "WF avg (t,Etime,Espace,Esymtime,Esymspace,Pltime,Plspace,TC) = %1.8e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
        *t, E_densityavg[0], E_densityavg[1], Esym_densityavg[0], Esym_densityavg[1], (NG - E_densityavg[0]),
        (NG - E_densityavg[1]), TC);
    if (swc == STORE) {
        idx[0] = idmeas - 1;
        idx[1] = 0;
        *data_storage_element(*ret, 1, idx) = *t;
        idx[1] = 1;
        *data_storage_element(*ret, 1, idx) = E_densityavg[0];
        idx[1] = 2;
        *data_storage_element(*ret, 1, idx) = E_densityavg[1];
        idx[1] = 3;
        *data_storage_element(*ret, 1, idx) = Esym_densityavg[0];
        idx[1] = 4;
        *data_storage_element(*ret, 1, idx) = Esym_densityavg[1];
        idx[1] = 5;
        *data_storage_element(*ret, 1, idx) = TC;
    }
#endif
}

data_storage_array *WF_update_and_measure(WF_integrator_type wft, suNg_field *V, double *tmax, double *eps, double *delta,
                                          int nmeas, storage_switch swc) {
    data_storage_array *ret = NULL;

    int k = 1;
    double t = 0.;
    double epsilon = *eps;
    double dt = *tmax / nmeas;
    double epsilon_new = *eps;
    double epsilon_meas = 0;

#ifdef WITH_GPU
    error(1, 1, "WF_update_and_measure [wilsonflow.c]", "The wilsonflow is not implemented on GPU. Compile without WITH_GPU.");
#endif

    WF_measure_and_store(V, swc, &ret, nmeas, k, &t);

    while (t < *tmax) {
        if (t + epsilon > (double)k * dt) { epsilon = (double)k * dt - t; }

        switch (wft) {
        case EUL:
            WilsonFlow1(V, epsilon);
            t += epsilon;
            break;

        case RK3:
            WilsonFlow3(V, epsilon);
            t += epsilon;
            break;

        case RK3_ADAPTIVE:
            while (!WilsonFlow3_adaptative(V, &epsilon, &epsilon_new, delta)) {
                epsilon = epsilon_new;
            }

            t += epsilon;
            epsilon = epsilon_new;
            break;
        }

        if ((double)k * dt - t < epsilon) {
            if ((double)k * dt - t < 1.e-10) {
                k++;

                WF_measure_and_store(V, swc, &ret, nmeas, k, &t);
                if (epsilon_meas > epsilon) { epsilon = epsilon_meas; }
                if (epsilon > dt) { epsilon = dt; }
            } else {
                epsilon_meas = epsilon;
                epsilon = (double)k * dt - t;
            }
        }
    }
    return ret;
}