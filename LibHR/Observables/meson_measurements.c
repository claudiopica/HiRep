/*******************************************************************************
 *                                                                              *
 * Wrapper functions for different type of measurements                         *
 * Copyright (c) 2013 Rudy Arthur, Ari Hietanen                                 *
 *                                                                              *
 *******************************************************************************/
#include "observables.h"
#include "libhr_core.h"
#include "io.h"
#include "memory.h"
#include "utils.h"
#include "update.h"
#include "random.h"
#include "inverters.h"

static void fix_T_bc(int tau) {
    int index;
    int ix, iy, iz;
    suNf *u;
    if (--tau < 0) { tau += GLB_T; }
    lprintf("meson_measurements", 15, "Setting Dirichlet boundary conidtion at global time slice %d (local %d)\n", tau,
            tau - zerocoord[0]);
    if ((zerocoord[0] - 1 <= tau && zerocoord[0] + T > tau) || (zerocoord[0] == 0 && tau == GLB_T - 1)) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    if (tau == zerocoord[0] - 1 || (zerocoord[0] == 0 && tau == GLB_T - 1)) {
                        index = ipt_ext(0, ix, iy, iz);
                    } else {
                        index = ipt_ext(T_BORDER + tau - zerocoord[0], ix, iy, iz);
                    }
                    if (index != -1) {
                        u = pu_gauge_f(index, 0);
                        _suNf_zero(*u);
                    }
                }
            }
        }
    }
    lprintf("meson_measurements", 50, "Boundaries set!\n");
}

static void flip_T_bc(int tau) {
    int index;
    int ix, iy, iz;
    suNf *u;
    tau -= 1;
    if (tau < 0) { tau += GLB_T; }
    lprintf("meson_measurements", 15, "Flipping the boundary at global time slice %d\n", tau);
    fflush(stdout);
    if ((zerocoord[0] - 1 <= tau && zerocoord[0] + T > tau) || (zerocoord[0] == 0 && tau == GLB_T - 1)) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    if ((tau == zerocoord[0] - 1) || (zerocoord[0] == 0 && tau == GLB_T - 1)) {
                        index = ipt_ext(0, ix, iy, iz);
                    } else {
                        index = ipt_ext(T_BORDER + tau - zerocoord[0], ix, iy, iz);
                    }
                    if (index != -1) {
                        u = pu_gauge_f(index, 0);
                        _suNf_minus(*u, *u);
                    }
                }
            }
        }
    }
    lprintf("meson_measurements", 50, "Flipping DONE!\n");
}

/********************************
 *	Point Sources		*
 *********************************/

#define corr_ind(px, py, pz, n_mom, tc, nm, cm) \
    ((px) * (n_mom) * (n_mom) * (24) * (nm) + (py) * (n_mom) * (24) * (nm) + (pz) * (24) * (nm) + ((cm) * (24)) + (tc))

void measure_spectrum_pt(int tau, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                         data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm * NF, &glattice);
    for (int i = 0; i < 4 * nm * NF; i++) {
        zero_spinor_field(prop + i);
    }

    // init data storage here
    if (swc == STORE) {
        int idx[5] = { nm, pow(n_mom, 3), 16, GLB_T, 2 };
        *ret = allocate_data_storage_array(1);
        allocate_data_storage_element(*ret, 0, 5, idx); // ( 1 ) * (nmom^3*ngamma*GLB_T * 2 reals )
        lprintf("MAIN", 0, "data_storage_element allocated !\n");
    }

    init_propagator_eo(nm, m, precision);
    int k;
    lprintf("MAIN", 0, "Point Source at (%d,0,0,0) \n", tau);
    for (k = 0; k < NF; ++k) {
        create_point_source(source, tau, k);
        calc_propagator(prop + 4 * k, source, 4); // 4 for spin components
        if (n_mom > 1) {
            measure_point_mesons_momenta(meson_correlators, prop + 4 * k, source, nm, tau, n_mom);
        } else {
            measure_mesons(meson_correlators, prop + 4 * k, source, nm, tau);
        }
    }
    measure_conserved_currents(cvc_correlators, prop, source, nm, tau);

    if (swc == STORE) {
        double norm = -(1. / GLB_VOL3);
        meson_observable *motmp = meson_correlators;

        int iG = 0;
        while (motmp != NULL) {
            global_sum(motmp->corr_re, motmp->corr_size);
            global_sum(motmp->corr_im, motmp->corr_size);
            for (int i = 0; i < motmp->corr_size; i++) {
                motmp->corr_re[i] *= norm;
                motmp->corr_im[i] *= norm;
            }
            int ip = 0;
            for (int px = 0; px < n_mom; ++px) {
                for (int py = 0; py < n_mom; ++py) {
                    for (int pz = 0; pz < n_mom; ++pz) {
                        for (int im = 0; im < nm; im++) {
                            if (motmp->ind1 == motmp->ind2) {
                                for (int t = 0; t < GLB_T; ++t) {
                                    lprintf("MAIN", 0, "(px,py,pz) = (%d,%d,%d), im = %d, ip =%d, iG=%d, t=%d  %3.10e\n", px,
                                            py, pz, im, ip, iG, t, motmp->corr_re[corr_ind(px, py, pz, n_mom, t, nm, im)]);
                                    int idx[5] = { im, ip, iG, t, 0 };
                                    *data_storage_element(*ret, 0, idx) =
                                        motmp->corr_re[corr_ind(px, py, pz, n_mom, t, nm, im)];
                                    idx[4] = 1;
                                    *data_storage_element(*ret, 0, idx) =
                                        motmp->corr_im[corr_ind(px, py, pz, n_mom, t, nm, im)];
                                }
                            }
                        }
                        ip += 1;
                    }
                }
            }
            iG += 1;
            motmp = motmp->next;
        }
    }

    print_mesons(meson_correlators, 1., conf_num, nm, m, GLB_T, n_mom, "DEFAULT_POINT");
    print_mesons(cvc_correlators, 1., conf_num, nm, m, GLB_T, n_mom, "DEFAULT_POINT");

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
}

void measure_spectrum_pt_ext(int tau, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                             data_storage_array **ret) {
    int k, l;
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop_p = alloc_spinor_field(8 * nm, &glattice);
    spinor_field *prop_a = prop_p + 4 * nm;

    for (int i = 0; i < 8 * nm; i++) {
        zero_spinor_field(prop_p + i);
    }

    init_propagator_eo(nm, m, precision);
    lprintf("MAIN", 10, "Point Source at (%d,0,0,0) \n", tau);
    for (k = 0; k < NF; ++k) {
        create_point_source(source, tau, k);
        calc_propagator(prop_p, source, 4); // 4 for spin components
        if (n_mom > 0) {
            measure_point_mesons_momenta_ext(meson_correlators, prop_p, source, nm, tau, n_mom, 0);
        } else {
            measure_mesons_ext(meson_correlators, prop_p, source, nm, tau, 0);
        }
        flip_T_bc(tau);
        calc_propagator(prop_a, source, 4); // 4 for spin components
        flip_T_bc(tau);
        for (l = 0; l < 4 * nm; ++l) {
            add_assign_spinor_field(&prop_p[l], &prop_a[l]);
            mul_spinor_field(&prop_p[l], 0.5, &prop_p[l]);
        }
        if (n_mom > 1) {
            measure_point_mesons_momenta_ext(meson_correlators, prop_p, source, nm, tau, n_mom, 1);
            for (l = 0; l < 4 * nm; ++l) {
                mul_add_assign_spinor_field(&prop_p[l], -1., &prop_a[l]);
            }
            measure_point_mesons_momenta_ext(meson_correlators, prop_p, source, nm, tau, n_mom, 2);
        } else {
            measure_mesons_ext(meson_correlators, prop_p, source, nm, tau, 1);
            for (l = 0; l < 4 * nm; ++l) {
                mul_add_assign_spinor_field(&prop_p[l], -1., &prop_a[l]);
            }
            measure_mesons_ext(meson_correlators, prop_p, source, nm, tau, 2);
        }
    }
    print_mesons(meson_correlators, 1., conf_num, nm, m, 3 * GLB_T, n_mom, "EXTENDED_POINT");
    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop_p);
    free_spinor_field(prop_a);
}

void measure_spectrum_pt_fixedbc(int tau, int dt, int nm, double *m, int n_mom, int conf_num, double precision,
                                 storage_switch swc, data_storage_array **ret) {
    int k;
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);
    suNf_field *u_gauge_old = alloc_suNf_field(&glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }
    char label[100];
    copy_suNf_field(u_gauge_old, u_gauge_f);
    init_propagator_eo(nm, m, precision);
    fix_T_bc(tau - dt); // Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.
    lprintf("MAIN", 0, "Point Source at (%d,0,0,0) \n", tau);
    for (k = 0; k < NF; ++k) {
        create_point_source(source, tau, k);
        calc_propagator(prop, source, 4); // 4 for spin components
        if (n_mom > 0) {
            measure_point_mesons_momenta(meson_correlators, prop, source, nm, tau, n_mom);
        } else {
            measure_mesons(meson_correlators, prop, source, nm, tau);
        }
    }
    sprintf(label, "DIRICHLET_POINT dt=%d", dt);
    print_mesons(meson_correlators, 1., conf_num, nm, m, GLB_T, n_mom, label);
    copy_suNf_field(u_gauge_f, u_gauge_old);
    free_spinor_field(source);
    free_spinor_field(prop);
    free_suNf_field(u_gauge_old);
    free_propagator_eo();
}

/********************************
 *	SEMWall Sources		*
 *********************************/

void measure_diquark_semwall_background(int nm, double *m, int nhits, int conf_num, double precision, double Q, int n,
                                        storage_switch swc, data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glat_even);
    spinor_field *prop_u = alloc_spinor_field(4 * nm, &glattice);
    spinor_field *prop_d = alloc_spinor_field(4 * nm, &glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop_u + i);
    }
    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop_d + i);
    }

    int tau, k;
    suNg_field *u_gauge_old = alloc_suNg_field(&glattice);
    copy_suNg_field(u_gauge_old, u_gauge);

    error(nm != 1, 1, "measure_diquark_semwall_background", "nm cannot be different from 1 !\n");

    init_propagator_eo(nm, m, precision);

    for (k = 0; k < nhits; ++k) {
        tau = create_diluted_source_equal_eo(source);
        // apply background and calculate first prop
        apply_background_field_zdir(u_gauge, Q, n);
        represent_gauge_field();
        calc_propagator_eo(prop_u, source, 4); // 4 for spin dilution
            // apply background and calculate second prop
        copy_suNg_field(u_gauge, u_gauge_old);

        apply_background_field_zdir(u_gauge, -Q, n);
        represent_gauge_field();
        calc_propagator_eo(prop_d, source, 4); // 4 for spin dilution

        measure_diquarks(meson_correlators, prop_u, prop_d, source, nm, tau);

        copy_suNg_field(u_gauge, u_gauge_old);
        represent_gauge_field();
    }
    print_mesons(meson_correlators, nhits * GLB_VOL3 / 2., conf_num, nm, m, GLB_T, 1, "DEFAULT_DIQUARK_SEMWALL_BACKGROUND");

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop_u);
    free_spinor_field(prop_d);
    free_suNg_field(u_gauge_old);
}

void measure_spectrum_semwall(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                              data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);

    // init data storage here
    if (swc == STORE) {
        int idx[4] = { nm, 16, GLB_T, 2 };
        *ret = allocate_data_storage_array(1);
        allocate_data_storage_element(*ret, 0, 4, idx); // ( 1 ) * (nmom^3*ngamma*GLB_T * 2 reals )
        lprintf("MAIN", 0, "data_storage_element allocated !\n");
    }
    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }

    int tau, k;
    init_propagator_eo(nm, m, precision);
    for (k = 0; k < nhits; ++k) {
        tau = create_diluted_source_equal_eo(source);

        calc_propagator_eo(prop, source, 4); // 4 for spin dilution

        measure_mesons(meson_correlators, prop, source, nm, tau);
    }

    if (swc == STORE) {
        double norm = -(1. / (nhits * GLB_VOL3 / 2.)) / GLB_VOL3;
        meson_observable *motmp = meson_correlators;

        int iG = 0;
        while (motmp != NULL) {
            global_sum(motmp->corr_re, motmp->corr_size);
            global_sum(motmp->corr_im, motmp->corr_size);
            for (int i = 0; i < motmp->corr_size; i++) {
                motmp->corr_re[i] *= norm;
                motmp->corr_im[i] *= norm;
            }
            for (int im = 0; im < nm; im++) {
                if (motmp->ind1 == motmp->ind2) {
                    for (int t = 0; t < GLB_T; ++t) {
                        lprintf("MAIN", 0, " im = %d, iG=%d, t=%d  %3.10e\n", im, iG, t,
                                motmp->corr_re[corr_ind(0, 0, 0, 1, t, nm, im)]);
                        int idx[4] = { im, iG, t, 0 };
                        *data_storage_element(*ret, 0, idx) = motmp->corr_re[corr_ind(0, 0, 0, 1, t, nm, im)];
                        idx[3] = 1;
                        *data_storage_element(*ret, 0, idx) = motmp->corr_im[corr_ind(0, 0, 0, 1, t, nm, im)];
                    }
                }
            }
            iG += 1;
            motmp = motmp->next;
        }
    }

    print_mesons(meson_correlators, nhits * GLB_VOL3 / 2., conf_num, nm, m, GLB_T, 1, "DEFAULT_SEMWALL");
    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
}

void measure_spectrum_semwall_ext(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                  data_storage_array **ret) {
    int k, l, tau;
    spinor_field *source = alloc_spinor_field(4, &glat_even);
    spinor_field *prop_p = alloc_spinor_field(8 * nm, &glattice);
    spinor_field *prop_a = prop_p + 4 * nm;

    for (int i = 0; i < 8 * nm; i++) {
        zero_spinor_field(prop_p + i);
    }

    int dilution = 4; // 4 for spin dilution
    init_propagator_eo(nm, m, precision);
    for (k = 0; k < nhits; ++k) {
        tau = create_diluted_source_equal_eo(source);
        lprintf("MEASURE_SPECTRUM_SEMWALL_EXT", 10, "SEM wall source (noise) at time slice %d.\n", tau);
        calc_propagator_eo(prop_p, source, dilution);
        measure_mesons_ext(meson_correlators, prop_p, source, nm, tau, 0);
        flip_T_bc(tau);
        calc_propagator_eo(prop_a, source, dilution);
        flip_T_bc(tau);
        for (l = 0; l < dilution * nm; ++l) {
            add_assign_spinor_field(&prop_p[l], &prop_a[l]);
            mul_spinor_field(&prop_p[l], 0.5, &prop_p[l]);
        }
        measure_mesons_ext(meson_correlators, prop_p, source, nm, tau, 1);
        for (l = 0; l < dilution * nm; ++l) {
            sub_assign_spinor_field(&prop_p[l], &prop_a[l]);
        }
        measure_mesons_ext(meson_correlators, prop_p, source, nm, tau, 2);
    }
    print_mesons(meson_correlators, 1. * nhits * GLB_VOL3 / 2., conf_num, nm, m, 3 * GLB_T, 1, "EXTENDED_SEMWALL");
    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop_p);
}

void measure_spectrum_semwall_fixedbc(int dt, int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                      data_storage_array **ret) {
    int tau, k;
    spinor_field *source = alloc_spinor_field(4, &glat_even);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);
    suNf_field *u_gauge_old = alloc_suNf_field(&glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }

    char label[100];
    copy_suNf_field(u_gauge_old, u_gauge_f);
    init_propagator_eo(nm, m, precision);
    for (k = 0; k < nhits; ++k) {
        tau = create_diluted_source_equal_eo(source);
        fix_T_bc(tau - dt);
        calc_propagator_eo(prop, source, 4); // 4 for spin dilution
        measure_mesons(meson_correlators, prop, source, nm, tau);
        copy_suNf_field(u_gauge_f, u_gauge_old);
    }
    sprintf(label, "DIRICHLET_SEMWALL dt=%d", dt);
    print_mesons(meson_correlators, nhits * GLB_VOL3 / 2., conf_num, nm, m, GLB_T, 1, label);
    free_spinor_field(source);
    free_spinor_field(prop);
    free_suNf_field(u_gauge_old);
    free_propagator_eo();
}

/****************************************
 *	Gauge Fixed Wall Sources	*
 *****************************************/
void measure_spectrum_gfwall(int nm, double *m, int conf_num, double precision, storage_switch swc, data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);
    suNg_field *u_gauge_old = alloc_suNg_field(&glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }
    int tau, k;
    tau = 0;
    copy_suNg_field(u_gauge_old, u_gauge);
    // Fix the Gauge
    double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
                          1.8, // overrelax
                          10000, // maxit
                          1e-12, // tolerance
                          u_gauge // gauge
    );
    lprintf("GFWALL", 0, "Gauge fixed action  %1.6f\n", act);
    double p2 = calc_plaq(u_gauge);
    lprintf("TEST", 0, "fixed_gauge plaq %1.6f\n", p2);
    full_plaquette();
    represent_gauge_field();

    init_propagator_eo(nm, m, precision);
    for (k = 0; k < NF; ++k) {
        create_gauge_fixed_wall_source(source, tau, k);
        calc_propagator(prop, source, 4); // 4 for spin dilution
        measure_mesons(meson_correlators, prop, source, nm, tau);
    }
    print_mesons(meson_correlators, GLB_VOL3, conf_num, nm, m, GLB_T, 1, "DEFAULT_GFWALL");

    copy_suNg_field(u_gauge, u_gauge_old);
    represent_gauge_field();

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
    free_suNg_field(u_gauge_old);
}

void measure_spectrum_gfwall_fixedbc(int dt, int nm, double *m, int conf_num, double precision, storage_switch swc,
                                     data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);
    suNg_field *u_gauge_old = alloc_suNg_field(&glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }
    int tau, k;
    tau = 0;
    copy_suNg_field(u_gauge_old, u_gauge);

    // Fix the Gauge
    double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
                          1.8, // overrelax
                          10000, // maxit
                          1e-12, // tolerance
                          u_gauge // gauge
    );
    lprintf("GFWALL", 0, "Gauge fixed action  %1.6f\n", act);
    double p2 = calc_plaq(u_gauge);
    lprintf("TEST", 0, "fixed_gauge plaq %1.6f\n", p2);
    full_plaquette();
    represent_gauge_field();

    fix_T_bc(tau - dt); // Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.

    init_propagator_eo(nm, m, precision);
    for (k = 0; k < NF; ++k) {
        create_gauge_fixed_wall_source(source, tau, k);
        calc_propagator(prop, source, 4); // 4 for spin dilution
        measure_mesons(meson_correlators, prop, source, nm, tau);
    }
    print_mesons(meson_correlators, GLB_VOL3, conf_num, nm, m, GLB_T, 1, "DIRICHLET_GFWALL");

    copy_suNg_field(u_gauge, u_gauge_old);
    represent_gauge_field();

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
    free_suNg_field(u_gauge_old);
}

/****************************************
 *	Disconnected Measurements	*
 *****************************************/

void measure_spectrum_discon_semwall(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                     data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }
    int k, beta;
    char label[100];
    init_propagator_eo(nm, m, precision);
    for (k = 0; k < nhits; ++k) {
        create_noise_source_equal_eo(source);
        for (beta = 0; beta < 4; beta++) {
            source[beta].type = &glat_even;
        }
        calc_propagator(prop, source, 4); // 4 for spin dilution
        for (beta = 0; beta < 4; beta++) {
            source[beta].type = &glattice;
        }
        measure_mesons(discon_correlators, prop, source, nm, 0);
        sprintf(label, "src %d DISCON_SEMWALL", k);
        print_mesons(discon_correlators, GLB_VOL3 / 2., conf_num, nm, m, GLB_T, 1, label);
    }
    // print_mesons(discon_correlators,nhits*GLB_VOL3/2.,conf_num,nm,m,GLB_T,1,"DISCON_SEMWALL");
    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
}

void measure_spectrum_discon_gfwall(int nm, double *m, int conf_num, double precision, storage_switch swc,
                                    data_storage_array **ret) {
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);
    suNg_field *u_gauge_old = alloc_suNg_field(&glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }

    int tau, k;
    tau = 0;

    copy_suNg_field(u_gauge_old, u_gauge);
    // Fix the Gauge
    double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
                          1.8, // overrelax
                          10000, // maxit
                          1e-12, // tolerance
                          u_gauge // gauge
    );
    lprintf("GFWALL", 0, "Gauge fixed action  %1.6f\n", act);
    double p2 = calc_plaq(u_gauge);
    lprintf("TEST", 0, "fixed_gauge plaq %1.6f\n", p2);
    full_plaquette();
    represent_gauge_field();
    init_propagator_eo(nm, m, precision);

    for (tau = 0; tau < GLB_T; ++tau) {
        for (k = 0; k < NF; ++k) {
            create_gauge_fixed_wall_source(source, tau, k);
            calc_propagator(prop, source, 4); // 4 for spin dilution
            create_point_source(source, tau, k); // to get the contraction right
            measure_mesons(discon_correlators, prop, source, nm, 0);
        }
    }
    print_mesons(discon_correlators, GLB_T, conf_num, nm, m, GLB_T, 1, "DISCON_GFWALL");

    copy_suNg_field(u_gauge, u_gauge_old);
    represent_gauge_field();

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
    free_suNg_field(u_gauge_old);
}

void measure_spectrum_discon_volume(int nm, double *m, int conf_num, double precision, int dil, storage_switch swc,
                                    data_storage_array **ret) {
    // Spin diluted
    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4 * nm, &glattice);

    for (int i = 0; i < 4 * nm; i++) {
        zero_spinor_field(prop + i);
    }

    init_propagator_eo(nm, m, precision);
    int p;
    for (p = 0; p < dil; p++) {
        create_diluted_volume_source(source, p, dil);
        calc_propagator(prop, source, 4); // spin dilution
        measure_mesons(discon_correlators, prop, source, nm, 0);
    }
    print_mesons(discon_correlators, 1., conf_num, nm, m, GLB_T, 1, "DISCON_VOLUME");

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
}

/****************************************
 *	Form Factor Measurements	*
 *****************************************/
void measure_formfactor_pt(int ti, int tf, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                           data_storage_array **ret) {
    spinor_field *source;
    spinor_field *source_seq;
    spinor_field *prop_i;
    spinor_field *prop_seq;
    int k;

    source = alloc_spinor_field(4, &glattice);
    prop_i = alloc_spinor_field(4 * NF, &glattice);
    prop_seq = alloc_spinor_field(4 * NF, &glattice);
    source_seq = alloc_spinor_field(4 * NF, &glattice);

    for (int i = 0; i < 4 * NF; i++) {
        zero_spinor_field(prop_i + i);
        zero_spinor_field(prop_seq + i);
    }

    init_propagator_eo(1, m, precision); // 1 for number of masses
    int pt[4];
    generate_random_point(pt);
    pt[0] = pt[1] = pt[2] = pt[3] = 0;
    for (k = 0; k < NF; ++k) {
        // create_point_source(source,ti,k);
        create_point_source_loc(source, ti, pt[1], pt[2], pt[3], k);
        calc_propagator(prop_i + 4 * k, source, 4); // 4 for spin components
    }
    create_sequential_source(source_seq, tf, prop_i);
    calc_propagator(prop_seq, source_seq, 4 * NF);

    measure_formfactors(prop_seq, prop_i, source_seq, nm, ti, tf, n_mom, pt); // eats two propagators
    print_formfactor(conf_num, nm, m, n_mom, "DEFAULT_FF_POINT", tf - ti);
    free_spinor_field(source);
    free_spinor_field(source_seq);
    free_spinor_field(prop_i);
    free_spinor_field(prop_seq);
    free_propagator_eo();
}

void measure_formfactor_fixed(int ti, int tf, int dt, int nm, double *m, int n_mom, int conf_num, double precision,
                              storage_switch swc, data_storage_array **ret) {
    spinor_field *source;
    spinor_field *prop_i;
    spinor_field *source_seq;
    spinor_field *prop_seq;

    int k;
    char label[100];
    suNf_field *u_gauge_old = alloc_suNf_field(&glattice);
    copy_suNf_field(u_gauge_old, u_gauge_f); // Save the gaugefield

    double p = avr_plaquette();
    lprintf("MESON_MEASUREMENTS", 0, "<P> = %g\n", p);
    fix_T_bc(ti - dt); // Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.
    p = avr_plaquette();
    lprintf("MESON_MEASUREMENTS", 0, "<P> = %g\n", p);

    source = alloc_spinor_field(4, &glattice);
    source_seq = alloc_spinor_field(4 * NF, &glattice);
    prop_i = alloc_spinor_field(4 * NF, &glattice);
    prop_seq = alloc_spinor_field(4 * NF, &glattice);

    for (int i = 0; i < 4 * NF; i++) {
        zero_spinor_field(prop_i + i);
        zero_spinor_field(prop_seq + i);
    }

    init_propagator_eo(1, m, precision); // 1 for number of masses
    int pt[4];
    generate_random_point(pt);
    // for(k=0;k<4;k++) pt[k] = 0;
    lprintf("MESON_MEASUREMENTS", 0, "Source at (%d,%d,%d,%d)\n", ti, pt[1], pt[2], pt[3]);
    for (k = 0; k < NF; ++k) {
        // create_point_source(source,ti,k);
        create_point_source_loc(source, ti, pt[1], pt[2], pt[3], k);
        calc_propagator(prop_i + 4 * k, source, 4); // 4 for spin components
    }
    create_sequential_source(source_seq, tf, prop_i); // prop_i = S(x,0);
    calc_propagator(prop_seq, source_seq, 4 * NF); // prop_seq = S(y,x) S(x,0) delta(x, (2,0,0,0) )

    measure_formfactors(prop_seq, prop_i, source, nm, ti, tf, n_mom, pt); // eats two propagators

    sprintf(label, "DIRICHLET_FF_POINT dt=%d", dt);
    print_formfactor(conf_num, nm, m, n_mom, label, tf - ti);
    copy_suNf_field(u_gauge_f, u_gauge_old); // Restore the gaugefield

    free_spinor_field(source);
    free_spinor_field(source_seq);
    free_spinor_field(prop_i);
    free_spinor_field(prop_seq);

    free_suNf_field(u_gauge_old);
    free_propagator_eo();
}
