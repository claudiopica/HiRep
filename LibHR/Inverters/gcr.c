/***************************************************************************\
* Copyright (c) 2023 Sofie Martins                                          *
* All rights reserved.                                                      *
\***************************************************************************/

// Written after the following paper:
// ===================================
// Solution of the Dirac equation in lattice QCD using a domain decomposition method
// Martin Luscher (CERN)
// e-Print: hep-lat/0310048 [hep-lat]
// DOI: 10.1016/S0010-4655(03)00486-7
// Published in: Comput.Phys.Commun. 156 (2004), 209-220

#include "inverters.h"
#include "io.h"
#include "memory.h"
#include "utils.h"
#include <assert.h>
#include "update.h"

static spinor_field *rho, *chi, *xi, *tmp;
static int nmr = 4;
static int ncy = 5;
static int maxit_gcr = 15;
static int maxit_sapgcr = 27;

void set_nmr(int newval) {
    nmr = newval;
}

void set_ncy(int newval) {
    ncy = newval;
}

void set_maxit_gcr(int newval) {
    maxit_gcr = newval;
}

void set_maxit_sapgcr(int newval) {
    maxit_sapgcr = newval;
}

hr_complex *a, *b, *c, *alpha;

int sapgcr(mshift_par *par, spinor_operator M, spinor_field *eta, spinor_field *phi) {
    int iter = 0;

    a = (hr_complex *)malloc((maxit_sapgcr + 1) * (maxit_sapgcr + 1) * sizeof(hr_complex));
    b = (hr_complex *)malloc((maxit_sapgcr + 1) * sizeof(hr_complex));
    c = (hr_complex *)malloc((maxit_sapgcr + 1) * sizeof(hr_complex));
    alpha = (hr_complex *)malloc((maxit_sapgcr + 1) * sizeof(hr_complex));

    rho = alloc_spinor_field(4 + 2 * maxit_sapgcr, eta->type);
    tmp = rho + 1;
    chi = tmp + 1;
    xi = chi + maxit_sapgcr + 1;

    double const innorm2 = sqnorm(eta);

    zero(phi);
    copy(rho, eta);

    do {
        int k;
        for (k = 0; k <= maxit_sapgcr; k++) {
            SAP_prec(nmr, ncy, &BiCGstab, par, M, rho, &xi[k]);
            M(&chi[k], &xi[k]);

            for (int l = 0; l < k; l++) {
                a[l + (maxit_sapgcr + 1) * k] = prod(&chi[l], &chi[k]);
                mulc_add_assign(&chi[k], -a[l + (maxit_sapgcr + 1) * k], &chi[l]);
            }

            b[k] = sqrt(sqnorm(&chi[k]));
            mul(&chi[k], 1. / b[k], &chi[k]);
            c[k] = prod(&chi[k], rho);
            mulc_add_assign(rho, -c[k], &chi[k]);

            if (sqnorm(rho) < par->err2 * innorm2) { break; }

            iter++;
        }

        if (k > maxit_sapgcr) { k--; }

        for (int l = k; l >= 0; l--) {
            hr_complex sum = 0;
            for (int i = l + 1; i <= k; i++) {
                sum += a[l + (maxit_sapgcr + 1) * i] * alpha[i];
            }
            alpha[l] = (c[l] - sum) / b[l];
        }

        for (int l = 0; l <= k; l++) {
            mulc_add_assign(phi, alpha[l], &xi[l]);
        }

        if (k < maxit_sapgcr) { break; }

        M(tmp, phi);
        sub(rho, eta, tmp);

    } while ((par->max_iter == 0 || iter < par->max_iter));

    free(a);
    free(b);
    free(c);
    free(alpha);

    free_spinor_field(rho);
    return iter;
}

int gcr(mshift_par *par, spinor_operator M, spinor_field *eta, spinor_field *phi) {
    int iter = 0;

    a = (hr_complex *)malloc((maxit_sapgcr + 1) * (maxit_sapgcr + 1) * sizeof(hr_complex));
    b = (hr_complex *)malloc((maxit_sapgcr + 1) * sizeof(hr_complex));
    c = (hr_complex *)malloc((maxit_sapgcr + 1) * sizeof(hr_complex));
    alpha = (hr_complex *)malloc((maxit_sapgcr + 1) * sizeof(hr_complex));

    rho = alloc_spinor_field(4 + 2 * maxit_sapgcr, eta->type);
    tmp = rho + 1;
    chi = tmp + 1;
    xi = chi + maxit_sapgcr + 1;

    double const innorm2 = sqnorm(eta);

    M(tmp, phi);
    sub(rho, eta, tmp);

    do {
        int k;
        for (k = 0; k <= maxit_sapgcr; k++) {
            copy(&xi[k], rho); // no preconditioning
            M(&chi[k], &xi[k]);

            for (int l = 0; l < k; l++) {
                a[l + (maxit_sapgcr + 1) * k] = prod(&chi[l], &chi[k]);
                mulc_add_assign(&chi[k], -a[l + (maxit_sapgcr + 1) * k], &chi[l]);
            }

            b[k] = sqrt(sqnorm(&chi[k]));
            mul(&chi[k], 1. / b[k], &chi[k]);
            c[k] = prod(&chi[k], rho);
            mulc_add_assign(rho, -c[k], &chi[k]);

            if (sqnorm(rho) < par->err2 * innorm2) { break; }

            iter++;
        }

        if (k > maxit_sapgcr) { k--; }

        for (int l = k; l >= 0; l--) {
            hr_complex sum = 0;
            for (int i = l + 1; i <= k; i++) {
                sum += a[l + (maxit_sapgcr + 1) * i] * alpha[i];
            }
            alpha[l] = (c[l] - sum) / b[l];
        }

        for (int l = 0; l <= k; l++) {
            mulc_add_assign(phi, alpha[l], &xi[l]);
        }

        if (k < maxit_sapgcr) { break; }

        M(tmp, phi);
        sub(rho, eta, tmp);
    } while ((par->max_iter == 0 || iter < par->max_iter));

    free(a);
    free(b);
    free(c);
    free(alpha);

    free_spinor_field(rho);
    return iter;
}
