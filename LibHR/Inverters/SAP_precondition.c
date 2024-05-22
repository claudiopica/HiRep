#include "inverters.h"
#include "libhr_core.h"
#include "memory.h"
#include "utils.h"
#include "update.h"

#define _SWITCH_TO_BLACK(_type) \
    do {                        \
        (_type)->SAP = 1;       \
        (_type)->parity = 0;    \
    } while (0)

#define _SWITCH_TO_BLACK_BOUNDARY(_type) \
    do {                                 \
        (_type)->SAP = 2;                \
        (_type)->parity = 0;             \
    } while (0)

#define _SWITCH_TO_RED(_type) \
    do {                      \
        (_type)->SAP = 1;     \
        (_type)->parity = 1;  \
    } while (0)

#define _SWITCH_TO_RED_BOUNDARY(_type) \
    do {                               \
        (_type)->SAP = 2;              \
        (_type)->parity = 1;           \
    } while (0)

#define _SWITCH_TO_GLOBAL(_type) \
    do {                         \
        (_type)->SAP = 0;        \
    } while (0)

int MINRES_SAP(int nmr, spinor_operator_flt M, spinor_field_flt *eta, spinor_field_flt *psi, spinor_field_flt *r) {
    if (!eta->type->SAP || eta->type->parity == PARITY) {
        spinor_field_flt *Mpsi, *p;
        Mpsi = alloc(Mpsi, 2, eta->type);
        p = Mpsi + 1;
        hr_complex alpha;

        zero(psi);
        //M(Mpsi, psi);
        //sub(r, eta, Mpsi);
        copy(r, eta);

        for (int it = 0; it < nmr; it++) {
            M(p, r);
            alpha = prod(p, r);
            alpha /= sqnorm(p);
            mulc_add_assign(psi, alpha, r);
            mulc_add_assign(r, -alpha, p);
        }
    }
    return nmr;
}

int SAP_prec(int nmr, int ncy, inverter_ptr inv, mshift_par *par, spinor_operator M, spinor_field *eta, spinor_field *psi) {
#if defined(DPHI_FLT) && defined(WITH_GPU) && defined(WITH_MPI)
    spinor_field_flt *rho, *Mp, *xi, *res, *eta_flt, *psi_flt;
    int cgiter = 0;
    hr_complex alpha;

    rho = alloc(rho, 6, eta->type);
    xi = rho + 1;
    Mp = xi + 1;
    res = Mp + 1;
    eta_flt = res + 1;
    psi_flt = eta_flt + 1;

    assign_sd2s_gpu(eta_flt, eta);

    zero(psi_flt);

    while (ncy--) {
        _SWITCH_TO_GLOBAL(eta->type);
        D_flt(Mp, psi_flt);
        sub(rho, eta_flt, Mp);
        zero(xi);

        _SWITCH_TO_BLACK(eta->type);
        for (int it = 0; it < nmr; it++) {
            D_flt(Mp, rho);
            if (eta->type->parity == PARITY) {
                alpha = prod(Mp, rho);
                alpha /= sqnorm(Mp);
                mulc_add_assign(xi, alpha, rho);
                mulc_add_assign(rho, -alpha, Mp);
            }
        }
        add_assign(psi_flt, xi);

        _SWITCH_TO_GLOBAL(eta->type);
        D_flt(Mp, psi_flt);
        sub(rho, eta_flt, Mp);
        zero(xi);

        _SWITCH_TO_RED(eta->type);
        for (int it = 0; it < nmr; it++) {
            D_flt(Mp, rho);
            if (eta->type->parity == PARITY) {
                alpha = prod(Mp, rho);
                alpha /= sqnorm(Mp);
                mulc_add_assign(xi, alpha, rho);
                mulc_add_assign(rho, -alpha, Mp);
            }
        }
        add_assign(psi_flt, xi);
    }

    _SWITCH_TO_GLOBAL(eta->type);

    assign_s2sd_gpu(psi, psi_flt);

    // Free temporary spinors
    free_field(rho);
#else
    error(1, 1, __func__, "SAP not implemented for these compilation variables.\n");
#endif
    return cgiter;
}
