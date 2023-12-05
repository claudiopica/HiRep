#include "inverters.h"
#include "libhr_core.h"
#include "memory.h"

void SAP_prec(int nu, inverter_ptr inv, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {
    spinor_field *res;
    spinor_field *stmp1;
    spinor_field *stmp2;

    res = alloc_spinor_field(3, in->type);
    stmp1 = res + 1;
    stmp2 = res + 2;

    while (nu--) {
        // Compute residue: res = in - M * out
        M(stmp1, out);
        sub_spinor_field(res, in, stmp1);
        zero_spinor_field(stmp1);
        empty_buffers(stmp1);

        // Invert black: stmp1 = M^-1 * res
        if (out->type == &glattice) {
            res->type = &glat_black;
            stmp1->type = &glat_black;
            inv(par, M, res, stmp1);
            res->type = &glattice;
            stmp1->type = &glattice;
        } else if (out->type == &glat_even) {
            res->type = &glat_even_black;
            stmp1->type = &glat_even_black;
            inv(par, M, res, stmp1);
            res->type = &glat_even;
            stmp1->type = &glat_even;
        } else {
            res->type = &glat_odd_black;
            stmp1->type = &glat_odd_black;
            inv(par, M, res, stmp1);
            res->type = &glat_odd;
            stmp1->type = &glat_odd;
        }

        // Update residue: res = res - M * stmp1
        M(stmp2, stmp1);
        sub_assign_spinor_field(res, stmp2);

        // Update solution: out = out + stmp1
        add_assign_spinor_field(out, stmp1);
        zero_spinor_field(stmp1);
        empty_buffers(stmp1);

        // Invert red: stmp1 = M^-1 * res
        if (out->type == &glattice) {
            res->type = &glat_red;
            stmp1->type = &glat_red;
            inv(par, M, res, stmp1);
            res->type = &glattice;
            stmp1->type = &glattice;
        } else if (out->type == &glat_even) {
            res->type = &glat_even_red;
            stmp1->type = &glat_even_red;
            inv(par, M, res, stmp1);
            res->type = &glat_even;
            stmp1->type = &glat_even;
        } else {
            res->type = &glat_odd_red;
            stmp1->type = &glat_odd_red;
            inv(par, M, res, stmp1);
            res->type = &glat_odd;
            stmp1->type = &glat_odd;
        }

        // Update solution: out = out + stmp1
        add_assign_spinor_field(out, stmp1);
    }

    // Free temporary spinors
    free_spinor_field(res);
}
