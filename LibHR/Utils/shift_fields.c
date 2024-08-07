/*
Routine to geometrically shift the fields (both spinor and/or gauge) 
The routine is MPI friendly and allows for any shift size.
If any input is left null the routine acts only the other field.
*/

#include "utils.h"
#include "libhr_core.h"
#include "error.h"
#include "memory.h"
#include "inverters.h"
#include "io.h"

void shift_fields(int *shift, spinor_field *sin, suNg_field *uin, spinor_field *sout, suNg_field *uout) {
    for (int i = 0; i < 4; i++) {
        error(shift[i] < 0, 1, "shift_fields [shift_fields.c]",
              "Error the shift vector components must all be greater than zero.");
    }
#if defined(BC_T_OPEN) || defined(BC_T_SF_ROTATED) || defined(BC_T_SF)
    if (shift[0] != 0) {
        shift[0] = 0;
        lprintf("SHIFT_FIELDS", 0, "Ignoring the time shift due to incompatible BC\n");
    }
#endif

    int total_shift = shift[0] + shift[1] + shift[2] + shift[3];
    if (total_shift == 0) {
        if (sin != NULL) { copy_spinor_field_cpu(sout, sin); }

        if (uin != NULL) { copy_suNg_field_cpu(uout, uin); }
    }

    if (sin != NULL) {
        if (sin->type == &glattice) {
            _TWO_SPINORS_MATCHING(sin, sout);
        } else {
            if (total_shift % 2 == 0) {
                _TWO_SPINORS_MATCHING(sin, sout);
            } else {
                error((sin)->type == (sout)->type, 1, "shift_fields [shift_fields.c]",
                      "Odd shift implies that spinor type must not match!");
            }
        }
    }

    suNg_field *ubuf[2], *utmp[2];
    spinor_field *sbuf[2], *stmp[2];

    int dd, x0, x1, x2, x3, mu, ipin, ipout;

    if (uin != NULL) {
        ubuf[0] = alloc_suNg_field(&glattice);
        ubuf[1] = alloc_suNg_field(&glattice);
    } else {
        ubuf[0] = ubuf[1] = NULL;
    }

    if (sin != NULL) {
        sbuf[0] = alloc_spinor_field(1, &glattice);
        sbuf[1] = alloc_spinor_field(1, &glattice);
        zero_spinor_field_cpu(sbuf[0]);
        zero_spinor_field_cpu(sbuf[1]);
    } else {
        sbuf[0] = sbuf[1] = NULL;
    }

    utmp[0] = uin;
    utmp[1] = ubuf[1];

    stmp[0] = sin;
    stmp[1] = sbuf[1];

    if (uin != NULL) {
        start_sendrecv_suNg_field(uin);
        complete_sendrecv_suNg_field(uin);
    }
    if (sin != NULL) {
        start_sendrecv_spinor_field(sin);
        complete_sendrecv_spinor_field(sin);
    }

    for (int i = 0; i < total_shift; i++) {
        if (i == total_shift - 1 && total_shift != 1) {
            utmp[1] = uout;
            stmp[1] = sout;
        }

        if (i < shift[0]) {
            dd = 0;
        } else if (i < shift[0] + shift[1]) {
            dd = 1;
        } else if (i < shift[0] + shift[1] + shift[2]) {
            dd = 2;
        } else {
            dd = 3;
        }

        for (x0 = 0; x0 < T_EXT; x0++) {
            for (x1 = 0; x1 < X_EXT; x1++) {
                for (x2 = 0; x2 < Y_EXT; x2++) {
                    for (x3 = 0; x3 < Z_EXT; x3++) {
                        ipin = ipt_ext(x0, x1, x2, x3);
                        ipout = iup(ipin, dd);
                        if (ipin != -1 && ipout != -1) {
                            if (uin != NULL) {
                                for (mu = 0; mu < 4; mu++) {
                                    *((utmp[1]->ptr) + coord_to_index(ipout, mu)) =
                                        *((utmp[0]->ptr) + coord_to_index(ipin, mu));
                                }
                            }
                            if (sin != NULL) {
                                if (stmp[0]->type->master_shift <= ipin &&
                                    stmp[0]->type->master_shift + stmp[0]->type->gsize_spinor > ipin &&
                                    stmp[1]->type->master_shift <= ipout &&
                                    stmp[1]->type->master_shift + stmp[1]->type->gsize_spinor > ipout) {
                                    *_FIELD_AT(stmp[1], ipout) = *_FIELD_AT(stmp[0], ipin);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (uin != NULL) {
            start_sendrecv_suNg_field(utmp[1]);
            complete_sendrecv_suNg_field(utmp[1]);
        }
        if (sin != NULL) {
            start_sendrecv_spinor_field(stmp[1]);
            complete_sendrecv_spinor_field(stmp[1]);
        }

        if (i % 2 == 0) {
            utmp[0] = ubuf[1];
            utmp[1] = ubuf[0];

            stmp[0] = sbuf[1];
            stmp[1] = sbuf[0];
        } else {
            utmp[0] = ubuf[0];
            utmp[1] = ubuf[1];

            stmp[0] = sbuf[0];
            stmp[1] = sbuf[1];
        }
    }
    if (uin != NULL) {
        if (total_shift == 1) { copy_suNg_field(uout, utmp[0]); }

        free_suNg_field(ubuf[0]);
        free_suNg_field(ubuf[1]);
    }
    if (sin != NULL) {
        if (total_shift == 1) { copy_spinor_field(sout, stmp[0]); }

        free_spinor_field(sbuf[0]);
        free_spinor_field(sbuf[1]);
    }
}