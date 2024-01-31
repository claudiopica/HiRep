/*******************************************************************************
*
* Check of GPU boundary conditions (antiperiodic)
* NOCOMPILE= !WITH_GPU
*
*******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int errors = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    suNf_field *tmp = alloc(tmp, 1, &glattice);

    copy_from_gpu_suNf_field(u_gauge_f);
    apply_BCs_on_represented_gauge_field();

    // Now CPU operation

#ifdef BC_T_ANTIPERIODIC
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
#endif
#ifdef BC_X_ANTIPERIODIC
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
#endif
#ifdef BC_Y_ANTIPERIODIC
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
#endif
#ifdef BC_Z_ANTIPERIODIC
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
#endif

    copy_cpu(tmp, u_gauge_f);
    copy_to_gpu_suNf_field(tmp);
    sub_assign(tmp, u_gauge_f);
    double diffnorm = max(tmp);
    errors = check_diff_norm(diffnorm, 1e-14);

    finalize_process();
    return errors;
}
