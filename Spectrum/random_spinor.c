/*******************************************************************************
 *
 * Generate random cnfg
 *
 *******************************************************************************/

#include "libhr.h"
#include <string.h>

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

char input_filename[255] = "input_file";
char output_filename[255] = "random_spinor.out";

void read_cmdline(int argc, char *argv[]) {
    int i, ai = 0, ao = 0;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-i") == 0) {
            ai = i + 1;
        } else if (strcmp(argv[i], "-o") == 0) {
            ao = i + 1;
        }
    }

    if (ao != 0) { strcpy(output_filename, argv[ao]); }
    if (ai != 0) { strcpy(input_filename, argv[ai]); }
}

int main(int argc, char *argv[]) {
    char cnfg_filename[256];
    char tmp[256];
    spinor_field *sp;

    /* setup process id and communications */
    read_cmdline(argc, argv);
    setup_process(&argc, &argv);

    read_input(rlx_var.read, input_filename);

    /* logger setup */
    /* disable logger for MPI processes != 0 */
    if (PID != 0) { logger_disable(); }
    logger_setlevel(0, 40);
    sprintf(tmp, ">%s", output_filename);
    logger_stdout(tmp);
    sprintf(tmp, "err_%d", PID);
    freopen(tmp, "w", stderr);

    lprintf("MAIN", 0, "PId =  %d [world_size: %d]\n\n", PID, WORLD_SIZE);

    /* read & broadcast parameters */

    lprintf("MAIN", 0, "RLXD [%d,%d]\n", rlx_var.rlxd_level, rlx_var.rlxd_seed);
    rlxd_init(rlx_var.rlxd_level, rlx_var.rlxd_seed);

    lprintf("MAIN", 0, "Gauge group: SU(%d)\n", NG);
    lprintf("MAIN", 0, "Fermion representation: " REPR_NAME " [dim=%d]\n", NF);

    /* setup communication geometry */
    if (geometry_init() == 1) {
        finalize_process();
        return 0;
    }

    /* setup lattice geometry */
    geometry_mpi_eo();
    /* test_geometry_mpi_eo(); */

    init_BCs(NULL);

    /* alloc fields */
    sp = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(sp);
    apply_BCs_on_spinor_field(sp);

    sprintf(cnfg_filename, "sp_%dx%dx%dx%dNc%dNf%d", GLB_T, GLB_X, GLB_Y, GLB_Z, NG, NF);
    /*  write_spinor_field_eo_lexi(cnfg_filename,sp);*/
    write_spinor_field(cnfg_filename, sp);

    free_BCs();

    free_spinor_field_f(sp);

    return 0;
}
