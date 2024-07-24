/***************************************************************************\
 * Copyright (c) 2023, Sofie Martins                                         *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/******************************************************************************
*
* NOCOMPILE= WITH_GPU
* NOCOMPILE= WITH_MPI
* NOCOMPILE= !WITH_NEW_GEOMETRY
*
******************************************************************************/

/**
 * @file stitch.c
 * @brief Stitch configurations from smaller lattices together to larger configs.
 */

#include "libhr.h"
#include <string.h>

// TODO: add this to suN.h
#define _suNg_star(u, v)                 \
    _complex_star((u).c[0], (v).c[0]);   \
    _complex_star((u).c[1], (v).c[1]);   \
    _complex_star((u).c[2], (v).c[2]);   \
    _complex_star((u).c[3], (v).c[3]);   \
    _complex_star((u).c[4], (v).c[4]);   \
    _complex_star((u).c[5], (v).c[5]);   \
    _complex_star((u).c[6], (v).c[6]);   \
    _complex_star((u).c[7], (v).c[7]);   \
    _complex_star((u).c[8], (v).c[8]);   \
    _complex_star((u).c[9], (v).c[9]);   \
    _complex_star((u).c[10], (v).c[10]); \
    _complex_star((u).c[11], (v).c[11]); \
    _complex_star((u).c[12], (v).c[12]); \
    _complex_star((u).c[13], (v).c[13]); \
    _complex_star((u).c[14], (v).c[14]); \
    _complex_star((u).c[15], (v).c[15])

int GLB_T_OUT, GLB_X_OUT, GLB_Y_OUT, GLB_Z_OUT, CP_T;

#define init_input_config(varname)                                         \
    {                                                                      \
        .read = {                                                          \
            { "file_in", "file_in = %s", STRING_T, &((varname).file[0]) }, \
        }                                                                  \
    }

#define init_output_config(varname)                                         \
    {                                                                       \
        .read = {                                                           \
            { "GLB_T_OUT", "GLB_T_OUT = %d", INT_T, &GLB_T_OUT },           \
            { "GLB_X_OUT", "GLB_X_OUT = %d", INT_T, &GLB_X_OUT },           \
            { "GLB_Y_OUT", "GLB_Y_OUT = %d", INT_T, &GLB_Y_OUT },           \
            { "GLB_Z_OUT", "GLB_Z_OUT = %d", INT_T, &GLB_Z_OUT },           \
            { "CP_transform", "CP_transform = %d", INT_T, &CP_T },          \
            { "file_in", "file_out = %s", STRING_T, &((varname).file[0]) }, \
        }                                                                   \
    }

typedef struct {
    char file[256];
    input_record_t read[2];
} input_file_var;

typedef struct {
    char file[256];
    input_record_t read[7];
} output_file_var;

// save ipt from smaller input lattice
static int *ipt_old;

// We need more general lexi and ipt functions
// for general boxes [T,X,Y,Z]
#define _lexi_convert(_T, _X, _Y, _Z, _t, _x, _y, _z) \
    (((((_t) % (_T)) * (_X) + ((_x) % (_X))) * (_Y) + ((_y) % (_Y))) * (_Z) + ((_z) % (_Z)))
#define ipt_box(_T, _X, _Y, _Z, _t, _x, _y, _z) ipt_old[_lexi_convert(_T, _X, _Y, _Z, _t, _x, _y, _z)];

static int p_transformed_index(int t, int x, int y, int z, int T, int X, int Y, int Z, int p1, int p2, int p3, int mu) {
    int x_loc = x;
    int y_loc = y;
    int z_loc = z;
    if (mu == 0) {
        if (p1) { x_loc = X - 1 - (x % X); }
        if (p2) { y_loc = Y - 1 - (y % Y); }
        if (p3) { z_loc = Z - 1 - (z % Z); }
    } else {
        if (p1) {
            x_loc = X - 2 - (x % X);
            if (x_loc < 0) { x_loc = X - 1; }
        }
        if (p2) {
            y_loc = Y - 2 - (y % Y);
            if (y_loc < 0) { y_loc = Y - 1; }
        }
        if (p3) {
            z_loc = Z - 2 - (z % Z);
            if (z_loc) { z_loc = Z - 1; }
        }
    }
    return ipt_box(T, X, Y, Z, t, x_loc, y_loc, z_loc);
}

int main(int argc, char *argv[]) {
    setup_process(&argc, &argv);
    input_file_var input_var = init_input_config(input_var);
    output_file_var output_var = init_output_config(output_var);

    char *inputfile = get_input_filename();

    read_input(input_var.read, inputfile);
    read_input(output_var.read, inputfile);

    setup_gauge_fields();

    lprintf("INFO", 0, "Input Config dimensions given: (T, X, Y, Z) = (%d, %d, %d, %d)\n", GLB_T, GLB_X, GLB_Y, GLB_Z);
    lprintf("INFO", 0, "Input filename: %s\n", input_var.file);
    lprintf("INFO", 0, "Output Config dimensions given: (T, X, Y, Z) = (%d, %d, %d, %d)\n", GLB_T_OUT, GLB_X_OUT, GLB_Y_OUT,
            GLB_Z_OUT);
    lprintf("INFO", 0, "Output filename: %s\n", output_var.file);

    read_gauge_field(input_var.file);

    lprintf("SANITY", 0, "gfield sqnorm: %0.6e\n", sqnorm(u_gauge));

    int T_old = GLB_T;
    int X_old = GLB_X;
    int Y_old = GLB_Y;
    int Z_old = GLB_Z;
    size_t old_field_size = 4 * glattice.gsize_gauge;

    ipt_old = (int *)malloc(GLB_VOLUME * sizeof(int));
    memcpy(ipt_old, ipt, GLB_VOLUME * sizeof(int));

    suNg *tmp = malloc(old_field_size * sizeof(suNg));

    memcpy(tmp, u_gauge->ptr, 4 * glattice.gsize_gauge * sizeof(suNg));

    GLB_T = GLB_T_OUT;
    GLB_X = GLB_X_OUT;
    GLB_Y = GLB_Y_OUT;
    GLB_Z = GLB_Z_OUT;

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    free_gfield_f(u_gauge_f);
#endif
    GLB_VOL3 = ((long int)GLB_X) * ((long int)GLB_Y) * ((long int)GLB_Z);
    GLB_VOLUME = GLB_VOL3 * ((long int)GLB_T);

    T = GLB_T;
    X = GLB_X;
    Y = GLB_Y;
    Z = GLB_Z;

    VOL3 = GLB_VOL3;
    VOLUME = GLB_VOLUME;

    X_EXT = X;
    Y_EXT = Y;
    Z_EXT = Z;
    T_EXT = T;

#ifdef WITH_NEW_GEOMETRY
    define_geometry();
#else
    geometry_mpi_eo();
#endif
    setup_gauge_fields();
    random_u(u_gauge);

    for (int t = 0; t < T; t++) {
        for (int x = 0; x < X; x++) {
            for (int y = 0; y < Y; y++) {
                for (int z = 0; z < Z; z++) {
                    int b0 = t / T_old;
                    int b1 = x / X_old;
                    int b2 = y / Y_old;
                    int b3 = z / Z_old;
                    int p1 = b1 % 2;
                    int p2 = b2 % 2;
                    int p3 = b3 % 2;
                    int sublat_parity = 0;
                    if (CP_T) { sublat_parity = (b0 + b1 + b2 + b3) % 2; }

                    for (int mu = 0; mu < 4; mu++) {
                        int idx = p_transformed_index(t, x, y, z, T_old, X_old, Y_old, Z_old, p1, p2, p3, mu);
                        int idx_new = ipt(t, x, y, z);
                        suNg in = *_4FIELD_AT_PTR(tmp, idx, mu, 0);
                        suNg out;
                        out = in;
                        if (mu == 0) { out = in; }
                        if ((mu == 1) && p1) { _suNg_minus(out, in); }
                        if (mu == 2 && p2) { _suNg_minus(out, in); }
                        if (mu == 3 && p3) { _suNg_minus(out, in); }
                        *_4FIELD_AT(u_gauge, idx_new, mu) = out;
                    }
                }
            }
        }
    }

    double plaq = avr_plaquette();
    lprintf("PLAQ", 0, "Check output plaquette: %0.15e\n", plaq);
    write_gauge_field(output_var.file);

    finalize_process();
    return 0;
}