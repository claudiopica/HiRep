#ifndef MEASURE_RENORMALIZATION_H
#define MEASURE_RENORMALIZATION_H

#include "spinor_field.h"
#include "Utils/data_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Same storage_switch/data_storage_array convention as measure_spectrum_pt
 * (see Include/Observables/spectrum.h, LibHR/Observables/meson_measurements.c).
 * When swc==STORE, *ret is allocated (1 element) and filled with the
 * already global-summed, 1/GLB_VOLUME-normalised 3-point functions,
 * shaped [channel, mass_idx, row, col, 0-or-1] (last dim: 0=real, 1=imag).
 * Look up a channel's index by name with renorm_channel_index() (e.g.
 * "Sin", "Sout", "id", "g5" -- see measure_renormalization.c for the full
 * list) rather than a hardcoded number, then read a matrix entry as:
 *   int idx_re[5] = {channel, mass_idx, row, col, 0};
 *   int idx_im[5] = {channel, mass_idx, row, col, 1};
 *   hr_complex val = *data_storage_element(*ret, 0, idx_re)
 *                  + I * (*data_storage_element(*ret, 0, idx_im));
 * Caller owns *ret and must free_data_storage(*ret) when done.
 * When swc==DONTSTORE, *ret is left untouched (pass NULL) and behaviour
 * is identical to before this parameter existed. */
void measure_renormalization(spinor_field *psi_in, spinor_field *psi_out, int nm, int pt_in, int px_in, int py_in, int pz_in,
                             int pt_out, int px_out, int py_out, int pz_out, storage_switch swc, data_storage_array **ret);

/* Index of channel `name` (e.g. "Sin", "Sout", "id", "g5", "g0" ...) into
 * the array measure_renormalization()'s STORE output is indexed by, or -1
 * if `name` isn't a known channel. The channel list itself stays private
 * to measure_renormalization.c -- this is the only way to get an index
 * from outside that file. */
int renorm_channel_index(const char *name);

void print_renormalization(int conf, int nm, double *mass, char *label, int pt_in, int px_in, int py_in, int pz_in, int pt_out,
                           int px_out, int py_out, int pz_out);

#ifdef __cplusplus
}
#endif
#endif //MEASURE_RENORMALIZATION_H
