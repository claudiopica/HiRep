#ifndef SOURCES_H
#define SOURCES_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

//sources.c
void create_point_source(spinor_field *source, int tau, int color);
void create_full_point_source(spinor_field *source, int tau);
void create_point_source_loc(spinor_field *source, int t, int x, int y, int z, int color);
int create_diluted_source_equal_eo(spinor_field *source);
void create_diluted_source_equal_atau_eo(spinor_field *source, int tau);
int create_diluted_source_equal(spinor_field *source);
void create_diluted_source_equal_spinorfield1(spinor_field *source, int tau);
void create_diluted_source_equal_atau(spinor_field *source, int tau);
void create_noise_source_equal_eo(spinor_field *source);
void create_noise_source_equal_oe(spinor_field *source);
void create_diluted_source_equal_atau_col(spinor_field *source, int tau, int col);
void create_noise_source_equal_col_dil(spinor_field *source, int col);
void create_gauge_fixed_wall_source(spinor_field *source, int tau, int color);
void create_gauge_fixed_momentum_source(spinor_field *source, int pt, int px, int py, int pz, int color);
void create_sequential_source(spinor_field *source, int tf, spinor_field *prop);
void create_sequential_source_stoch(spinor_field *source, int tf, spinor_field *prop);
void restrict_timeslice(spinor_field *source, int tf, spinor_field *prop);
void create_diluted_volume_source(spinor_field *source, int parity_component, int mod);
void create_z2_volume_source(spinor_field *source);
void add_momentum(spinor_field *out, spinor_field *in, int px, int py, int pz);
void zero_even_or_odd_site_spinorfield(spinor_field *source, int nspinor, int eo);

#ifdef __cplusplus
	}
#endif
#endif //SOURCES_H
