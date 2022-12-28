/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef MOREIO_H
#define MOREIO_H

#include "input_par.h"
#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

void read_gauge_field_mpieo_BE(char filename[]);
void write_gauge_field_mpieo_BE(char filename[]);

void read_gauge_field_mpieo_LE(char filename[]);
void write_gauge_field_mpieo_LE(char filename[]);

void read_gauge_field_eolexi_BE(char filename[]);
void write_gauge_field_eolexi_BE(char filename[]);

void read_gauge_field_eolexi_LE(char filename[]);
void write_gauge_field_eolexi_LE(char filename[]);

void read_gauge_field_milc(char filename[]);
void read_gauge_field_milc_no3row(char filename[]);
void read_gauge_field_ascii(char filename[]);
void read_gauge_field_fortran(char filename[]);

void read_gauge_field_openQCD(char filename[]);
void read_gauge_field_openQCD_SF(char filename[]);
void write_gauge_field_openQCD(char filename[]);
void write_gauge_field_hirep_pbc_to_obc(char filename[]);
void write_gauge_field_hirep_pbc_to_sf(char filename[]);

/* void write_spinor_field_eo_lexi(char filename[],spinor_field *sp); */
void write_spinor_field(char filename[],spinor_field *sp);
void read_spinor_field(char filename[],spinor_field *sp);

#ifdef __cplusplus
    }
#endif
#endif
