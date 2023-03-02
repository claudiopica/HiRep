#ifndef WORK_SPACE_H
#define WORK_SPACE_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

//work_space.c
/* Workspace database*/
int iup_wrk(int site, int dir);
int idn_wrk(int site, int dir);
suNg *pu_gauge_wrk(int site, int dir);
suNg_field *u_gauge_wrk();
void reset_wrk_pointers();
void set_wrk_space(int i);
void set_wrk_space_and_pointers(int i, suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out);
int reserve_wrk_space();
int reserve_wrk_space_with_pointers(suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out);
void release_wrk_space(int id_release);
void free_wrk_space();

#ifdef __cplusplus
}
#endif
#endif //WORK_SPACE_H
