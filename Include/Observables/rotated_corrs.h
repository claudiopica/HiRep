// Header file for:
// - rotated_corrs_gp.c
// - rotated_corrs_lp.c
// - rotated_corrs_gm.c
// - rotated_corrs_lm.c

#ifndef ROTATED_CORRS_H
#define ROTATED_CORRS_H

#ifdef ROTATED_SF

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

//rotated_corrs_gp.c
typedef struct
{
  hr_complex ***g1_ij, ***g2_ij, ***g3_ij, ***g4_ij, *g1, *g2, *g3, *g4, **M;

  hr_complex *l11, *l12, *l13;
  hr_complex *l21, *l22, *l23;
  hr_complex *l31, *l32, *l33;
  hr_complex *l41, *l42, *l43;

  hr_complex ***l11_ij, ***l12_ij, ***l13_ij;
  hr_complex ***l21_ij, ***l22_ij, ***l23_ij;
  hr_complex ***l31_ij, ***l32_ij, ***l33_ij;
} chisf_mem;

chisf_mem *init_rotated_corr_mem(void);
void rotated_gXuup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gXddp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gXudp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gXdup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtuup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtddp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtdup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtudp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1uup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1ddp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1udp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1dup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);

//rotated_corrs_lp.c
void rotated_lXuup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_lXddp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_lXdup(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_lXudp(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);

//rotated_corrs_gm.c
void rotated_gXuum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gXddm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gXudm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gXdum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtuum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtddm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtdum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_gvtudm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1uum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1ddm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1udm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_g1dum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);

//rotated_corrs_lm.c
void rotated_lXuum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_lXddm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_lXudm(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);
void rotated_lXdum(chisf_mem *corr_mem, suNf_spinor *chi, spinor_field *prop_uu, spinor_field *prop_dd);



#ifdef __cplusplus
	}
#endif
#endif //defined(ROTATED_SF)
#endif //ROTATED_CORRS_H
