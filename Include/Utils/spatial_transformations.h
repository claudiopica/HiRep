#ifndef SPATIAL_TRANSFORMATIONS_H
#define SPATIAL_TRANSFORMATIONS_H
#ifdef __cplusplus
extern "C" {
#endif

/* Spatial Trasformations*/

//global variables
extern int *active_slices_list;
extern int *glbT_to_active_slices;
extern int n_active_slices;

/* Spatial blocking */
typedef enum { NEW_SBLK = 1, CONT_SBLK = 0 } eval_spat_block;

void initialize_spatial_active_slices(int *tlist);
void free_spatial_active_slices();
int spatial_blocking_wrkspace(eval_spat_block eval, unsigned int level);
int single_level_spatial_blocking_wrkspace(int wrk_in);
/* Spatial rotation*/
void assign_spatial_rotated_wrkspace(int *map, int idx_wrkspace);
/* Spatial APE smearing*/
int spatial_APE_smear_wrkspace(double *smear_val, int wrkspace_in);

#ifdef __cplusplus
}
#endif
#endif //SPATIAL_TRANSFORMATIONS_H
