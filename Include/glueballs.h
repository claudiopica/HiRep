#ifndef GLUEBALLS_H
#define GLUEBALLS_H
#include "hr_complex.h"
#include "logger.h"

int **direct_spatial_rotations();
int **inverse_spatial_rotations();
void request_spatial_paths_evaluation();
void eval_all_glueball_ops(int t, double complex *numerical_op);
void measure_1pt_glueballs(int nblocking, double *smear_val, double complex *gb_storage);
void report_op_group_setup();

typedef struct
{
    int t1;
    int t2;
    int id;
    int n_pairs;
} cor_points;

typedef struct
{
    int n_entries;
    cor_points *list;
    int n_corrs;
} cor_list;

void evaluate_1pt_functions(cor_list *lcor, int nblocking, double complex *gb_storage);
    

#define total_n_glue_op 16
#endif
