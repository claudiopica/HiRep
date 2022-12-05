#ifndef GLUEBALLS_H
#define GLUEBALLS_H
#include "hr_complex.h"
#include "logger.h"
#include "suN.h"

int **direct_spatial_rotations();
int **inverse_spatial_rotations();
void request_spatial_paths_evaluation();
void eval_all_glueball_ops(int t, hr_complex *numerical_op);
void measure_1pt_glueballs(int nblockingstart, int nblockingend, double *smear_val, hr_complex *gb_storage);
void eval_all_torellon_ops(int t, hr_complex *numerical_op, hr_complex ** polyf);
void measure_1pt_torellons(double *smear_val, hr_complex *tor_storage, hr_complex **pf);
void report_gb_group_setup();
void report_tor_group_setup();

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

typedef struct
{
    suNg *p;
    int ix;
    hr_complex tr;
} wilson_lines;

wilson_lines *polyleg(int ix, int d);

void collect_1pt_glueball_functions(cor_list *lcor, int nblocking, hr_complex *gb_storage);
void collect_1pt_torellon_functions(cor_list *lcor, hr_complex *tor_storage, hr_complex ** polyf);
    

#define total_n_glue_op 15
#define total_n_tor_op 8
#define npoly_dist 1
#endif
