#ifndef GLUEBALLS_H
#define GLUEBALLS_H
#include "hr_complex.h"
#include "logger.h"

int **direct_spatial_rotations();
int **inverse_spatial_rotations();
void request_spatial_paths_evaluation();
/*void eval_time_momentum_glueball_paths(int t, int px, int py, int pz);*/
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

void evaluate_correlators(cor_list *lcor, int nblocking, double complex *gb_storage);
    

#define dim_p_m1_m1_m1_Ir_1 1
#define n_OP_oneTr_p_m1_m1_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_m1_m1_Ir_1_C_1 0
#define dim_p_m1_m1_m1_Ir_2 1
#define n_OP_oneTr_p_m1_m1_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_m1_m1_Ir_2_C_1 0
#define dim_p_m1_m1_m1_Ir_3 2
#define n_OP_oneTr_p_m1_m1_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_m1_m1_Ir_3_C_1 0
#define dim_p_m1_m1_0_Ir_1 1
#define n_OP_oneTr_p_m1_m1_0_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_m1_0_Ir_1_C_1 0
#define dim_p_m1_m1_0_Ir_2 1
#define n_OP_oneTr_p_m1_m1_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_m1_0_Ir_2_C_1 0
#define dim_p_m1_m1_0_Ir_3 1
#define n_OP_oneTr_p_m1_m1_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_m1_0_Ir_3_C_1 0
#define dim_p_m1_m1_0_Ir_4 1
#define n_OP_oneTr_p_m1_m1_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_m1_m1_0_Ir_4_C_1 0
#define dim_p_m1_m1_1_Ir_1 1
#define n_OP_oneTr_p_m1_m1_1_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_m1_1_Ir_1_C_1 0
#define dim_p_m1_m1_1_Ir_2 1
#define n_OP_oneTr_p_m1_m1_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_m1_1_Ir_2_C_1 0
#define dim_p_m1_m1_1_Ir_3 2
#define n_OP_oneTr_p_m1_m1_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_m1_1_Ir_3_C_1 0
#define dim_p_m1_0_m1_Ir_1 1
#define n_OP_oneTr_p_m1_0_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_0_m1_Ir_1_C_1 0
#define dim_p_m1_0_m1_Ir_2 1
#define n_OP_oneTr_p_m1_0_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_0_m1_Ir_2_C_1 0
#define dim_p_m1_0_m1_Ir_3 1
#define n_OP_oneTr_p_m1_0_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_0_m1_Ir_3_C_1 0
#define dim_p_m1_0_m1_Ir_4 1
#define n_OP_oneTr_p_m1_0_m1_Ir_4_C_m1 0
#define n_OP_oneTr_p_m1_0_m1_Ir_4_C_1 0
#define dim_p_m1_0_0_Ir_1 1
#define n_OP_oneTr_p_m1_0_0_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_0_0_Ir_1_C_1 0
#define dim_p_m1_0_0_Ir_2 1
#define n_OP_oneTr_p_m1_0_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_0_0_Ir_2_C_1 0
#define dim_p_m1_0_0_Ir_3 2
#define n_OP_oneTr_p_m1_0_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_0_0_Ir_3_C_1 0
#define dim_p_m1_0_0_Ir_4 1
#define n_OP_oneTr_p_m1_0_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_m1_0_0_Ir_4_C_1 0
#define dim_p_m1_0_0_Ir_5 1
#define n_OP_oneTr_p_m1_0_0_Ir_5_C_m1 0
#define n_OP_oneTr_p_m1_0_0_Ir_5_C_1 0
#define dim_p_m1_0_1_Ir_1 1
#define n_OP_oneTr_p_m1_0_1_Ir_1_C_m1 1
void OP_oneTr_p_m1_0_1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_m1_0_1_Ir_1_C_1 2
void OP_oneTr_p_m1_0_1_Ir_1_C_1(double complex * numop);
#define dim_p_m1_0_1_Ir_2 1
#define n_OP_oneTr_p_m1_0_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_0_1_Ir_2_C_1 0
#define dim_p_m1_0_1_Ir_3 1
#define n_OP_oneTr_p_m1_0_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_0_1_Ir_3_C_1 0
#define dim_p_m1_0_1_Ir_4 1
#define n_OP_oneTr_p_m1_0_1_Ir_4_C_m1 0
#define n_OP_oneTr_p_m1_0_1_Ir_4_C_1 0
#define dim_p_m1_1_m1_Ir_1 1
#define n_OP_oneTr_p_m1_1_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_1_m1_Ir_1_C_1 0
#define dim_p_m1_1_m1_Ir_2 1
#define n_OP_oneTr_p_m1_1_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_1_m1_Ir_2_C_1 0
#define dim_p_m1_1_m1_Ir_3 2
#define n_OP_oneTr_p_m1_1_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_1_m1_Ir_3_C_1 0
#define dim_p_m1_1_0_Ir_1 1
#define n_OP_oneTr_p_m1_1_0_Ir_1_C_m1 2
void OP_oneTr_p_m1_1_0_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_m1_1_0_Ir_1_C_1 2
void OP_oneTr_p_m1_1_0_Ir_1_C_1(double complex * numop);
#define dim_p_m1_1_0_Ir_2 1
#define n_OP_oneTr_p_m1_1_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_1_0_Ir_2_C_1 0
#define dim_p_m1_1_0_Ir_3 1
#define n_OP_oneTr_p_m1_1_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_1_0_Ir_3_C_1 0
#define dim_p_m1_1_0_Ir_4 1
#define n_OP_oneTr_p_m1_1_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_m1_1_0_Ir_4_C_1 0
#define dim_p_m1_1_1_Ir_1 1
#define n_OP_oneTr_p_m1_1_1_Ir_1_C_m1 0
#define n_OP_oneTr_p_m1_1_1_Ir_1_C_1 0
#define dim_p_m1_1_1_Ir_2 1
#define n_OP_oneTr_p_m1_1_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_m1_1_1_Ir_2_C_1 0
#define dim_p_m1_1_1_Ir_3 2
#define n_OP_oneTr_p_m1_1_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_m1_1_1_Ir_3_C_1 0
#define dim_p_0_m1_m1_Ir_1 1
#define n_OP_oneTr_p_0_m1_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_0_m1_m1_Ir_1_C_1 0
#define dim_p_0_m1_m1_Ir_2 1
#define n_OP_oneTr_p_0_m1_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_m1_m1_Ir_2_C_1 0
#define dim_p_0_m1_m1_Ir_3 1
#define n_OP_oneTr_p_0_m1_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_m1_m1_Ir_3_C_1 0
#define dim_p_0_m1_m1_Ir_4 1
#define n_OP_oneTr_p_0_m1_m1_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_m1_m1_Ir_4_C_1 0
#define dim_p_0_m1_0_Ir_1 1
#define n_OP_oneTr_p_0_m1_0_Ir_1_C_m1 0
#define n_OP_oneTr_p_0_m1_0_Ir_1_C_1 0
#define dim_p_0_m1_0_Ir_2 1
#define n_OP_oneTr_p_0_m1_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_m1_0_Ir_2_C_1 0
#define dim_p_0_m1_0_Ir_3 2
#define n_OP_oneTr_p_0_m1_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_m1_0_Ir_3_C_1 0
#define dim_p_0_m1_0_Ir_4 1
#define n_OP_oneTr_p_0_m1_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_m1_0_Ir_4_C_1 0
#define dim_p_0_m1_0_Ir_5 1
#define n_OP_oneTr_p_0_m1_0_Ir_5_C_m1 0
#define n_OP_oneTr_p_0_m1_0_Ir_5_C_1 0
#define dim_p_0_m1_1_Ir_1 1
#define n_OP_oneTr_p_0_m1_1_Ir_1_C_m1 2
void OP_oneTr_p_0_m1_1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_0_m1_1_Ir_1_C_1 2
void OP_oneTr_p_0_m1_1_Ir_1_C_1(double complex * numop);
#define dim_p_0_m1_1_Ir_2 1
#define n_OP_oneTr_p_0_m1_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_m1_1_Ir_2_C_1 0
#define dim_p_0_m1_1_Ir_3 1
#define n_OP_oneTr_p_0_m1_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_m1_1_Ir_3_C_1 0
#define dim_p_0_m1_1_Ir_4 1
#define n_OP_oneTr_p_0_m1_1_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_m1_1_Ir_4_C_1 0
#define dim_p_0_0_m1_Ir_1 1
#define n_OP_oneTr_p_0_0_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_0_0_m1_Ir_1_C_1 0
#define dim_p_0_0_m1_Ir_2 1
#define n_OP_oneTr_p_0_0_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_0_m1_Ir_2_C_1 0
#define dim_p_0_0_m1_Ir_3 2
#define n_OP_oneTr_p_0_0_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_0_m1_Ir_3_C_1 0
#define dim_p_0_0_m1_Ir_4 1
#define n_OP_oneTr_p_0_0_m1_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_0_m1_Ir_4_C_1 0
#define dim_p_0_0_m1_Ir_5 1
#define n_OP_oneTr_p_0_0_m1_Ir_5_C_m1 0
#define n_OP_oneTr_p_0_0_m1_Ir_5_C_1 0
#define dim_p_0_0_0_Ir_1 1
#define n_OP_oneTr_p_0_0_0_Ir_1_C_m1 1
void OP_oneTr_p_0_0_0_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_0_0_0_Ir_1_C_1 2
void OP_oneTr_p_0_0_0_Ir_1_C_1(double complex * numop);
#define dim_p_0_0_0_Ir_2 1
#define n_OP_oneTr_p_0_0_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_2_C_1 0
#define dim_p_0_0_0_Ir_3 2
#define n_OP_oneTr_p_0_0_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_3_C_1 0
#define dim_p_0_0_0_Ir_4 3
#define n_OP_oneTr_p_0_0_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_4_C_1 0
#define dim_p_0_0_0_Ir_5 3
#define n_OP_oneTr_p_0_0_0_Ir_5_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_5_C_1 0
#define dim_p_0_0_0_Ir_6 1
#define n_OP_oneTr_p_0_0_0_Ir_6_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_6_C_1 0
#define dim_p_0_0_0_Ir_7 1
#define n_OP_oneTr_p_0_0_0_Ir_7_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_7_C_1 0
#define dim_p_0_0_0_Ir_8 2
#define n_OP_oneTr_p_0_0_0_Ir_8_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_8_C_1 0
#define dim_p_0_0_0_Ir_9 3
#define n_OP_oneTr_p_0_0_0_Ir_9_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_9_C_1 0
#define dim_p_0_0_0_Ir_10 3
#define n_OP_oneTr_p_0_0_0_Ir_10_C_m1 0
#define n_OP_oneTr_p_0_0_0_Ir_10_C_1 0
#define dim_p_0_0_1_Ir_1 1
#define n_OP_oneTr_p_0_0_1_Ir_1_C_m1 1
void OP_oneTr_p_0_0_1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_0_0_1_Ir_1_C_1 2
void OP_oneTr_p_0_0_1_Ir_1_C_1(double complex * numop);
#define dim_p_0_0_1_Ir_2 1
#define n_OP_oneTr_p_0_0_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_0_1_Ir_2_C_1 0
#define dim_p_0_0_1_Ir_3 2
#define n_OP_oneTr_p_0_0_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_0_1_Ir_3_C_1 0
#define dim_p_0_0_1_Ir_4 1
#define n_OP_oneTr_p_0_0_1_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_0_1_Ir_4_C_1 0
#define dim_p_0_0_1_Ir_5 1
#define n_OP_oneTr_p_0_0_1_Ir_5_C_m1 0
#define n_OP_oneTr_p_0_0_1_Ir_5_C_1 0
#define dim_p_0_1_m1_Ir_1 1
#define n_OP_oneTr_p_0_1_m1_Ir_1_C_m1 2
void OP_oneTr_p_0_1_m1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_0_1_m1_Ir_1_C_1 2
void OP_oneTr_p_0_1_m1_Ir_1_C_1(double complex * numop);
#define dim_p_0_1_m1_Ir_2 1
#define n_OP_oneTr_p_0_1_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_1_m1_Ir_2_C_1 0
#define dim_p_0_1_m1_Ir_3 1
#define n_OP_oneTr_p_0_1_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_1_m1_Ir_3_C_1 0
#define dim_p_0_1_m1_Ir_4 1
#define n_OP_oneTr_p_0_1_m1_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_1_m1_Ir_4_C_1 0
#define dim_p_0_1_0_Ir_1 1
#define n_OP_oneTr_p_0_1_0_Ir_1_C_m1 2
void OP_oneTr_p_0_1_0_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_0_1_0_Ir_1_C_1 2
void OP_oneTr_p_0_1_0_Ir_1_C_1(double complex * numop);
#define dim_p_0_1_0_Ir_2 1
#define n_OP_oneTr_p_0_1_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_1_0_Ir_2_C_1 0
#define dim_p_0_1_0_Ir_3 2
#define n_OP_oneTr_p_0_1_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_1_0_Ir_3_C_1 0
#define dim_p_0_1_0_Ir_4 1
#define n_OP_oneTr_p_0_1_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_1_0_Ir_4_C_1 0
#define dim_p_0_1_0_Ir_5 1
#define n_OP_oneTr_p_0_1_0_Ir_5_C_m1 0
#define n_OP_oneTr_p_0_1_0_Ir_5_C_1 0
#define dim_p_0_1_1_Ir_1 1
#define n_OP_oneTr_p_0_1_1_Ir_1_C_m1 2
void OP_oneTr_p_0_1_1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_0_1_1_Ir_1_C_1 2
void OP_oneTr_p_0_1_1_Ir_1_C_1(double complex * numop);
#define dim_p_0_1_1_Ir_2 1
#define n_OP_oneTr_p_0_1_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_0_1_1_Ir_2_C_1 0
#define dim_p_0_1_1_Ir_3 1
#define n_OP_oneTr_p_0_1_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_0_1_1_Ir_3_C_1 0
#define dim_p_0_1_1_Ir_4 1
#define n_OP_oneTr_p_0_1_1_Ir_4_C_m1 0
#define n_OP_oneTr_p_0_1_1_Ir_4_C_1 0
#define dim_p_1_m1_m1_Ir_1 1
#define n_OP_oneTr_p_1_m1_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_1_m1_m1_Ir_1_C_1 0
#define dim_p_1_m1_m1_Ir_2 1
#define n_OP_oneTr_p_1_m1_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_m1_m1_Ir_2_C_1 0
#define dim_p_1_m1_m1_Ir_3 2
#define n_OP_oneTr_p_1_m1_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_m1_m1_Ir_3_C_1 0
#define dim_p_1_m1_0_Ir_1 1
#define n_OP_oneTr_p_1_m1_0_Ir_1_C_m1 2
void OP_oneTr_p_1_m1_0_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_1_m1_0_Ir_1_C_1 2
void OP_oneTr_p_1_m1_0_Ir_1_C_1(double complex * numop);
#define dim_p_1_m1_0_Ir_2 1
#define n_OP_oneTr_p_1_m1_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_m1_0_Ir_2_C_1 0
#define dim_p_1_m1_0_Ir_3 1
#define n_OP_oneTr_p_1_m1_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_m1_0_Ir_3_C_1 0
#define dim_p_1_m1_0_Ir_4 1
#define n_OP_oneTr_p_1_m1_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_1_m1_0_Ir_4_C_1 0
#define dim_p_1_m1_1_Ir_1 1
#define n_OP_oneTr_p_1_m1_1_Ir_1_C_m1 0
#define n_OP_oneTr_p_1_m1_1_Ir_1_C_1 0
#define dim_p_1_m1_1_Ir_2 1
#define n_OP_oneTr_p_1_m1_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_m1_1_Ir_2_C_1 0
#define dim_p_1_m1_1_Ir_3 2
#define n_OP_oneTr_p_1_m1_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_m1_1_Ir_3_C_1 0
#define dim_p_1_0_m1_Ir_1 1
#define n_OP_oneTr_p_1_0_m1_Ir_1_C_m1 1
void OP_oneTr_p_1_0_m1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_1_0_m1_Ir_1_C_1 2
void OP_oneTr_p_1_0_m1_Ir_1_C_1(double complex * numop);
#define dim_p_1_0_m1_Ir_2 1
#define n_OP_oneTr_p_1_0_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_0_m1_Ir_2_C_1 0
#define dim_p_1_0_m1_Ir_3 1
#define n_OP_oneTr_p_1_0_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_0_m1_Ir_3_C_1 0
#define dim_p_1_0_m1_Ir_4 1
#define n_OP_oneTr_p_1_0_m1_Ir_4_C_m1 0
#define n_OP_oneTr_p_1_0_m1_Ir_4_C_1 0
#define dim_p_1_0_0_Ir_1 1
#define n_OP_oneTr_p_1_0_0_Ir_1_C_m1 1
void OP_oneTr_p_1_0_0_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_1_0_0_Ir_1_C_1 2
void OP_oneTr_p_1_0_0_Ir_1_C_1(double complex * numop);
#define dim_p_1_0_0_Ir_2 1
#define n_OP_oneTr_p_1_0_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_0_0_Ir_2_C_1 0
#define dim_p_1_0_0_Ir_3 2
#define n_OP_oneTr_p_1_0_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_0_0_Ir_3_C_1 0
#define dim_p_1_0_0_Ir_4 1
#define n_OP_oneTr_p_1_0_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_1_0_0_Ir_4_C_1 0
#define dim_p_1_0_0_Ir_5 1
#define n_OP_oneTr_p_1_0_0_Ir_5_C_m1 0
#define n_OP_oneTr_p_1_0_0_Ir_5_C_1 0
#define dim_p_1_0_1_Ir_1 1
#define n_OP_oneTr_p_1_0_1_Ir_1_C_m1 1
void OP_oneTr_p_1_0_1_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_1_0_1_Ir_1_C_1 2
void OP_oneTr_p_1_0_1_Ir_1_C_1(double complex * numop);
#define dim_p_1_0_1_Ir_2 1
#define n_OP_oneTr_p_1_0_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_0_1_Ir_2_C_1 0
#define dim_p_1_0_1_Ir_3 1
#define n_OP_oneTr_p_1_0_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_0_1_Ir_3_C_1 0
#define dim_p_1_0_1_Ir_4 1
#define n_OP_oneTr_p_1_0_1_Ir_4_C_m1 0
#define n_OP_oneTr_p_1_0_1_Ir_4_C_1 0
#define dim_p_1_1_m1_Ir_1 1
#define n_OP_oneTr_p_1_1_m1_Ir_1_C_m1 0
#define n_OP_oneTr_p_1_1_m1_Ir_1_C_1 0
#define dim_p_1_1_m1_Ir_2 1
#define n_OP_oneTr_p_1_1_m1_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_1_m1_Ir_2_C_1 0
#define dim_p_1_1_m1_Ir_3 2
#define n_OP_oneTr_p_1_1_m1_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_1_m1_Ir_3_C_1 0
#define dim_p_1_1_0_Ir_1 1
#define n_OP_oneTr_p_1_1_0_Ir_1_C_m1 2
void OP_oneTr_p_1_1_0_Ir_1_C_m1(double complex * numop);
#define n_OP_oneTr_p_1_1_0_Ir_1_C_1 2
void OP_oneTr_p_1_1_0_Ir_1_C_1(double complex * numop);
#define dim_p_1_1_0_Ir_2 1
#define n_OP_oneTr_p_1_1_0_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_1_0_Ir_2_C_1 0
#define dim_p_1_1_0_Ir_3 1
#define n_OP_oneTr_p_1_1_0_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_1_0_Ir_3_C_1 0
#define dim_p_1_1_0_Ir_4 1
#define n_OP_oneTr_p_1_1_0_Ir_4_C_m1 0
#define n_OP_oneTr_p_1_1_0_Ir_4_C_1 0
#define dim_p_1_1_1_Ir_1 1
#define n_OP_oneTr_p_1_1_1_Ir_1_C_m1 0
#define n_OP_oneTr_p_1_1_1_Ir_1_C_1 0
#define dim_p_1_1_1_Ir_2 1
#define n_OP_oneTr_p_1_1_1_Ir_2_C_m1 0
#define n_OP_oneTr_p_1_1_1_Ir_2_C_1 0
#define dim_p_1_1_1_Ir_3 2
#define n_OP_oneTr_p_1_1_1_Ir_3_C_m1 0
#define n_OP_oneTr_p_1_1_1_Ir_3_C_1 0
#define total_n_glue_op 30
#endif
