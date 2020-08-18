#ifndef _UTILS_
#define _UTILS_
#include <fstream>
#include <string>
#include "type.h"
#include <complex>

void read_2pt_def(par *apar);

void read_1pt_map(par *apar);

void read_data(double *cor_dataout, std::complex<double> *vev_dataout, par *apar);

void print_parameters(const par *apar);

void create_simmetrize_corr_jk(par *apar, double *cor_jk, double *cor, std::complex<double> *vev);

void deactivate_op(par *apar, int opid);

void deactivate_pull_corr_jk(par *apar, double *cor_jk, double pull, int dt);

void normalize_corr_jk(par *apar, double *cor_jk, int dt);

void get_cm1c_matrix(const par *apar, double *jk_mat, int jkstep, double *mat_inverse, int point, double *jk_corr_mat_t);

void check_list_dir_and_data_file(par *apar);

void print_diag_cor_avg(par *apar, int nop, double *cor_jk);

double *generate_cut_matrix(par *apar, double *cor_jk);

#endif // _UTILS_
