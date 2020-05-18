#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <list>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>
#include <iomanip>
#include <complex>

#include "type.h"
#include "utils.h"
#include "read_arg.h"
#include "matrix_algebra.h"
#include "fit.h"

int main(int argc, char **argv)
{
  par apar;

  std::cout << "[INFO][MAIN] Reading input parameters ... " << std::endl;

  read_arg(&apar, argc, argv);

  print_parameters(&apar);

  int vev_size = apar.nt * apar.numop * apar.numbinjk;
  std::complex<double> *vev_data = NULL;

  vev_data = new std::complex<double>[vev_size];
  for (int i = 0; i < vev_size; i++)
    vev_data[i] = 0.0;

  int cor_size = apar.numop * apar.numop * apar.ndt * apar.numbinjk;
  double *cor_data = new double[cor_size];
  for (int i = 0; i < cor_size; i++)
    cor_data[i] = 0.0;

  //reads the data and organize them in bin
  read_data(cor_data, vev_data, &apar);

  double *cor_jk = new double[cor_size];
  for (int i = 0; i < cor_size; i++)
    cor_jk[i] = 0.0;
  cor_size = apar.numop * apar.numop * apar.ndt;
  double *cor_avg = new double[cor_size];
  for (int i = 0; i < cor_size; i++)
    cor_avg[i] = 0.0;

  // create and symmetrize the jackknife bins
  create_simmetrize_corr_jk(&apar, cor_jk, cor_data, vev_data);

  delete[] cor_data;

  //normalize the correlator to 1 at dt
  normalize_corr_jk(&apar, cor_jk, 2);

  print_diag_cor_avg(&apar, apar.numop, cor_jk);

  // deactivate operators that are too noisy (pull smaller than) at some given timeslice dt
  deactivate_pull_corr_jk(&apar, cor_jk, 2.0, 2);

  //generating the reduced op matrix
  double *cor_jk_cut = generate_cut_matrix(&apar, cor_jk);

  print_diag_cor_avg(&apar, apar.activeop.size(), cor_jk_cut);

  int n_activeop = apar.activeop.size();
  int n_evec = std::min(apar.n_states, n_activeop);

  double *eigenval = new double[n_activeop];
  double *eigenval_avg = new double[n_activeop];
  double *eigenval_avg_err = new double[n_activeop];
  double *jk_inverse = new double[n_activeop * n_activeop];
  double *jk_cm1c_dt = new double[n_activeop * n_activeop];
  double *eigenvec = new double[n_activeop * n_evec];
  double *eigenvec_avg = new double[n_activeop * n_evec];
  double *diag_corr = new double[apar.ndt * apar.numbinjk * n_evec];

  for (int i = 0; i < n_activeop; i++)
    eigenval_avg[i] = eigenval_avg_err[i] = 0.0;

  for (int i = 0; i < n_activeop * n_evec; i++)
    eigenvec_avg[i] = 0.0;

  int i_inv = apar.corrdef[apar.p_inv].index;
  int i_diag = apar.corrdef[apar.p_diag].index;

  for (int jkstep = 0; jkstep < apar.numbinjk; jkstep++)
  {

    double *mat_t0 = cor_jk_cut + n_activeop * n_activeop * (i_inv + jkstep * apar.ndt);

    inverse(mat_t0, jk_inverse, n_activeop);

    get_cm1c_matrix(&apar, cor_jk_cut, jkstep, jk_inverse, i_diag, jk_cm1c_dt);

    get_eigenvalues(jk_cm1c_dt, n_activeop, eigenval);

    get_eigenvectors(jk_cm1c_dt, n_activeop, n_evec, eigenval, eigenvec);

    for (int i = 0; i < n_activeop; i++)
    {
      eigenval_avg[i] += eigenval[i] / apar.numbinjk;
      eigenval_avg_err[i] += eigenval[i] * eigenval[i] / apar.numbinjk;
      eigenvec_avg[i] += eigenvec[i] / apar.numbinjk;
    }

    for (int nt = 0; nt < apar.ndt; nt++)
    {
      get_cm1c_matrix(&apar, cor_jk_cut, jkstep, jk_inverse, nt, jk_cm1c_dt);

      for (int nev = 0; nev < n_evec; nev++)
      {

        double sum = 0.;
        for (int j = 0; j < n_activeop; j++)
          for (int k = 0; k < n_activeop; k++)
            sum += eigenvec[nev + j * n_evec] * jk_cm1c_dt[j + n_activeop * k] * eigenvec[nev + k * n_evec];

        diag_corr[nt + apar.ndt * (jkstep + apar.numbinjk * nev)] = sum;
      }
    }
  }

  for (int i = 0; i < n_activeop; i++)
    eigenval_avg_err[i] = sqrt((apar.numbinjk - 1.0) / apar.numbinjk) * sqrt(eigenval_avg_err[i] - eigenval_avg[i] * eigenval_avg[i]);

  std::cout << std::endl
            << "[INFO][MAIN] Eigenvalues:" << std::endl;

  for (int i = 0; i < n_activeop; i++)
    std::cout << "\t" << i << " " << eigenval_avg[i] << " " << eigenval_avg_err[i] << std::endl;

  std::cout << std::endl
            << std::endl;

  std::cout << "[INFO][MAIN] Evaluating the masses corresponding to the last " << n_evec << " eigenvalues.\n"
            << std::endl;

  for (int nev = 0; nev < n_evec; nev++)
  {
    std::cout << "#---------------------------------------------------------------" << std::endl;
    std::cout << "[INFO][MAIN] Eigenvalue:"
              << "\t"
              << " " << eigenval_avg[nev] << " +/- " << eigenval_avg_err[nev] << std::endl
              << std::endl;

    fitexp(&apar, diag_corr + (apar.ndt * apar.numbinjk) * nev);
  }

  return 0;
}
