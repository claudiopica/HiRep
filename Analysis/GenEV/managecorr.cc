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
#include <iomanip> // std::setprecision
#include "type.h"
#include "utils.h"
#include "listrange.h"
#include "read_arg.h"
#include "matrix_algebra.h"
#include "fit.h"
using namespace std;

oplist::iterator iter_range1;
oplist::iterator iter_range2;

int main(int argc, char **argv)
{
  par apar;

  cout << "[INFO][MAIN] Reading input parameters ... " << endl;

  read_arg(&apar, argc, argv);

  print_parameters(&apar);

  oplist mylist(&apar);

  int cor_size = apar.ndt * apar.numop * apar.numop;
  double *avg_cor = new double[cor_size];
  for (int i = 0; i < cor_size; i++)
    avg_cor[i] = 0.0;

  cor_size *= apar.numbinjk;
  double *jk_cor = new double[cor_size];
  for (int i = 0; i < cor_size; i++)
    jk_cor[i] = 0.0;

  int vev_size = apar.numop * apar.numbinjk;
  double *jk_vev = NULL;

  if (apar.vflag == 1)
  {
    jk_vev = new double[vev_size];
    for (int i = 0; i < vev_size; i++)
      jk_vev[i] = 0.0;
  }

  int vev_counter = 0, cor_counter = 0;

  create_cor_file_name(&line, &cor_name, apar.rapr_name);   /* create the name of the correlator file based on the IRREP name */
  read_data(TWOPT, &cor_counter, jk_cor, &apar); /* reads the measurements on the correlator file and fills the correlator matrix */

  if (apar.vflag == 1)
    read_data(ONEPT, &vev_counter, jk_vev, &apar);
//Antonio

  create_simmetrize_corr_avg_jk(&apar, avg_cor, jk_cor, jk_vev); /* create normalize and symmetrize the jackknife bins */

  double *avg_cor_mat_t0 = new double[apar.numop * apar.nblock * apar.numop * apar.nblock];

  set_invert_matrix(&apar, avg_cor, avg_cor_mat_t0);

  if (apar.Gflag == 1)
  {
    cout << "[INFO][MAIN] Starting the cut procedure. " << endl
         << endl;
    while (true)
    {
      cout << "[INFO][MAIN] Op base under test: " << endl;
      mylist.report();

      double *mat = mylist.cut(avg_cor_mat_t0, 1);
      int matsize = mylist.size();

      double mineig = MinEigval(mat, matsize);

      if (fabs(mineig) > 1.e-14)
      {
        mylist.back().is_set = true;
        cout << "[INFO][MAIN] Cut accepted." << endl
             << endl;
      }
      else
        cout << "[INFO][MAIN] Cut discarded." << endl
             << endl;

      delete[] mat;
      if (!mylist.add_op())
        break;
    }
    cout << "[INFO][MAIN] End of the cut procedure. " << endl;
  }
  // else { proceed with the cut selected from the command line }
  // note: do we want to write out the maximal cut somewhere? remember that this cut usually depends on the number of bins
  delete[] avg_cor_mat_t0;

  double *eigenval = new double[mylist.size()];
  double *eigenval_avg = new double[mylist.size()];
  double *eigenval_avg_err = new double[mylist.size()];

  double *jk_inverse = new double[mylist.size() * mylist.size()];
  double *jk_cm1c_dt = new double[mylist.size() * mylist.size()];

  int n_evec = min(mylist.size(), apar.n_states);

  double *eigenvec = new double[mylist.size() * n_evec];
  double *eigenvec_avg = new double[mylist.size() * n_evec];

  for (int i = 0; i < mylist.size(); i++)
    eigenval_avg[i] = eigenval_avg_err[i] = 0.0;

  for (int i = 0; i < mylist.size() * n_evec; i++)
    eigenvec_avg[i] = 0.0;

  double *diag_corr = new double[apar.ntmax * apar.numbinjk * n_evec];

  for (int jkstep = 0; jkstep < apar.numbinjk; jkstep++)
  {
    get_inverse_matrix(jkstep, &apar, &mylist, jk_cor, jk_inverse);

    get_cm1c_matrix(jkstep, &apar, &mylist, jk_cor, jk_inverse, apar.p_diag, jk_cm1c_dt);

    get_eigenvalues(jk_cm1c_dt, mylist.size(), eigenval);

    get_eigenvectors(jk_cm1c_dt, mylist.size(), n_evec, eigenval, eigenvec);

    for (int i = 0; i < mylist.size(); i++)
    {
      eigenval_avg[i] += eigenval[i] / apar.numbinjk;
      eigenval_avg_err[i] += eigenval[i] * eigenval[i] / apar.numbinjk;
      eigenvec_avg[i] += eigenvec[i] / apar.numbinjk;
    }

    for (int nt = 0; nt < apar.ntmax; nt++)
    {
      get_cm1c_matrix(jkstep, &apar, &mylist, jk_cor, jk_inverse, nt, jk_cm1c_dt);

      for (int nev = 0; nev < n_evec; nev++)
      {

        double sum = 0.;
        for (int j = 0; j < mylist.size(); j++)
          for (int k = 0; k < mylist.size(); k++)
            sum += eigenvec[nev + j * n_evec] * jk_cm1c_dt[j + mylist.size() * k] * eigenvec[nev + k * n_evec];

        diag_corr[nt + apar.ntmax * jkstep + (apar.ntmax * apar.numbinjk) * (n_evec - nev - 1)] = sum;
      }
    }
  }

  for (int i = 0; i < mylist.size(); i++)
    eigenval_avg_err[i] = sqrt((apar.numbinjk - 1.0) / apar.numbinjk) * sqrt(eigenval_avg_err[i] - eigenval_avg[i] * eigenval_avg[i]);
  cout << endl
       << "[INFO][MAIN] Eigenvalues:" << endl;
  for (int i = 0; i < mylist.size(); i++)
    cout << "\t" << i << " " << eigenval_avg[i] << " " << eigenval_avg_err[i] << endl;
  cout << endl
       << endl;

  cout << "[INFO][MAIN] Evaluating the masses corresponding to the last " << n_evec << " eigenvalues." << endl
       << endl;

  for (int nev = 0; nev < n_evec; nev++)
  {
    cout << "#---------------------------------------------------------------" << endl;
    cout << "[INFO][MAIN] Eigenvalue:"
         << "\t"
         << " " << eigenval_avg[n_evec - nev - 1] << " +/- " << eigenval_avg_err[n_evec - nev - 1] << endl
         << endl;

    fitexp(&apar, diag_corr + (apar.ntmax * apar.numbinjk) * nev);
  }

  delete[] jk_cor;
  delete[] avg_cor;
  delete[] jk_cm1c_dt;
  delete[] diag_corr;

  return 0;
}
