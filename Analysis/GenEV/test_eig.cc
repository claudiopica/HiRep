/**
   Test the eigenvalue routines

 **/

#include <random>
#include <iostream>

#include <complex>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "matrix_algebra.h"
#include "utils.h"

using namespace std;

void eigen_from_gsl(double *m, int dim, gsl_vector_complex *eval, gsl_matrix_complex *evec);

int main()
{
  int dim = 4;
  int nev = 4;

  double *rm = new double[dim * dim];

  double *e_val = new double[dim];
  double *e_vec = new double[dim * nev];

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(1, 6);

  cout << "Testing eigenvalue routines random matrix order " << dim << endl;
  for (int i = 0; i < dim; ++i)
    for (int j = i; j < dim; ++j)
    {
      double dice_roll = distribution(generator); // generates number in the range 1..6
      rm[i + dim * j] = dice_roll;
      rm[j + dim * i] = dice_roll;

      if (i == j)
      {
        rm[j + dim * i] += 20;
      }
    }

  /** compute the eigenvalues  **/

  get_eigenvalues(rm, dim, e_val);

  get_eigenvectors(rm, dim, nev, e_val, e_vec);

  /** compare with results from GSL routiness **/
  gsl_vector_complex *eval = gsl_vector_complex_alloc(dim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(dim, dim);

  eigen_from_gsl(rm, dim, eval, evec);

  for (int i = 0; i < nev; ++i)
  {
    gsl_complex z = gsl_vector_complex_get(eval, dim - 1 - i);

    cout << "EIGENVALUE[" << i << "]= " << e_val[i] << "\t\t"
         << "GSL_EIG[" << i << "]= " << GSL_REAL(z) << " " << GSL_IMAG(z) << endl;

    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, dim - 1 - i);

    printf("Our routine's eigenvector\t\tGSL eigenvector \n");
    for (int j = 0; j < dim; ++j)
    {
      gsl_complex z =
          gsl_vector_complex_get(&evec_i.vector, j);

      printf("%g\t\t %g + %gi\n", e_vec[i + j * nev], GSL_REAL(z), GSL_IMAG(z));
    }
  }

  /** free up memory **/
  delete[] rm;
  delete[] e_vec;
  delete[] e_val;

  return 0;
}

/**
    Use the gun science library to compute the eigenvalues and eiegenvectors.

 **/

void eigen_from_gsl(double *m, int dim, gsl_vector_complex *eval, gsl_matrix_complex *evec)
{
  cout << "Testing eigenvalues with GSL \n";
  cout << "GSL == GNU science library  \n";

  gsl_matrix_view m_gsl = gsl_matrix_view_array(m, dim, dim);

  gsl_eigen_nonsymmv_workspace *gwk = gsl_eigen_nonsymmv_alloc(dim);

  int check = gsl_eigen_nonsymmv(&m_gsl.matrix, eval, evec, gwk);

  if (check != 0)
  {
    cout << "Error with call to gsl_eigen_nonsymmv" << endl;
    exit(0);
  }

  check = gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  if (check != 0)
  {
    cout << "Error with call to gsl_eigen_nonsymmmv_sort" << endl;
    exit(0);
  }
  gsl_eigen_nonsymmv_free(gwk);
}
