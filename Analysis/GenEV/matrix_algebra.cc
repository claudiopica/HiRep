#include <iostream>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include "matrix_algebra.h"
#define TINY 1.0e-20
using namespace std;

void lubksb(double *a, int n, int *indx, double *b)
{
  int i, ii = -1, ip, j;
  double sum;

  for (i = 0; i < n; i++)
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii + 1)
      for (j = ii; j <= i; j++)
        sum -= a[i * n + j] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n - 1; i >= 0; i--)
  {
    sum = b[i];
    for (j = i + 1; j < n; j++)
      sum -= a[i * n + j] * b[j];
    b[i] = sum / a[i * n + i];
  }
}

void ludcmp(double *a, int n, int *indx, double *d)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  double *vv = new double[n];

  *d = 1.0;
  for (i = 0; i < n; i++)
  {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs(a[i * n + j])) > big)
        big = temp;
    if (big == 0.0)
      cerr << "[ERROR][LUDCMP] Singular matrix!!" << endl;
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a[i * n + j];
      for (k = 0; k < i; k++)
        sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)
    {
      sum = a[i * n + j];
      for (k = 0; k < j; k++)
        sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)
    {
      for (k = 0; k < n; k++)
      {
        dum = a[imax * n + k];
        a[imax * n + k] = a[j * n + k];
        a[j * n + k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j * n + j] == 0.0)
      a[j * n + j] = TINY;
    if (j != n)
    {
      dum = 1.0 / (a[j * n + j]);
      for (i = j + 1; i < n; i++)
        a[i * n + j] *= dum;
    }
  }
  delete[] vv;
}

double Det(double *a, int n)
{
  double *b = new double[n * n];
  for (int i = 0; i < n * n; i++)
    b[i] = a[i];

  double dd;
  int *idx = new int[n];

  ludcmp(b, n, idx, &dd);

  for (int i = 0; i < n; i++)
    dd *= b[i * (n + 1)];

  delete[] idx;
  delete[] b;

  return dd;
}

double MinEigval(double *a, int n)
{
  double *b = new double[n * n];
  for (int i = 0; i < n * n; i++)
    b[i] = a[i];

  double dd;
  int *idx = new int[n];

  ludcmp(b, n, idx, &dd);

  dd = 100;

  for (int i = 0; i < n; i++)
  {
    if (fabs(dd) > fabs(b[i * (n + 1)]))
      dd = b[i * (n + 1)];
  }
  delete[] idx;
  delete[] b;

  return dd;
}

void inverse(double *a, double *ainv, int n)
{
  double *b = new double[n * n];
  for (int i = 0; i < n * n; i++)
    b[i] = a[i];

  int i, j;

  double dd;
  int *idx = new int[n];

  double *c = new double[n];

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
      ainv[i * n + j] = 0.;
    ainv[i + i * n] = 1.;
  }

  ludcmp(b, n, idx, &dd);

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < n; i++)
      c[i] = ainv[i * n + j];
    lubksb(b, n, idx, c);
    for (i = 0; i < n; i++)
      ainv[i * n + j] = c[i];
  }

  delete[] c;
  delete[] idx;
  delete[] b;
}

// Ordering algorithms: needed to order the eigenvalues of the Correlator matrix

int *orders(const int n, const double *vec)
{
  int *order = new int[n];
  double maxvecnum = 1000000.0;
  double smaller;
  int idx;
  double *dat = new double[n];

  for (int i = 0; i < n; i++)
  {
    dat[i] = vec[i];
    order[i] = 0;
  }

  for (int i = 0; i < n; i++)
  {
    smaller = maxvecnum;
    idx = 0;
    for (int j = 0; j < n; j++)
    {
      if (smaller > dat[j])
      {
        smaller = dat[j];
        idx = j;
      }
    }
    order[i] = idx;
    dat[idx] = maxvecnum + 1;
  }
  delete[] dat;

  return order;
}

void sorting(const int n, const double *vec, double *ordered_vec)
{
  //    int * order = new int(n);
  double maxvecnum = 1000000.0;
  double smaller;
  int idx;
  double *dat = new double[n];

  for (int i = 0; i < n; i++)
  {
    dat[i] = vec[i];
  }

  for (int i = 0; i < n; i++)
  {
    smaller = maxvecnum;
    idx = 0;
    for (int j = 0; j < n; j++)
    {
      if (smaller > dat[j])
      {
        smaller = dat[j];
        idx = j;
      }
    }
    ordered_vec[i] = vec[idx];
    dat[idx] = maxvecnum + 1;
  }
  delete[] dat;
}

#define SWAP(g, h) \
  {                \
    float y = (g); \
    (g) = (h);     \
    (h) = y;       \
  }
void quicksort(int n, double *v)
{
  int i, j;
  double separator;

  if (n < 2)
    return;
  separator = v[n / 2];
  i = 0;
  j = n - 1;
  while (j > i)
  {
    while (v[i] > separator)
      i++;
    while (v[j] < separator)
      j--;
    if (i < j)
    {
      SWAP(v[i], v[j]);
      i++;
    }
  }
  quicksort(i, v);
  quicksort(n - i, v + i);
}

void elmhes(double *a, int n)
{
  int m, j, i;
  double y, x;

  for (m = 1; m < (n - 1); m++)
  {
    x = 0.0;
    i = m;
    for (j = m; j < n; j++)
    {
      if (fabs(a[j * n + (m - 1)]) > fabs(x))
      {
        x = a[j * n + (m - 1)];
        i = j;
      }
    }
    if (i != m)
    {
      for (j = (m - 1); j < n; j++)
        SWAP(a[i * n + j], a[m * n + j]);
      for (j = 0; j < n; j++)
        SWAP(a[j * n + i], a[j * n + m]);
    }
    if (x != 0.0)
    {
      for (i = (m + 1); i < n; i++)
      {
        if ((y = a[i * n + (m - 1)]) != 0.0)
        {
          y /= x;
          a[i * n + (m - 1)] = y;
          for (j = m; j < n; j++)
            a[i * n + j] -= y * a[m * n + j];
          for (j = 0; j < n; j++)
            a[j * n + m] += y * a[j * n + i];
        }
      }
    }
  }
}

static int imaxarg1, imaxarg2;
#define IMAX(a, b) (imaxarg1 = (a), imaxarg2 = (b), (imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void hqr(double *a, int n, double *wr, double *wi)
{
  int nn, m, l, k, j, its, i, mmin;
  double z, y, x, w, v, u, t, s, r = 0.0, q = 0.0, p = 0.0, anorm;

  anorm = 0.0;
  for (i = 0; i < n; i++)
    for (j = IMAX(i - 1, 0); j < n; j++)
      anorm += fabs(a[i * n + j]);
  nn = n - 1;
  t = 0.0;
  while (nn > -1)
  {
    its = 0;
    do
    {
      for (l = nn; l > 0; l--)
      {
        s = fabs(a[(l - 1) * n + l - 1]) + fabs(a[l * n + l]);
        if (s == 0.0)
          s = anorm;
        if ((fabs(a[l * n + l - 1]) + s) == s)
        {
          a[l * n + l - 1] = 0.0;
          break;
        }
      }
      x = a[nn * n + nn];
      if (l == nn)
      {
        wr[nn] = x + t;
        wi[nn--] = 0.0;
      }
      else
      {
        y = a[(nn - 1) * n + nn - 1];
        w = a[nn * n + nn - 1] * a[(nn - 1) * n + nn];
        if (l == (nn - 1))
        {
          p = 0.5 * (y - x);
          q = p * p + w;
          z = sqrt(fabs(q));
          x += t;
          if (q >= 0.0)
          {
            z = p + SIGN(z, p);
            wr[nn - 1] = wr[nn] = x + z;
            if (z)
              wr[nn] = x - w / z;
            wi[nn - 1] = wi[nn] = 0.0;
          }
          else
          {
            wr[nn - 1] = wr[nn] = x + p;
            wi[nn - 1] = -(wi[nn] = z);
          }
          nn -= 2;
        }
        else
        {
          if (its == 30)
          {
            cerr << "[ERROR][HQR] Too many iterations in hqr" << endl;
            exit(1);
          }
          if (its == 10 || its == 20)
          {
            t += x;
            for (i = 0; i < nn + 1; i++)
              a[i * n + i] -= x;
            s = fabs(a[nn * n + nn - 1]) + fabs(a[(nn - 1) * n + nn - 2]);
            y = x = 0.75 * s;
            w = -0.4375 * s * s;
          }
          ++its;
          for (m = (nn - 2); m >= l; m--)
          {
            z = a[m * n + m];
            r = x - z;
            s = y - z;
            p = (r * s - w) / a[(m + 1) * n + m] + a[m * n + m + 1];
            q = a[(m + 1) * n + m + 1] - z - r - s;
            r = a[(m + 2) * n + m + 1];
            s = fabs(p) + fabs(q) + fabs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m == l)
              break;
            u = fabs(a[m * n + m - 1]) * (fabs(q) + fabs(r));
            v = fabs(p) * (fabs(a[(m - 1) * n + m - 1]) + fabs(z) + fabs(a[(m + 1) * n + m + 1]));
            if ((u + v) == v)
              break;
          }
          for (i = m + 2; i <= nn; i++)
          {
            a[i * n + i - 2] = 0.0;
            if (i != (m + 2))
              a[i * n + i - 3] = 0.0;
          }
          for (k = m; k < nn; k++)
          {
            if (k != m)
            {
              p = a[k * n + k - 1];
              q = a[(k + 1) * n + k - 1];
              r = 0.0;
              if (k != (nn - 1))
                r = a[(k + 2) * n + k - 1];
              if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0)
              {
                p /= x;
                q /= x;
                r /= x;
              }
            }
            if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
            {
              if (k == m)
              {
                if (l != m)
                  a[k * n + k - 1] = -a[k * n + k - 1];
              }
              else
                a[k * n + k - 1] = -s * x;
              p += s;
              x = p / s;
              y = q / s;
              z = r / s;
              q /= p;
              r /= p;
              for (j = k; j < nn + 1; j++)
              {
                p = a[k * n + j] + q * a[(k + 1) * n + j];
                if (k != (nn - 1))
                {
                  p += r * a[(k + 2) * n + j];
                  a[(k + 2) * n + j] -= p * z;
                }
                a[(k + 1) * n + j] -= p * y;
                a[k * n + j] -= p * x;
              }
              mmin = nn < k + 3 ? nn : k + 3;
              for (i = l; i < mmin + 1; i++)
              {
                p = x * a[i * n + k] + y * a[i * n + k + 1];
                if (k != (nn))
                {
                  p += z * a[i * n + k + 2];
                  a[i * n + k + 2] -= p * r;
                }
                a[i * n + k + 1] -= p * q;
                a[i * n + k] -= p;
              }
            }
          }
        }
      }
    } while (l < nn - 1);
  }
}

/**
   Compute the eigenvectors of a matrix a
   dim   order of matrix
   
**/

void getvec(double *a, double *EIG, double *EIGVEC, int dim, int max_ev)
{
  int IOP;
  /** first   **/
  double *AA = new double[dim * dim];
  double *AA_inv = new double[dim * dim];
  double *BB = new double[dim * dim];
  double *CC = new double[dim * dim];
  double *AB = new double[dim * dim];
  double *BMAT = new double[dim * dim];
  double *VEC0 = new double[dim];
  double *VEC1 = new double[dim];
  double *ERR = new double[dim];

  for (IOP = 0; IOP < max_ev; ++IOP)
  {
    double ESHIFT = 0.001;
    int IOPM = IOP - 1;
    int IOPP = IOP + 1;

    if (IOPP > dim)
      IOPP = IOPM;
    if (IOPM < 1)
      IOPM = IOPP;

    double DIFM = fabs(EIG[IOPM] - EIG[IOP]) * 0.05;
    double DIFP = fabs(EIG[IOPP] - EIG[IOP]) * 0.05;

    if (ESHIFT > DIFM)
      ESHIFT = DIFM;
    if (ESHIFT > DIFP)
      ESHIFT = DIFP;

    double EIGG = EIG[IOP] + ESHIFT;
    double EIGH = EIG[IOP];

    double SUM;
    for (int J = 0; J < dim; ++J)
    {
      for (int I = 0; I < dim; ++I)
      {
        double SUM1 = a[I + J * dim];
        AA[I + J * dim] = SUM1;

        if (I == J)
          AA[I + J * dim] = AA[I + J * dim] - EIGG;

        BB[I + J * dim] = SUM1;
        CC[I + J * dim] = SUM1;
      }

    } /** end of the loop over J **/

    for (int I = 0; I < dim; ++I)
      BB[I + dim * I] = BB[I + dim * I] - EIGH;

    //	  cout << "DEBUG " << IOP << endl ;
    /** dump_matrix(AA, dim ) ; exit(0) ;   **/

    inverse(AA, AA_inv, dim);
    /** dump_matrix(AA_inv, dim ) ;   exit(0) ;   **/

    /***   *****/

    for (int J = 0; J < dim; ++J)
    {
      for (int I = 0; I < dim; ++I)
      {
        AB[I + J * dim] = 0.0;
        for (int K = 0; K < dim; ++K)
          AB[I + J * dim] += CC[I + K * dim] * AA_inv[K + J * dim];
      }
    }
    /** dump_matrix(AB, dim ) ;   exit(0) ;   **/

    /***   *****/

    double ANORM = 0.0;
    for (int I = 0; I < dim; ++I)
    {
      double RR = 0.77;
      int ISGN = +1;
      if (RR > 0.5)
        ISGN = -1;

      double VECC = (0.1 + 0.33) * ISGN;
      ANORM = ANORM + VECC * VECC;
      VEC0[I] = VECC;
    }
    ANORM = 1.0 / sqrt(ANORM);
    for (int I = 0; I < dim; ++I)
      VEC0[I] = VEC0[I] * ANORM;

    /**
	 Start inverse iteration.
      **/

    for (int ITER = 0; ITER < 12; ++ITER) /* Why is this a loop to ITER=12??*/
    {
      ANORM = 0.0;
      for (int I = 0; I < dim; ++I)
      {
        SUM = 0.0;
        for (int K = 0; K < dim; ++K)
          SUM = SUM + AA_inv[I + K * dim] * VEC0[K];

        ANORM = ANORM + SUM * SUM;
        VEC1[I] = SUM;

      } /** 13 in Fortran code **/

      ANORM = 1.0 / sqrt(ANORM);
      for (int I = 0; I < dim; ++I)
        VEC0[I] = VEC1[I] * ANORM;

    } /** end of iterations (12 fortran code) **/

    for (int I = 0; I < dim; ++I)
      EIGVEC[IOP + max_ev * I] = VEC0[I];

    ANORM = 0.0;
    for (int I = 0; I < dim; ++I)
    {
      SUM = 0.0;
      for (int K = 0; K < dim; ++K)
      {
        SUM = SUM + BB[I + K * dim] * VEC0[K];
      } /*** 22 fortran **/
      ANORM = ANORM + SUM * SUM;
    } /** 20 fortran **/
    ERR[IOP] = ANORM;

    /**
	 Compute the Rayleigh quotient

      **/

    SUM = 0.0;
    for (int I = 0; I < dim; ++I)
    {
      double VECI = VEC0[I];
      for (int K = 0; K < dim; ++K)
        SUM = SUM + CC[I + K * dim] * VEC0[K] * VECI;
    } /** 24 fortran code  ***/
      //     cout << "DEBUG Eigenvalue" << IOP << " computed " << SUM << " inputed " << EIG[IOP]  << endl ;

  } /** end of loop over IOP  ***/

  /****** free up memory  **/
  delete[] AA;
  delete[] AA_inv;
  delete[] AB;
  delete[] BB;
  delete[] CC;
  delete[] BMAT;
  delete[] VEC0;
  delete[] VEC1;
  delete[] ERR;
}

void get_eigenvalues(double *mat, int mat_size, double *evalreal)
{

  // input: jk_corr_mat_t -> matrix to diagonalize
  // output: evalreal[i] i={1,mat_size};

  double *eigi = new double[mat_size];

  double *mat_copy = new double[mat_size * mat_size];

  for (int i = 0; i < mat_size * mat_size; i++)
    mat_copy[i] = mat[i];

  elmhes(mat_copy, mat_size);
  hqr(mat_copy, mat_size, evalreal, eigi);
  quicksort(mat_size, evalreal);

  delete[] eigi;
  delete[] mat_copy;
}

// solo un'interfaccia
void get_eigenvectors(double *mat, int mat_size, int n_ev, double *evalreal, double *eigenvec)
{
  getvec(mat, evalreal, eigenvec, mat_size, n_ev);
}

#undef TINY
#undef SWAP
#undef IMAX

#undef SIGN
