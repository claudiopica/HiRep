/*
* NOCOMPILE= BASIC_SF
* NOCOMPILE= ROTATED_SF
* NOCOMPILE= FERMION_THETA
*/
/*******************************************************************************
*
* Checks of propagator, spinmatrix and the sequential sources
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gaugefix.h"
#include "spectrum.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "propagator.h"
#include "setup.h"

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

/* Mesons parameters */
typedef struct _input_mesons
{
  char mstring[256];
  double precision;
  int ti;
  int tf;
  int ff_fixed_point;
  int dt;
  /* for the reading function */
  input_record_t read[8];
} input_mesons;

input_mesons mes_var;

//Copied From meson_measurements
static void fix_T_bc(int tau)
{
  int index;
  int ix, iy, iz;
  suNf *u;
  if (--tau < 0)
    tau += GLB_T;
  lprintf("meson_measurements", 15, "Setting Dirichlet boundary conidtion at global time slice %d, %d\n", tau, T_BORDER);
  if ((zerocoord[0] - 1 <= tau && zerocoord[0] + T > tau) || (zerocoord[0] == 0 && tau == GLB_T - 1))
  {
    for (ix = 0; ix < X_EXT; ++ix)
      for (iy = 0; iy < Y_EXT; ++iy)
        for (iz = 0; iz < Z_EXT; ++iz)
        {
          if (((tau == zerocoord[0] - 1) || (zerocoord[0] == 0 && tau == GLB_T - 1)) && (NP_T > 1))
          {
            //printf("PID = %d, zc = %d, tau = %d\n", PID, zerocoord[0], tau);
            index = ipt_ext(0, ix, iy, iz);
          }
          else
          {
            index = ipt_ext(T_BORDER + tau - zerocoord[0], ix, iy, iz);
          }
          if (index != -1)
          {
            u = pu_gauge_f(index, 0);
            _suNf_zero(*u);
          }
        }
  }
  lprintf("meson_measurements", 50, "Boundaries set!\n");
}

//Check: g5 D^dag(x,0) g5 = D(0,x)
static int check_g5herm(spinor_field *prop1, int t1, spinor_field *prop2)
{

  int beta, a, ix1, ix2;

  suNf_propagator sp1, sp2, spdag;
  double complex tr = 0.;
#ifdef WITH_MPI
  MPI_Status st;
#endif

  int sender_rank = 0;
  if (t1 > zerocoord[0] && t1 < zerocoord[0] + T && COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0)
  {
    sender_rank = PID;
    ix1 = ipt(t1 - zerocoord[0], 0, 0, 0);
    for (a = 0; a < NF; ++a)
      for (beta = 0; beta < 4; beta++)
        _propagator_assign(sp1, *_FIELD_AT(&prop1[a * 4 + beta], ix1), a, beta); //S( (ti,0,0,0), (0,0,0,0) )
#ifdef WITH_MPI
    MPI_Send(&sp1, NF * 4 * 4 * sizeof(suNf_vector) / sizeof(double), MPI_DOUBLE, 0, 999, GLB_COMM);
#endif
  }

  global_sum_int(&sender_rank, 1);

  if (COORD[0] == 0 && COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0)
  {
#ifdef WITH_MPI
    MPI_Recv(&sp1, NF * 4 * 4 * sizeof(suNf_vector) / sizeof(double), MPI_DOUBLE, sender_rank, 999, GLB_COMM, &st);
#endif

    ix2 = ipt(0, 0, 0, 0);

    for (a = 0; a < NF; ++a)
    {
      for (beta = 0; beta < 4; beta++)
      {
        _propagator_assign(sp2, *_FIELD_AT(&prop2[a * 4 + beta], ix2), a, beta); //S( (0,0,0,0), (ti,0,0) )
      }
    }
    _propagator_dagger(spdag, sp2);
    _g5_propagator(sp2, spdag);
    _propagator_g5(spdag, sp2);
    lprintf("CK_G5HERM", 0, "Propagator1\n");
    //print_prop(sp1);
    lprintf("CK_G5HERM", 0, "g5 Propagator2^dagger g5 \n");
    //print_prop(spdag);
    _propagator_sub(sp2, sp1, spdag);
    lprintf("CK_G5HERM", 0, "Propagator1 - g5 Propagator2^dagger g5 \n");
    //print_prop(sp2);
    _propagator_trace(tr, sp2);
    lprintf("CK_G5HERM", 0, "Tr[ g5 Propagator1^dag g5 - Propagator2 ] = %g + I%g\n", creal(tr), cimag(tr));
  }

  if (cabs(tr) > 1.e-14)
    return 1;
  else
    return 0;
}

//source = g5 prop( x, 0 ) delta( x, (tf,0,0,0) )
void create_sequential_source_point(spinor_field *source, int tf, spinor_field *prop)
{

  int beta, a, ix;

  suNf_propagator sp0, sp1;

  for (beta = 0; beta < 4 * NF; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }

  ix = ipt(tf - zerocoord[0], 0, 0, 0);
  for (a = 0; a < NF; ++a)
  {
    for (beta = 0; beta < 4; beta++)
    {
      _propagator_assign(sp0, *_FIELD_AT(&prop[a * 4 + beta], ix), a, beta);
    }
  }
  _g5_propagator(sp1, sp0); //g5 Prop
  _propagator_transpose(sp0, sp1);

  for (a = 0; a < NF; ++a)
  {
    for (beta = 0; beta < 4; beta++)
    {
      *_FIELD_AT(&source[a * 4 + beta], ix) = sp0.c[a].c[beta];
    }
  }

  for (beta = 0; beta < 4 * NF; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}
//Check sequential Gamma Seq(0,0) = Gamma g5 S^dag(x,0) g5 g5 S(x,0)
static int check_sequential_point(spinor_field *prop_1, spinor_field *prop_2, spinor_field *prop_seq, int ti)
{

  lprintf("CK_SEQ", 0, "Only Works in serial!\n");
  int ix1 = ipt(0, 0, 0, 0);
  int ix2 = ipt(ti, 0, 0, 0);
  suNf_propagator sp1, sp2, sp3, sptmp1, sptmp2;
  int a, beta;
  double complex tr;

  for (a = 0; a < NF; ++a)
  {
    for (beta = 0; beta < 4; beta++)
    {
      _propagator_assign(sp1, *_FIELD_AT(&prop_seq[a * 4 + beta], ix1), a, beta); //S( (0,0,0,0), (2,0,0,0) ) g5 S( (2,0,0,0), (0,0,0,0) )
      _propagator_assign(sp2, *_FIELD_AT(&prop_2[a * 4 + beta], ix1), a, beta);   //S( (0,0,0,0), (2,0,0,0) )
      _propagator_assign(sp3, *_FIELD_AT(&prop_1[a * 4 + beta], ix2), a, beta);   //S( (2,0,0,0), (0,0,0,0) )
    }
  }

  _g5_propagator(sptmp1, sp3);
  _propagator_mul(sptmp2, sp2, sptmp1);

  _propagator_trace(tr, sptmp2);
  lprintf("CK_SEQ", 0, "S(0,x) g5 S(x,0) point = %g + I%g\n", creal(tr), cimag(tr));
  //print_prop(sptmp2);

  _propagator_trace(tr, sp1);
  lprintf("CK_SEQ", 0, "S(0,x) g5 S(x,0) seq   = %g + I%g\n", creal(tr), cimag(tr));
  //print_prop(sp1);

  _propagator_sub(sptmp1, sp1, sptmp2);
  _propagator_trace(tr, sptmp1);
  lprintf("CK_SEQ", 0, "point - seq   = %g + I%g\n", creal(tr), cimag(tr));
  //print_prop(sptmp1);

  if (cabs(tr) > 1.e-14)
    return 1;
  else
    return 0;
}

static int check_sequential(spinor_field *prop_seq, spinor_field *prop_1, int tau, int tf)
{

  lprintf("CK_SEQ", 0, "Only Works in serial!\n");

  double complex Corr[2][GLB_T];
  int ix, t, x, y, z, a, beta, tc;
  double complex tr;
  suNf_propagator sp0, sp1, spdag, sptmp;

  int i;
  for (i = 0; i < 2; i++)
  {
    for (t = 0; t < GLB_T; t++)
    {
      Corr[i][t] = 0;
    }
  }

  for (t = 0; t < T; t++)
  {
    tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
    for (x = 0; x < X; x++)
      for (y = 0; y < Y; y++)
        for (z = 0; z < Z; z++)
        {
          ix = ipt(t, x, y, z);

          for (a = 0; a < NF; ++a)
          {
            for (beta = 0; beta < 4; beta++)
            {
              _propagator_assign(sp0, *_FIELD_AT(&prop_seq[a * 4 + beta], ix), a, beta);
              _propagator_assign(sp1, *_FIELD_AT(&prop_1[a * 4 + beta], ix), a, beta);
            }
          }
          _propagator_dagger(spdag, sp1);

          //Pion
          _propagator_mul(sptmp, sp1, spdag);
          _propagator_trace(tr, sptmp);
          Corr[0][tc] += tr;

          //g5 Seq(0,0)
          if (t == 0 && x == 0 && y == 0 && z == 0)
          {
            _g5_propagator(sptmp, sp0);
            _propagator_trace(tr, sptmp);
            Corr[1][tc] += creal(tr);
          }
        } //END SPATIAL LOOP
  }       //END TIME LOOP

  lprintf("CORR", 0, "Pion: ");
  for (t = 0; t < T; t++)
  {
    lprintf("CORR", 0, "%g ", creal(Corr[0][t]) / (GLB_VOL3));
  }
  lprintf("CORR", 0, "\n");
  lprintf("CORR", 0, "SeqPion: ");
  for (t = 0; t < T; t++)
  {
    lprintf("CORR", 0, "%g ", creal(Corr[1][t] / (GLB_VOL3)));
  }
  lprintf("CORR", 0, "\n");

  if (cabs(Corr[0][tf] - Corr[1][0]) > 1e-14)
    return 1;
  else
    return 0;
}

int main(int argc, char *argv[])
{
  int k, tmp;
  int return_value = 0;
  double m[256];
  struct timeval start, end, etime;
  /* setup process id and communications */

  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  spinor_field *source;
  spinor_field *prop_1;
  spinor_field *prop_2;
  spinor_field *source_seq;
  spinor_field *prop_seq;

  source = alloc_spinor_field_f(4, &glattice);
  source_seq = alloc_spinor_field_f(4 * NF, &glattice);
  prop_1 = alloc_spinor_field_f(4 * NF, &glattice);
  prop_2 = alloc_spinor_field_f(4 * NF, &glattice);
  prop_seq = alloc_spinor_field_f(4 * NF, &glattice);

  for (k = 0; k < 4 * NF; k++)
  {
    spinor_field_zero_f(prop_1 + k);
    spinor_field_zero_f(prop_2 + k);
    spinor_field_zero_f(prop_seq + k);
  }

  mes_var.precision = 1e-24;
  mes_var.ti = GLB_T / 4;
  mes_var.tf = 3 * GLB_T / 4;
  m[0] = 20.12;

  lprintf("MAIN", 0, "ti = %d\n", mes_var.ti);
  lprintf("MAIN", 0, "tf = %d\n", mes_var.tf);
  lprintf("MAIN", 0, "m = %f\n", m[0]);
  mes_var.ff_fixed_point = 1;
  lprintf("MAIN", 0, "Inverter precision = %e\n", mes_var.precision);

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  represent_gauge_field();

  gettimeofday(&start, 0);

  init_propagator_eo(1, m, mes_var.precision); //1 for number of masses
  for (k = 0; k < NF; ++k)
  {
    create_point_source(source, 0, k);
    calc_propagator(prop_1 + 4 * k, source, 4); //4 for spin components
    create_point_source(source, mes_var.ti, k);
    calc_propagator(prop_2 + 4 * k, source, 4); //4 for spin components
  }
  tmp = check_g5herm(prop_1, mes_var.ti, prop_2);
  return_value += tmp;
  //  create_sequential_source_point(source_seq, mes_var.ti, prop_1);
  //  calc_propagator(prop_seq, source_seq, 4 * NF);
  //  tmp = check_sequential_point(prop_1, prop_2, prop_seq, mes_var.ti);
  //  return_value += tmp;
  //  create_sequential_source(source_seq, mes_var.tf, prop_1);
  //  calc_propagator(prop_seq, source_seq, 4 * NF);
  //  tmp = check_sequential(prop_seq, prop_1, 0, mes_var.tf);
  //  return_value += tmp;
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("MAIN", 0, "Random configuration analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

  free_propagator_eo();
  free_spinor_field_f(source);
  free_spinor_field_f(source_seq);
  free_spinor_field_f(prop_1);
  free_spinor_field_f(prop_2);
  free_spinor_field_f(prop_seq);

  lprintf("MAIN", 0, "return_value = %d\n", return_value);
  finalize_process();

  return return_value;
}
