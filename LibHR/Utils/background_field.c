/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *
* All rights reserved.                                                      *
\***************************************************************************/

#include "global.h"
#include "utils.h"
#include "suN.h"
#include "update.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "logger.h"
#include "hr_complex.h"
#include <math.h>
#include <stdlib.h>

#define PI 3.141592653589793238462643383279502884197

/*#define DEBUG_BACKGROUND*/

void apply_background_field_zdir(suNg_field *V, double Q, int n)
{
#if defined(GAUGE_SON) || defined(WITH_QUATERNIONS)
   error(0==0,1,"apply_background_field_zdir [background_field.c]",
         "This function cannot be used with SO(N) or quaternion representation");
#else
#ifdef DEBUG_BACKGROUND
  static suNg_field *Vtest_old = NULL;
  static suNg_field *Vtest_new = NULL;
  double complex tmp, tmp2;
  double diff_re, diff_im;

  Vtest_old = alloc_gfield(&glattice);
  Vtest_new = alloc_gfield(&glattice);
  suNg_field_copy(Vtest_old, V);

  start_gf_sendrecv(Vtest_old);
  complete_gf_sendrecv(Vtest_old);
#endif
  int index;
  double A = 0;
  double complex phase;
  suNg utmp;
  static suNg_field *V_old = NULL;
  int c[4];
  double E = 2. * PI * (double)n / (Q * (double)GLB_T * (double)GLB_Z);
  int x3, x4;

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          index = ipt(c[0], c[1], c[2], c[3]);
          x4 = c[0] + zerocoord[0];
          A = -E * (double)x4;
          phase = cos(Q * A) + I * sin(Q * A);

          utmp = *_4FIELD_AT(V, index, 3);
          _suNg_mulc(*_4FIELD_AT(V, index, 3), phase, utmp);
          if (x4 == GLB_T - 1)
          {
            x3 = c[3] + zerocoord[3];
            phase = cos(Q * E * GLB_T * x3) + I * sin(Q * E * GLB_T * x3);
            utmp = *_4FIELD_AT(V, index, 0);
            _suNg_mulc(*_4FIELD_AT(V, index, 0), phase, utmp);
          }

        } // end loop local volume

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);

  // test sum of the difference between old and new plaquette for each site. Real and imaginary parts are treated independantly
#ifdef DEBUG_BACKGROUND

  diff_re = 0.;
  diff_im = 0.;

  suNg_field_copy(Vtest_new, V);
  complete_gf_sendrecv(Vtest_new);

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {

          _complex_0(tmp);
          _complex_0(tmp2);

          index = ipt(c[0], c[1], c[2], c[3]);

          // calculate new plaquettes in direction 3,4
          cplaq(&tmp, index, 3, 0);

          // now divide by e^{iQE}
          phase = cos(Q * E) + I * sin(Q * E);
          _complex_div(tmp2, tmp, phase);

          // restore old gauge field
          suNg_field_copy(V, Vtest_old);

          start_gf_sendrecv(V);
          complete_gf_sendrecv(V);

          cplaq(&tmp, index, 3, 0);

          // compute absolute value of the difference between old and new for the real and imaginary parts independantly
          diff_re += fabs(creal(tmp2 - tmp));
          diff_im += fabs(cimag(tmp2 - tmp));

          // restore new gauge field
          suNg_field_copy(V, Vtest_new);

          start_gf_sendrecv(V);
          complete_gf_sendrecv(V);
        }
  global_sum(&diff_re, 1);
  global_sum(&diff_im, 1);

  free_gfield(Vtest_old);
  free_gfield(Vtest_new);

  lprintf("DEBUG", 0, " difference of the real and imaginary part of the old and new plaquette summed on the volume: %e %e \n", diff_re, diff_im);

#endif // DEBUG_BACKGROUND

  free_gfield(V_old);
#endif //SON or WITH_QUATERNIONS
}
