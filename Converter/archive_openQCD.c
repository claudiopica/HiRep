/*************************************************************************** \
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "suN.h"
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "communications.h"
#include "moreio.h"
#include "utils.h"
#include "observables.h"
#include "linear_algebra.h"

void read_gauge_field_openQCD(char filename[])
{
  FILE *fp = NULL;
  int g[4];
  int mu, i, j;
  struct timeval start, end, etime;
  double test[2 * NG * NG];
  int size[4];
  double readplaq;

  gettimeofday(&start, 0);

  error((fp = fopen(filename, "rb")) == NULL, 1, "read_gauge_field_openQCD",
        "Failed to open file for reading");

  error(fread_LE_int(size, 4, fp) != 4,
        1, "read_gauge_field_openQCD",
        "Failed to read lattice size from the header of the conf file");

  error(fread(&readplaq, sizeof(double), 1, fp) != 1, 1, "read_gauge_field_openQCD",
        "Failed to read the plaquette value from the header  of the conf file");

  int id;
  error(size[0] != GLB_T, 1, "read_gauge_field_openQCD", "Wrong lattice size");
  error(size[1] != GLB_X, 1, "read_gauge_field_openQCD", "Wrong lattice size");
  error(size[2] != GLB_Y, 1, "read_gauge_field_openQCD", "Wrong lattice size");
  error(size[3] != GLB_Z, 1, "read_gauge_field_openQCD", "Wrong lattice size");

  for (g[0] = 0; g[0] < GLB_T; g[0]++)
    for (g[1] = 0; g[1] < GLB_X; g[1]++)
      for (g[2] = 0; g[2] < GLB_Y; g[2]++)
        for (g[3] = 0; g[3] < GLB_Z; g[3]++)
          if ((g[0] + g[1] + g[2] + g[3]) % 2 == 1)
            for (mu = 0; mu < 4; mu++)
            {

              id = ipt(g[0], g[1], g[2], g[3]);

              error(fread(test, sizeof(double), 2 * NG * NG, fp) != 2 * NG * NG,
                    1, "read_gauge_field_openQCD",
                    "Failed to read header from file");
              for (j = 0; j < NG; j++)
                for (i = 0; i < NG; i++)
                {
                  int k = j + i * NG;
                  pu_gauge(id, mu)->c[k] = test[2 * k] + I * test[2 * k + 1];
                }

              id = idn(id, mu);

              error(fread(test, sizeof(double), 2 * NG * NG, fp) != 2 * NG * NG,
                    1, "read_gauge_field_openQCD",
                    "Failed to read header from file");
              for (j = 0; j < NG; j++)
                for (i = 0; i < NG; i++)
                {
                  int k = j + i * NG;
                  pu_gauge(id, mu)->c[k] = test[2 * k] + I * test[2 * k + 1];
                }
            }

  fclose(fp);

  double new_plaq = avr_plaquette();

  if (sqrt((NG * new_plaq - readplaq) * (NG * new_plaq - readplaq)) > 1.e-14)
    error(1, 1, "read_gauge_field_openQCD", "Wrong plaquette checksum");
  else
    lprintf("IO", 0, "Plaquette checksum matches\n Initial plaquette: %1.8e \n", new_plaq);

  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("IO", 0, "Configuration [%s] read [%ld sec %ld usec]\n", filename, etime.tv_sec, etime.tv_usec);
}

void read_gauge_field_openQCD_SF(char filename[])
{
  FILE *fp = NULL;
  int g[4];
  int mu, i, j;
  struct timeval start, end, etime;
  double test[2 * NG * NG];
  int size[4];
  double readplaq;

  gettimeofday(&start, 0);

  error((fp = fopen(filename, "rb")) == NULL, 1, "read_gauge_field_openQCD",
        "Failed to open file for reading");

  error(fread_LE_int(size, 4, fp) != 4,
        1, "read_gauge_field_openQCD",
        "Failed to read lattice size from the header of the conf file");

  error(fread(&readplaq, sizeof(double), 1, fp) != 1, 1, "read_gauge_field_openQCD",
        "Failed to read the plaquette value from the header  of the conf file");

  int id;
  error(size[0] + 2 != GLB_T, 1, "read_gauge_field_openQCD", "Wrong lattice size");
  error(size[1] != GLB_X, 1, "read_gauge_field_openQCD", "Wrong lattice size");
  error(size[2] != GLB_Y, 1, "read_gauge_field_openQCD", "Wrong lattice size");
  error(size[3] != GLB_Z, 1, "read_gauge_field_openQCD", "Wrong lattice size");

  for (g[0] = 0; g[0] < GLB_T - 2; g[0]++)
    for (g[1] = 0; g[1] < GLB_X; g[1]++)
      for (g[2] = 0; g[2] < GLB_Y; g[2]++)
        for (g[3] = 0; g[3] < GLB_Z; g[3]++)
          if ((g[0] + g[1] + g[2] + g[3]) % 2 == 1)
            for (mu = 0; mu < 4; mu++)
            {

              id = ipt(g[0] + 1, g[1], g[2], g[3]);
              error(fread(test, sizeof(double), 2 * NG * NG, fp) != 2 * NG * NG,
                    1, "read_gauge_field_openQCD",
                    "Failed to read header from file");
              for (j = 0; j < NG; j++)
                for (i = 0; i < NG; i++)
                {
                  int k = j + i * NG;
                  pu_gauge(id, mu)->c[k] = test[2 * k] + I * test[2 * k + 1];
                }

              if (g[0] == 0)
              {
                id = ipt(GLB_T - 1, g[1], g[2], g[3]);
              }
              id = idn(id, mu);

              error(fread(test, sizeof(double), 2 * NG * NG, fp) != 2 * NG * NG,
                    1, "read_gauge_field_openQCD",
                    "Failed to read header from file");
              for (j = 0; j < NG; j++)
                for (i = 0; i < NG; i++)
                {
                  int k = j + i * NG;
                  pu_gauge(id, mu)->c[k] = test[2 * k] + I * test[2 * k + 1];
                }
            }

  fclose(fp);

  double new_plaq = avr_plaquette();

  if (sqrt((NG * new_plaq - readplaq) * (NG * new_plaq - readplaq)) > 1.e-14)
    lprintf("WARNING", 0, " Plaquette DOES NOT match!  \n");
  else
    lprintf("IO", 0, "Plaquette checksum matches\n Initial plaquette: %1.8e \n", new_plaq);

  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("IO", 0, "Configuration [%s] read [%ld sec %ld usec]\n", filename, etime.tv_sec, etime.tv_usec);
}

void write_gauge_field_openQCD(char filename[])
{
  FILE *fp = NULL;
  int g[4];
  int mu, i, j;
  struct timeval start, end, etime;
  double test[2 * NG * NG];
  int size[4] = {GLB_T, GLB_X, GLB_Y, GLB_Z};
  double writeplaq = NG * avr_plaquette();

  gettimeofday(&start, 0);

  error((fp = fopen(filename, "wb")) == NULL, 1, "write_gauge_field_openQCD",
        "Failed to open file for writing");

  error(fwrite_LE_int(size, 4, fp) != 4,
        1, "write_gauge_field_openQCD",
        "Failed to write lattice size into the header of the conf file");

  error(fwrite_LE_double(&writeplaq, 1, fp) != 1, 1, "write_gauge_field_openQCD",
        "Failed to write the plaquette value into the header of the conf file");

  int id;

  for (g[0] = 0; g[0] < GLB_T; g[0]++)
    for (g[1] = 0; g[1] < GLB_X; g[1]++)
      for (g[2] = 0; g[2] < GLB_Y; g[2]++)
        for (g[3] = 0; g[3] < GLB_Z; g[3]++)
          if ((g[0] + g[1] + g[2] + g[3]) % 2 == 1)
            for (mu = 0; mu < 4; mu++)
            {

              id = ipt(g[0], g[1], g[2], g[3]);

              for (j = 0; j < NG; j++)
                for (i = 0; i < NG; i++)
                {
                  int k = j + i * NG;
                  test[2 * k] = creal(pu_gauge(id, mu)->c[k]);
                  test[2 * k + 1] = cimag(pu_gauge(id, mu)->c[k]);
                }

              error(fwrite_LE_double(test, 2 * NG * NG, fp) != 2 * NG * NG,
                    1, "read_gauge_field_openQCD",
                    "Failed to read header from file");

              id = idn(id, mu);

              for (j = 0; j < NG; j++)
                for (i = 0; i < NG; i++)
                {
                  int k = j + i * NG;
                  test[2 * k] = creal(pu_gauge(id, mu)->c[k]);
                  test[2 * k + 1] = cimag(pu_gauge(id, mu)->c[k]);
                }
              error(fwrite_LE_double(test, 2 * NG * NG, fp) != 2 * NG * NG,
                    1, "read_gauge_field_openQCD",
                    "Failed to read header from file");
            }

  fclose(fp);

  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("IO", 0, "Plaquette of the stored configuration: %1.8e\n", writeplaq);
  lprintf("IO", 0, "Configuration [%s] wrote [%ld sec %ld usec]\n", filename, etime.tv_sec, etime.tv_usec);
}

static void write_gauge_field_hirep(char filename[], double subs)
{
  FILE *fp = NULL;
  int g[4], p[4];
  double *buff = NULL;
  int zsize, rz;
  double plaq;
  struct timeval start, end, etime;

#ifndef ALLOCATE_REPR_GAUGE_FIELD
  complete_gf_sendrecv(u_gauge);
  apply_BCs_on_represented_gauge_field(); //Save the link variables with periodic boundary conditions
#endif

  plaq = avr_plaquette(); /* to use as a checksum in the header */

  double plaqt[GLB_T], plaqs[GLB_T];
  avr_plaquette_time(plaqt, plaqs);
  plaq = 0.;
  for (int kk = 1; kk < GLB_T - 1; kk++)
  {
    plaq += plaqt[kk] + plaqs[kk];
  }
  plaq += plaqt[0] + (plaqs[GLB_T - 1] + plaqs[0]) / 2.0;
  plaq /= 2.0 * (GLB_T - subs);

  int d[5] = {NG, GLB_T, GLB_X, GLB_Y, GLB_Z};
  error((fp = fopen(filename, "wb")) == NULL, 1, "write_gauge_field",
        "Failed to open file for writing");
  /* write NG and global size */
  error(fwrite_BE_int(d, (size_t)(5), fp) != (5),
        1, "write_gauge_field",
        "Failed to write gauge field geometry");
  /* write average plaquette */
  error(fwrite_BE_double(&plaq, (size_t)(1), fp) != (1),
        1, "write_gauge_field",
        "Failed to write gauge field plaquette");

  gettimeofday(&start, 0);

  zsize = GLB_Z / NP_Z;
  rz = GLB_Z - zsize * NP_Z;

  buff = malloc(sizeof(suNg) * 4 * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));

  g[3] = 0;
  for (g[0] = 0; g[0] < GLB_T; ++g[0])
  { /* loop over T, X and Y direction */
    for (g[1] = 0; g[1] < GLB_X; ++g[1])
    {
      for (g[2] = 0; g[2] < GLB_Y; ++g[2])
      {
        for (p[3] = 0; p[3] < NP_Z; ++p[3])
        { /* loop over processors in Z direction */
          int bsize;

          bsize = sizeof(suNg) / sizeof(double) * 4 * (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */

          /* fill link buffer */
          int lsite[4];
          suNg *cm;

          /* convert global to local coordinate */
          origin_coord(lsite);
          lsite[0] = g[0] - lsite[0];
          lsite[1] = g[1] - lsite[1];
          lsite[2] = g[2] - lsite[2];

          /* fill buffer */
          cm = (suNg *)buff;
          for (lsite[3] = 0; lsite[3] < Z; ++lsite[3])
          { /* loop on local Z */
            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
            suNg *pm = pu_gauge(ix, 0);
            *(cm++) = *(pm++); /* copy 4 directions */
            *(cm++) = *(pm++);
            *(cm++) = *(pm++);
            *(cm++) = *(pm);
          }

          /* write buffer to file */
          error(fwrite_BE_double(buff, (size_t)(bsize), fp) != (bsize),
                1, "write_gauge_field",
                "Failed to write gauge field to file");

        } /* end loop over processors in Z direction */
      }
    }
  } /* end loop over T, X and Y direction */

  fclose(fp);
  free(buff);

  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("IO", 0, "Configuration [%s] saved [%ld sec %ld usec]\n", filename, etime.tv_sec, etime.tv_usec);

#ifndef ALLOCATE_REPR_GAUGE_FIELD
  complete_gf_sendrecv(u_gauge);
  apply_BCs_on_represented_gauge_field(); //Restore the right boundary conditions
#endif
}

void write_gauge_field_hirep_pbc_to_obc(char filename[])
{
  write_gauge_field_hirep(filename,1.0);
}

void write_gauge_field_hirep_pbc_to_sf(char filename[])
{
  write_gauge_field_hirep(filename,2.0);
}