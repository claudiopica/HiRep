#include "global.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "communications.h"
#include "observables.h"
#include "string.h"

#define _MIN(x, y) ((x) > (y) ? (y) : (x))

#ifndef NDEBUG
#define MPIRET(type) int type =
#else
#define MPIRET(type)
#endif

#if defined(ROTATED_SF) || defined(BASIC_SF) || defined(BC_T_OPEN)

void polyakov()
{
  int x[4], i4d, i3d, size3d, mu;
  int loc[4] = {T, X, Y, Z};
  int mins;
  mins = _MIN(X, _MIN(Y, Z));
  suNg *p;
  suNg *lp;
  suNg *bp;
  suNg tmp;
  double complex poly[GLB_T];
  double adjpoly[GLB_T];
  double dtmp;
#ifdef WITH_MPI
  int np[4] = {NP_T, NP_X, NP_Y, NP_Z};
#endif /* WITH_MPI */

  size3d = T * X * Y * Z / mins;
  p = malloc(sizeof(suNg) * size3d);
  lp = malloc(sizeof(suNg) * size3d);
  bp = malloc(sizeof(suNg) * size3d);

  for (mu = 1; mu < 4; mu++)
  {
    size3d = T * X * Y * Z / loc[mu];
    /* local wilson lines */

    i3d = 0;
    for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
      for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
        for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
        {
          _suNg_unit(lp[i3d]);
          for (x[mu] = 0; x[mu] < loc[mu]; x[mu]++)
          {
            i4d = ipt(x[0], x[1], x[2], x[3]);
            _suNg_times_suNg(tmp, lp[i3d], *pu_gauge(i4d, mu));
            lp[i3d] = tmp;
          }
          i3d++;
        }

        /* global wilson lines */

#ifdef WITH_MPI
    if (COORD[mu] != 0)
    {
      MPI_Status status;

      MPIRET(mpiret)
      MPI_Recv((double*)(bp),                                     /* buffer */
               size3d * sizeof(suNg) / sizeof(double), /* lenght in units of doubles */
               MPI_DOUBLE,                             /* basic datatype */
               proc_dn(CID, mu),                       /* cid of origin */
               COORD[mu] - 1,                          /* tag of communication */
               cart_comm,                              /* use the cartesian communicator */
               &status);

#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        if (status.MPI_ERROR != MPI_SUCCESS)
        {
          MPI_Error_string(status.MPI_ERROR, mesg, &mesglen);
          lprintf("MPI", 0, "Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                  0,
                  status.MPI_SOURCE,
                  status.MPI_TAG,
                  mesg);
        }
        error(1, 1, "polyakov.c", "Cannot receive buf_gtf[1]");
      }
#endif /* NDEBUG */
    }
#endif /* WITH_MPI */

    if (COORD[mu] == 0)
    {
      i3d = 0;
      for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
        for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
          for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
          {
            p[i3d] = lp[i3d];
            i3d++;
          }
    }
    else
    {
      i3d = 0;
      for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
        for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
          for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
          {
            _suNg_times_suNg(p[i3d], bp[i3d], lp[i3d]);
            i3d++;
          }
    }

#ifdef WITH_MPI
    if (COORD[mu] != np[mu] - 1)
    {

      MPIRET(mpiret)
      MPI_Send((double*)(p),                                      /* buffer */
               size3d * sizeof(suNg) / sizeof(double), /* lenght in units of doubles */
               MPI_DOUBLE,                             /* basic datatype */
               proc_up(CID, mu),                       /* cid of destination */
               COORD[mu],                              /* tag of communication */
               cart_comm                               /* use the cartesian communicator */
      );

#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "polyakov.c", "Cannot send buf_gtf[0]");
      }
#endif /* NDEBUG */

      memset(p, 0, sizeof(suNg) * size3d);
    }
#endif /* WITH_MPI */

    /* trace and average */
    memset(poly, 0, sizeof(double complex) * GLB_T);
    memset(adjpoly, 0, sizeof(double) * GLB_T);

    i3d = 0;
    for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
      for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
        for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
        {

          _suNg_trace_re(dtmp, p[i3d]);
          poly[x[0] + zerocoord[0]] += dtmp;
          adjpoly[x[0] + zerocoord[0]] += dtmp * dtmp;
#ifndef GAUGE_SON
          _suNg_trace_im(dtmp, p[i3d]);
          poly[x[0] + zerocoord[0]] += I * dtmp;
          adjpoly[x[0] + zerocoord[0]] += dtmp * dtmp - 1;
#endif
          i3d++;
        }

    global_sum((double *)poly, 2 * GLB_T);
    global_sum((double *)adjpoly, GLB_T);

    for (x[0] = 0; x[0] < GLB_T; x[0]++)
    {

      lprintf("FUND_POLYAKOV", 0, "Spatial Polyakov (direction, t, (fund), (adj) )= %d %d ( %1.8e %1.8e I ) , ( %1.8e ) \n", mu, x[0],
              creal(poly[x[0]] / (GLB_VOL3 * NG) * loc[mu]),
              cimag(poly[x[0]] / (GLB_VOL3 * NG) * loc[mu]),
              adjpoly[x[0]] / (GLB_VOL3 * NG) * loc[mu]);
    }
  }
  free(p);
  free(lp);
  free(bp);
}

#else
void polyakov()
{
  int x[4], i4d, i3d, size3d;
  int loc[4] = {T, X, Y, Z};
  suNg *p;
  suNg *lp;
  suNg *bp;
  suNg tmp;
  double complex poly;
  double adjpoly;
  double dtmp;
#ifdef WITH_MPI
  int np[4] = {NP_T, NP_X, NP_Y, NP_Z};
#endif /* WITH_MPI */

  for (int mu = 0; mu < 4; mu++)
  {

    size3d = T * X * Y * Z / loc[mu];
    p = malloc(sizeof(suNg) * size3d);
    lp = malloc(sizeof(suNg) * size3d);
    bp = malloc(sizeof(suNg) * size3d);

    /* local wilson lines */

    i3d = 0;
    for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
      for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
        for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
        {
          _suNg_unit(lp[i3d]);
          for (x[mu] = 0; x[mu] < loc[mu]; x[mu]++)
          {
            i4d = ipt(x[0], x[1], x[2], x[3]);
            _suNg_times_suNg(tmp, lp[i3d], *pu_gauge(i4d, mu));
            lp[i3d] = tmp;
          }
          i3d++;
        }

        /* global wilson lines */

#ifdef WITH_MPI
    if (COORD[mu] != 0)
    {
      MPI_Status status;
      MPIRET(mpiret)
      MPI_Recv(bp,                                     /* buffer */
               size3d * sizeof(suNg) / sizeof(double), /* lenght in units of doubles */
               MPI_DOUBLE,                             /* basic datatype */
               proc_dn(CID, mu),                       /* cid of origin */
               COORD[mu] - 1,                          /* tag of communication */
               cart_comm,                              /* use the cartesian communicator */
               &status);
#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        if (status.MPI_ERROR != MPI_SUCCESS)
        {
          MPI_Error_string(status.MPI_ERROR, mesg, &mesglen);
          lprintf("MPI", 0, "Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                  0,
                  status.MPI_SOURCE,
                  status.MPI_TAG,
                  mesg);
        }
        error(1, 1, "polyakov.c", "Cannot receive buf_gtf[1]");
      }
#endif /* NDEBUG */
    }
#endif /* WITH_MPI */

    if (COORD[mu] == 0)
    {
      i3d = 0;
      for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
        for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
          for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
          {
            p[i3d] = lp[i3d];
            i3d++;
          }
    }
    else
    {
      i3d = 0;
      for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
        for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
          for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
          {
            _suNg_times_suNg(p[i3d], bp[i3d], lp[i3d]);
            i3d++;
          }
    }

#ifdef WITH_MPI
    if (COORD[mu] != np[mu] - 1)
    {
      MPIRET(mpiret)
      MPI_Send(p,                                      /* buffer */
               size3d * sizeof(suNg) / sizeof(double), /* lenght in units of doubles */
               MPI_DOUBLE,                             /* basic datatype */
               proc_up(CID, mu),                       /* cid of destination */
               COORD[mu],                              /* tag of communication */
               cart_comm                               /* use the cartesian communicator */
      );
#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "polyakov.c", "Cannot send buf_gtf[0]");
      }
#endif /* NDEBUG */
    }
#endif /* WITH_MPI */

    /* broadcast wilson lines */

#ifdef WITH_MPI

    if (COORD[mu] == np[mu] - 1)
    {
      MPI_Request comm_req[np[mu] - 1];
      int destCID = CID;
      for (int t = 0; t < np[mu] - 1; t++)
      {
        destCID = proc_up(destCID, mu);
        MPIRET(mpiret)
        MPI_Isend(p,                                      /* buffer */
                  size3d * sizeof(suNg) / sizeof(double), /* lenght in units of doubles */
                  MPI_DOUBLE,                             /* basic datatype */
                  destCID,                                /* cid of destination */
                  t,                                      /* tag of communication */
                  cart_comm,                              /* use the cartesian communicator */
                  &(comm_req[t])                          /* handle to communication request */
        );
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS)
        {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret, mesg, &mesglen);
          lprintf("MPI", 0, "ERROR: %s\n", mesg);
          error(1, 1, "polyakov.c", "Cannot start send polyakov");
        }
#endif /* NDEBUG */
      }

      MPI_Status status[np[mu] - 1];
      MPIRET(mpiret)
      MPI_Waitall(np[mu] - 1, comm_req, status);
#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen, k;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        for (k = 0; k < np[mu] - 1; ++k)
        {
          if (status[k].MPI_ERROR != MPI_SUCCESS)
          {
            MPI_Error_string(status[k].MPI_ERROR, mesg, &mesglen);
            lprintf("MPI", 0, "Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                    k,
                    status[k].MPI_SOURCE,
                    status[k].MPI_TAG,
                    mesg);
          }
        }
        error(1, 1, "polyakov.c", "Cannot complete communications");
      }
#endif /* NDEBUG */
    }
    else
    {

      int sCOORD[4], sCID;
      sCOORD[0] = COORD[0];
      sCOORD[1] = COORD[1];
      sCOORD[2] = COORD[2];
      sCOORD[3] = COORD[3];
      sCOORD[mu] = np[mu] - 1;
      MPIRET(mpiret)
      MPI_Cart_rank(cart_comm, sCOORD, &sCID);
#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "polyakov.c", "Cannot retrieve source CID");
      }
#endif /* NDEBUG */
      MPI_Status status;
      MPIRET(mpiret)
      MPI_Recv(p,                                      /* buffer */
               size3d * sizeof(suNg) / sizeof(double), /* lenght in units of doubles */
               MPI_DOUBLE,                             /* basic datatype */
               sCID,                                   /* cid of destination */
               COORD[mu],                              /* tag of communication */
               cart_comm,                              /* use the cartesian communicator */
               &status);
#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS)
      {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        if (status.MPI_ERROR != MPI_SUCCESS)
        {
          MPI_Error_string(status.MPI_ERROR, mesg, &mesglen);
          lprintf("MPI", 0, "Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                  0,
                  status.MPI_SOURCE,
                  status.MPI_TAG,
                  mesg);
        }
        error(1, 1, "polyakov.c", "Cannot receive polyakov");
      }
#endif /* NDEBUG */
    }
#endif /* WITH_MPI */

    /* trace and average */
    _complex_0(poly);
    adjpoly = 0.;

    i3d = 0;
    for (x[(mu + 1) % 4] = 0; x[(mu + 1) % 4] < loc[(mu + 1) % 4]; x[(mu + 1) % 4]++)
      for (x[(mu + 2) % 4] = 0; x[(mu + 2) % 4] < loc[(mu + 2) % 4]; x[(mu + 2) % 4]++)
        for (x[(mu + 3) % 4] = 0; x[(mu + 3) % 4] < loc[(mu + 3) % 4]; x[(mu + 3) % 4]++)
        {
          _suNg_trace_re(dtmp, p[i3d]);
          /*      lprintf("LOC_POLYAKOV",0,"%d %d %d %d %d %1.8e\n",mu,x[(mu+1)%4],x[(mu+2)%4],x[(mu+3)%4],i3d,dtmp); */
          poly += dtmp;
          adjpoly += dtmp * dtmp;
#ifndef GAUGE_SON
          _suNg_trace_im(dtmp, p[i3d]);
          poly += I * dtmp;
          adjpoly += dtmp * dtmp - 1;
#endif
          i3d++;
        }

    global_sum((double *)&poly, 2);
    global_sum((double *)&adjpoly, 1);

    poly /= NG * (GLB_VOLUME / loc[mu]);
    adjpoly /= NG * (GLB_VOLUME / loc[mu]);

    lprintf("FUND_POLYAKOV", 0, "Polyakov direction %d = %1.8e %1.8e\n", mu, creal(poly), cimag(poly));
    lprintf("ADJ_POLYAKOV", 0, "Polyakov direction %d = %1.8e\n", mu, adjpoly);

    free(p);
    free(lp);
    free(bp);
  }
}
#endif

#undef MPIRET
#undef _MIN