/***************************************************************************\
* Copyright (c)                                                             *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "utils.h"
#include "suN.h"
#include "memory.h"
#include "communications.h"
#include "logger.h"
#include "hr_complex.h"
#include "random.h"
#include <math.h>
#include <stdlib.h>

#define gauge_field_pointer(u, ix, mu) ((u->ptr) + coord_to_index(ix, mu))

static suNg_field *_g_1 = NULL;
static suNg_field *_g_2 = NULL;
static suNg_field *_g = NULL;

static int *_iup_1;
static int *_idn_1;
static int *_iup_2;
static int *_idn_2;
static int *_iup;
static int *_idn;
static int *_nnpointer_mem;

static int blocking_level;
static int *blocking_list;
static int nblocking;

int iup_sblk(int site, int dir)
{
    return _iup[(site)*4 + (dir)];
}

int idn_sblk(int site, int dir)
{
    return _idn[(site)*4 + (dir)];
}

suNg *pu_gauge_sblk(int site, int dir)
{
    return ((_g->ptr) + coord_to_index(site, dir));
}

void initialize_spatial_blocking(int *tlist)
{
    if (NP_X != 1 || NP_Y != 1 || NP_Z != 1)
        error(1, 1, "initialize_spatial_blocking [spatial_blocking.c]",
              "Error the blocking module can be used on with periodic spatial boundary conditions");

    if (_g_1 == NULL)
    {

        _g = u_gauge;
        _iup = iup;
        _idn = idn;

        int t;
        lprintf("SPATIAL BLOCKING", 0, "Initializing memory allocation for Spatial Blocking\n");

        _g_1 = alloc_gfield(&glattice);
        error((_g_1 == NULL), 1, "initialize_spatial_blocking [spatial_blocking.c]",
              "Cannot allocate memory");

        unit_u(_g_1);

        _g_2 = alloc_gfield(&glattice);
        error((_g_2 == NULL), 1, "initialize_spatial_blocking [spatial_blocking.c]",
              "Cannot allocate memory");

        unit_u(_g_2);

        size_t req_mem = 4 * glattice.gsize_gauge; /* for each _iup_1/out and _idn_1/out */

        _nnpointer_mem = malloc(4 * req_mem * sizeof(int));
        error((_nnpointer_mem == NULL), 1, "initialize_spatial_blocking [spatial_blocking.c]",
              "Cannot allocate memory");
        _iup_1 = _nnpointer_mem;
        _iup_2 = _iup_1 + req_mem;
        _idn_1 = _iup_2 + req_mem;
        _idn_2 = _idn_1 + req_mem;
        if (tlist != NULL)
        {
            nblocking = 0;
            for (t = 0; t < GLB_T; t++)
            {
                if (tlist[t] != 0 && t >= zerocoord[0] && t < zerocoord[0] + T)
                    nblocking++;
            }
            blocking_list = malloc(nblocking * sizeof(int));
            nblocking = 0;

            for (t = 0; t < GLB_T; t++)
            {
                if (tlist[t] != 0 && t >= zerocoord[0] && t < zerocoord[0] + T)
                {
                    blocking_list[nblocking] = t - zerocoord[0];
                    nblocking++;
                }
            }
            lprintf("SPATIAL BLOCKING", 0, "Partial Spatial lattice blocking (%d t-slices per proc)\n", nblocking);
        }
        else
        {
            lprintf("SPATIAL BLOCKING", 0, "Full Spatial lattice blocking\n");
            blocking_list = malloc(T * sizeof(int));
            nblocking = T;

            for (t = 0; t < T; t++)
                blocking_list[t] = t;
        }
    }
    else
        lprintf("SPATIAL BLOCKING", 0, "Already initialized\n");
}

void free_spatial_blocking()
{
    if (_g_1 != NULL)
    {
        free_gfield(_g_1);

        free_gfield(_g_2);

        free(_nnpointer_mem);

        free(blocking_list);
    }
}

static void slice_spatial_blocking(int t, suNg_field *g_in, int *iup_in, int *idn_in, suNg_field *g_out, int *iup_out, int *idn_out)
{
    int ix, iy, iz, mu, mid, midn;
    suNg *u, *w, *v;

    for (ix = 0; ix < X; ix++)
        for (iy = 0; iy < Y; iy++)
            for (iz = 0; iz < Z; iz++)
            {
                mid = ipt(t, ix, iy, iz);
                for (mu = 1; mu < 4; mu++)
                {

                    u = gauge_field_pointer(g_out, mid, mu);

                    v = gauge_field_pointer(g_in, mid, mu);

                    midn = iup_in[4 * mid + mu];
                    iup_out[(mid)*4 + (mu)] = iup_in[midn * 4 + mu];
                    w = gauge_field_pointer(g_in, midn, mu);

                    _suNg_times_suNg(*u, *v, *w);

                    midn = idn_in[4 * mid + mu];
                    idn_out[4 * mid + mu] = idn_in[4 * midn + mu];
                }
            }
}

void spatial_blocking(unsigned int level, int force_eval)
{
    int t;

    error((_g_1 == NULL), 1, "spatial_blocking [spatial_blocking.c]", "Uninitialized function, run initialize_spatial_blocking");

    if (force_eval || blocking_level > level)
    {
        _g = u_gauge;
        _iup = iup;
        _idn = idn;
        blocking_level = 0;
    }

    if (blocking_level == level)
        return;

    if (level != blocking_level + 1)
        spatial_blocking(level - 1, 0);

    for (t = 0; t < nblocking; t++)
        slice_spatial_blocking(blocking_list[t], _g, _iup, _idn, _g_1, _iup_1, _idn_1);

    _g = _g_1;
    _iup = _iup_1;
    _idn = _idn_1;

    _g_1 = _g_2;
    _iup_1 = _iup_2;
    _idn_1 = _idn_2;

    _g_2 = _g;
    _iup_2 = _iup;
    _idn_2 = _idn;

    blocking_level++;
}

#undef gauge_field_pointer