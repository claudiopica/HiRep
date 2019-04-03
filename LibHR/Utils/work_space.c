#include "global.h"
#include "suN.h"
#include "memory.h"
#include <stdlib.h>
#include <string.h>
#include "logger.h"

static suNg_field **_g_wrk = NULL;
static int **_iup_wrk;
static int **_idn_wrk;

suNg_field *_g = NULL;
int *_iup;
int *_idn;

static int *_wrk_reserved = NULL;
static int n_alloc = 0;
static int n_reserved = 0;

int iup_wrk(int site, int dir)
{
    return _iup[(site)*4 + (dir)];
}

int idn_wrk(int site, int dir)
{
    return _idn[(site)*4 + (dir)];
}

suNg *pu_gauge_wrk(int site, int dir)
{
    return ((_g->ptr) + coord_to_index(site, dir));
}

suNg_field *u_gauge_wrk()
{
    return _g;
}

void reset_wrk_pointers()
{
    _iup = iup;
    _idn = idn;
    _g = u_gauge;
}

void set_wrk_space(int i)
{
    if (_wrk_reserved[i])
    {
        _iup = _iup_wrk[i];
        _idn = _idn_wrk[i];
        _g = _g_wrk[i];
    }
    else
        error(0, 1, "set_wrk_space", "Invalid work pointer index");
}

void set_wrk_space_and_pointers(int i, suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out)
{
    if (_wrk_reserved[i])
    {
        _iup = _iup_wrk[i];
        _idn = _idn_wrk[i];
        _g = _g_wrk[i];
        *g_wrk_out = _g_wrk[i];
        *i_up_wrk_out = _iup_wrk[i];
        *i_dn_wrk_out = _idn_wrk[i];
    }
    else
        error(0, 1, "set_wrk_space_and_pointers", "Invalid work pointer index");
}

int reserve_wrk_space()
{
    int j;
    if (n_alloc - n_reserved == 0)
    {
        suNg_field **_g_wrk_new = malloc(sizeof(suNg_field *) * (n_alloc + 1));
        int **_iup_wrk_new = malloc(sizeof(int *) * (n_alloc + 1));
        int **_idn_wrk_new = malloc(sizeof(int *) * (n_alloc + 1));
        int *_wrk_reserved_new = malloc(sizeof(int) * (n_alloc + 1));

        if (n_alloc != 0)
        {
            memcpy(_g_wrk_new, _g_wrk, sizeof(suNg_field *) * n_alloc);
            memcpy(_iup_wrk_new, _iup_wrk, sizeof(int *) * n_alloc);
            memcpy(_idn_wrk_new, _idn_wrk, sizeof(int *) * n_alloc);
            memcpy(_wrk_reserved_new, _wrk_reserved, sizeof(int) * n_alloc);
            free(_g_wrk);
            free(_iup_wrk);
            free(_idn_wrk);
            free(_wrk_reserved);
        }

        _g_wrk = _g_wrk_new;
        _iup_wrk = _iup_wrk_new;
        _idn_wrk = _idn_wrk_new;
        _wrk_reserved = _wrk_reserved_new;

        _g_wrk[n_alloc] = alloc_gfield(&glattice);

        error((_g_wrk[n_alloc] == NULL), 1, "reserve_wrk_space [work_space.c]",
              "Cannot allocate memory");

        size_t req_mem = 4 * glattice.gsize_gauge;

        _iup_wrk[n_alloc] = malloc(2 * req_mem * sizeof(int));
        error((_iup_wrk[n_alloc] == NULL), 1, "reserve_wrk_space [work_space.c]",
              "Cannot allocate memory");
        _idn_wrk[n_alloc] = _iup_wrk[n_alloc] + req_mem;

        memcpy(_iup_wrk[n_alloc], iup, req_mem * sizeof(int));
        memcpy(_idn_wrk[n_alloc], idn, req_mem * sizeof(int));

        _wrk_reserved[n_alloc] = 0;

        n_alloc++;
    }

    j = 0;

    for (j = 0; j < n_alloc; j++)
    {
        if (_wrk_reserved[j] == 0)
            break;
    }

    _wrk_reserved[j] = 1;
    _iup = _iup_wrk[j];
    _idn = _idn_wrk[j];
    _g = _g_wrk[j];

    n_reserved++;

    return j;
}

int reserve_wrk_space_with_pointers(suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out)
{
    int j = reserve_wrk_space();
    *g_wrk_out = _g_wrk[j];
    *i_up_wrk_out = _iup_wrk[j];
    *i_dn_wrk_out = _idn_wrk[j];
    return j;
}

void release_wrk_space(int id_release)
{
    _wrk_reserved[id_release] = 0;
    n_reserved--;
}

void free_wrk_space()
{
    int j;
    if (n_alloc != 0)
    {
        for (j = 0; j < n_alloc; j++)
        {
            free_gfield(_g_wrk[j]);
            free(_iup_wrk[j]);
        }

        free(_g_wrk);
        free(_iup_wrk);
        free(_idn_wrk);
        free(_wrk_reserved);
        n_reserved = n_alloc = 0;
    }
}
