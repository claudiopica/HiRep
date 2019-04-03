/***************************************************************************\
 *                                    * 
 \***************************************************************************/

#include <ctype.h>
#include "global.h"
#include "logger.h"
#include "suN_utils_multilevel.h"
#include "random.h"
#include "io.h"
#include "representation.h"
#include "update.h"
#include "memory.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>

static input_pg_ml pg_var_ml = init_input_pg_ml(pg_var_ml);

static void mk_gconf_name(char *name, pg_flow_ml *gf, int id)
{
    sprintf(name, "%s_%dx%dx%dx%dnc%db%.6fn%d",
            gf->run_name, GLB_T, GLB_X, GLB_Y, GLB_Z, NG,
            gf->pg_v->beta, id);
}

static char *add_dirname(char *dirname, char *filename)
{
    static char buf[256];
    strcpy(buf, dirname);
    return strcat(buf, filename);
}

static void slower(char *str)
{
    while (*str)
    {
        *str = (char)(tolower(*str));
        ++str;
    }
}

static int parse_gstart(pg_flow_ml *gf)
{

    int t, x, y, z, ng;
    double beta;
    int ret = 0;
    char buf[256];

    ret = sscanf(gf->g_start, "%[^_]_%dx%dx%dx%dnc%db%lfn%d",
                 buf, &t, &x, &y, &z, &ng, &beta, &gf->start);

    if (ret == 8)
    { /* we have a correct file name */
        /* increase gf->start: this will be the first conf id */
        gf->start++;

        /* do some check */
        if (t != GLB_T || x != GLB_X || y != GLB_Y || z != GLB_Z)
        {
            lprintf("WARNING", 0, "Size read from config name (%d,%d,%d,%d) is different from the lattice size!\n", t, x, y, z);
        }
        if (ng != NG)
        {
            lprintf("WARNING", 0, "Gauge group read from config name (NG=%d) is not the one used in this code!\n", ng);
        }
        if (strcmp(buf, gf->run_name) != 0)
        {
            lprintf("WARNING", 0, "Run name [%s] doesn't match conf name [%s]!\n", gf->run_name, buf);
        }

        lprintf("FLOW", 0, "Starting from conf [%s]\n", gf->g_start);

        return 0;
    }

    gf->start = 1; /* reset gf->start */

    /* try if it match a unit or random start */
    strcpy(buf, gf->g_start);
    slower(buf);
    ret = strcmp(buf, "unit");
    if (ret == 0)
    {
        lprintf("FLOW", 0, "Starting a new run from a unit conf!\n");
        return 1;
    }
    ret = strcmp(buf, "random");
    if (ret == 0)
    {
        lprintf("FLOW", 0, "Starting a new run from a random conf!\n");
        return 2;
    }

    lprintf("ERROR", 0, "Invalid starting gauge conf specified [%s]\n", gf->g_start);
    error(1, 1, "parse_gstart " __FILE__, "invalid config name");

    return -1;
}

static int parse_lastconf(pg_flow_ml *gf)
{

    int ret = 0;
    int addtostart = 0;

    if (gf->last_conf[0] == '+')
    {
        addtostart = 1;
    }
    if (addtostart)
    {
        ret = sscanf(gf->last_conf, "+%d", &(gf->end));
    }
    else
    {
        ret = sscanf(gf->last_conf, "%d", &(gf->end));
    }
    if (ret == 1)
    {
        if (addtostart)
            gf->end += gf->start;
        else
            gf->end++;
        return 0;
    }

    lprintf("ERROR", 0, "Invalid last conf specified [%s]\n", gf->last_conf);
    error(1, 1, "parse_lastconf " __FILE__, "invalid last config name");

    return -1;
}

int init_mc_ml(pg_flow_ml *gf, char *ifile)
{
    int start_t;

    strcpy(gf->g_start, "invalid");
    strcpy(gf->run_name, "run_name");
    strcpy(gf->last_conf, "invalid");
    strcpy(gf->conf_dir, "./");
    gf->save_freq = 0;
    gf->therm = 0;
    gf->pg_v = &pg_var_ml;

    read_input(pg_var_ml.read, ifile);

    pg_var_ml.ml_niteration = malloc(sizeof(int) * pg_var_ml.ml_levels);
    pg_var_ml.ml_nskip = malloc(sizeof(int) * pg_var_ml.ml_levels);

    const char sep[2] = ",";
    char *token;
    /* get the first token */
    token = strtok(pg_var_ml.cml_niteration, sep);
    /* walk through other tokens */
    for (int i = 0; i < pg_var_ml.ml_levels; i++)
    {
        error(token == NULL, 1, "init_mc_ml " __FILE__, "Missing one level of number of iterartions");
        pg_var_ml.ml_niteration[i] = atoi(token);
        token = strtok(NULL, sep);
    }

    token = strtok(pg_var_ml.cml_nskip, sep);
    /* walk through other tokens */
    for (int i = 0; i < pg_var_ml.ml_levels; i++)
    {
        error(token == NULL, 1, "init_mc_ml " __FILE__, "Missing one level of skip iterartions");
        pg_var_ml.ml_nskip[i] = atoi(token);
        token = strtok(NULL, sep);
    }

    read_input(gf->read, ifile);

    /* fix conf_dir: put a / at the end of it */
    start_t = strlen(gf->conf_dir);
    if (gf->conf_dir[start_t - 1] != '/')
    {
        gf->conf_dir[start_t] = '/';
        gf->conf_dir[start_t + 1] = '\0';
    }

    /* set initial configuration */
    start_t = parse_gstart(gf);
    /* set last conf id */
    parse_lastconf(gf);

    /* init gauge field */
    switch (start_t)
    {
    case 0:
        read_gauge_field(add_dirname(gf->conf_dir, gf->g_start));
        gf->therm = 0;
        break;
    case 1:
        unit_u(u_gauge);
        break;
    case 2:
        random_u(u_gauge);
        break;
    }

    apply_BCs_on_fundamental_gauge_field();
    represent_gauge_field();

    return 0;
}

int save_conf(pg_flow_ml *gf, int id)
{
    char buf[256];

    mk_gconf_name(buf, gf, id);
#if NG == 2 && !defined(WITH_QUATERNIONS)
    write_gauge_field_su2q(add_dirname(gf->conf_dir, buf));
#else
    write_gauge_field(add_dirname(gf->conf_dir, buf));
#endif

    return 0;
}

int end_mc_ml()
{
    free_BCs();

    /* free memory */
    free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
    free_gfield_f(u_gauge_f);
#endif

    return 0;
}

#undef repr_name
