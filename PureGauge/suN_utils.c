/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "logger.h"
#include "suN_utils.h"
#include "random.h"
#include "io.h"
#include "representation.h"
#include "update.h"
#include "memory.h"
#include "utils.h"

#include <stdio.h>
#include <string.h>


input_pg pg_var = init_input_pg(pg_var);

/* short 3-letter name to use in gconf name */
#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif

static void mk_gconf_name(char *name, pg_flow *gf, int id) {
  sprintf(name,"%s_%dx%dx%dx%dnc%db%.6fn%d",
           gf->run_name,GLB_T,GLB_X,GLB_Y,GLB_Z,NG,
           gf->pg_v->beta,id);
}

static char *add_dirname(char *dirname, char *filename) {
  static char buf[256];
  strcpy(buf,dirname);
  return strcat(buf,filename);
}

#include <ctype.h>

static void slower(char *str) {
  while (*str) {
    *str=(char)(tolower(*str));
    ++str;
  }
}

static int parse_gstart(pg_flow *gf) {

  int t, x, y, z, ng;
  double beta;
  int ret=0;
  char buf[256];

  ret=sscanf(gf->g_start,"%[^_]_%dx%dx%dx%dnc%db%lfn%d",
      buf,&t,&x,&y,&z,&ng,&beta,&gf->start);

  if(ret==8) { /* we have a correct file name */
    /* increase gf->start: this will be the first conf id */
    gf->start++;

    /* do some check */
    if(t!=GLB_T || x!=GLB_X || y!=GLB_Y || z!=GLB_Z) {
      lprintf("WARNING", 0, "Size read from config name (%d,%d,%d,%d) is different from the lattice size!\n",t,x,y,z);
    }
    if(ng!=NG) {
      lprintf("WARNING", 0, "Gauge group read from config name (NG=%d) is not the one used in this code!\n",ng);
    }
    if(strcmp(buf,gf->run_name)!=0) {
      lprintf("WARNING", 0, "Run name [%s] doesn't match conf name [%s]!\n",gf->run_name,buf);
    }

    lprintf("FLOW",0,"Starting from conf [%s]\n",gf->g_start);

    return 0;
  }

  gf->start=1; /* reset gf->start */

  /* try if it match a unit or random start */
  strcpy(buf,gf->g_start);
  slower(buf);
  ret=strcmp(buf,"unit");
  if (ret==0) {
    lprintf("FLOW",0,"Starting a new run from a unit conf!\n");
    return 1;
  }
  ret=strcmp(buf,"random");
  if (ret==0) {
    lprintf("FLOW",0,"Starting a new run from a random conf!\n");
    return 2;
  }
  
  lprintf("ERROR",0,"Invalid starting gauge conf specified [%s]\n",gf->g_start);
  error(1,1,"parse_gstart " __FILE__,"invalid config name");

  return -1;
}

static int parse_lastconf(pg_flow *gf) {

  int ret=0;
  int addtostart=0;

  if (gf->last_conf[0]=='+') { addtostart=1; }
  if(addtostart) {
      ret=sscanf(gf->last_conf,"+%d",&(gf->end));
  } else {
      ret=sscanf(gf->last_conf,"%d",&(gf->end));
  }
  if (ret==1) {
    if (addtostart) gf->end+=gf->start;
    else gf->end++;
    return 0;
  }

  lprintf("ERROR",0,"Invalid last conf specified [%s]\n",gf->last_conf);
  error(1,1,"parse_lastconf " __FILE__,"invalid last config name");

  return -1;
}

int init_mc(pg_flow *gf, char *ifile) {

  int start_t;

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  /* flow defaults */
  strcpy(gf->g_start,"invalid");
  strcpy(gf->run_name,"run_name");
  strcpy(gf->last_conf,"invalid");
  strcpy(gf->conf_dir,"./");
  gf->save_freq=0;
  gf->meas_freq=0;
  gf->therm=0;
  gf->pg_v=&pg_var;

  read_input(pg_var.read,ifile);
  read_input(gf->read,ifile);

  /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.,0.,0.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 1.,
    .SF_BCs = 0
  };
  init_BCs(&BCs_pars);

  /* fix conf_dir: put a / at the end of it */
  start_t=strlen(gf->conf_dir);
  if (gf->conf_dir[start_t-1]!='/') {
    gf->conf_dir[start_t]='/';
    gf->conf_dir[start_t+1]='\0';
  }

  /* set initial configuration */
  start_t=parse_gstart(gf);
  /* set last conf id */
  parse_lastconf(gf);

  
  /* init gauge field */
  switch(start_t) {
    case 0:
      read_gauge_field(add_dirname(gf->conf_dir,gf->g_start));
      gf->therm=0;
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
  
  /* init PG */
  lprintf("MAIN",0,"beta = %2.4f\n",pg_var.beta);
  lprintf("MAIN",0,"Number of thermalization cycles: %d\n",gf->therm);
  lprintf("MAIN",0,"Number of hb-or iterations per cycle: %d\n",pg_var.nit);
  lprintf("MAIN",0,"Number of heatbaths per iteration: %d\n",pg_var.nhb);   
  lprintf("MAIN",0,"Number or overrelaxations per iteration: %d\n",pg_var.nor);

  return 0;

}

int save_conf(pg_flow *gf, int id) {
  char buf[256];
  
  mk_gconf_name(buf, gf, id);
#if NG==2 && !defined(WITH_QUATERNIONS)
  write_gauge_field_su2q(add_dirname(gf->conf_dir,buf));
#else
  write_gauge_field(add_dirname(gf->conf_dir,buf));
#endif
  
  return 0;
}

int end_mc() {
  free_BCs();
  
  /* free memory */
  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  return 0;
}


#undef repr_name

