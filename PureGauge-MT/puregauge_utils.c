/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "logger.h"
#include "puregauge_utils.h"
#include "random.h"
#include "io.h"
#include "update.h"
#include "memory.h"
#include "utils.h"

#include <stdio.h>
#include <string.h>


input_puregauge puregauge_var = init_input_puregauge(puregauge_var);

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

/* build configuration name */
static void mk_gconf_name(char *name, puregauge_flow *rf, int id) {
  sprintf(name,"%s_%dx%dx%dx%dnc%db%.6fn%d",
           rf->run_name,GLB_T,GLB_X,GLB_Y,GLB_Z,NG,
           rf->puregauge_v->puregauge_p.beta,
           id);
}

/* add dirname to filename and return it */
static char *add_dirname(char *dirname, char *filename) {
  static char buf[256];
  strcpy(buf,dirname);
  return strcat(buf,filename);
}

#include <ctype.h>

/* convert string to lowercase */
static void slower(char *str) {
  while (*str) {
    *str=(char)(tolower(*str));
    ++str;
  }
}

/* read g_start string and decide what the initial config should be.
 * It also fill the content of the variable start in puregauge_flow with the
 * id of the first config to be used.
 * returns:
 * 0 => g_start is a valid file name which should be used as starting config
 * 1 => use unit gauge config
 * 2 => use a random config
 */
static int parse_gstart(puregauge_flow *rf) {

  int t, x, y, z, ng;
  double beta;
  int ret=0;
  char buf[256];

  ret=sscanf(rf->g_start,"%[^_]_%dx%dx%dx%dnc%db%lfn%d",
      buf,&t,&x,&y,&z,&ng,&beta,&rf->start);

  if(ret==8) { /* we have a correct file name */
    /* increase rf->start: this will be the first conf id */
    rf->start++;

    /* do some check */
    if(t!=GLB_T || x!=GLB_X || y!=GLB_Y || z!=GLB_Z) {
      lprintf("WARNING", 0, "Size read from config name (%d,%d,%d,%d) is different from the lattice size!\n",t,x,y,z);
    }
    if(ng!=NG) {
      lprintf("WARNING", 0, "Gauge group read from config name (NG=%d) is not the one used in this code!\n",ng);
    }
    /* dare un warning ?
       if(nf!=) {
       lprintf("TESTCP",0,"Invalid gauge group in file name!\n");
       error(1,1,"parse_gstart " __FILE__,"invalid gauge group");
       }
       */
    if(strcmp(buf,rf->run_name)!=0) {
      lprintf("WARNING", 0, "Run name [%s] doesn't match conf name [%s]!\n",rf->run_name,buf);
    }

    lprintf("FLOW",0,"Starting from conf [%s]\n",rf->g_start);

    return 0;
  }

  rf->start=1; /* reset rf->start */

  /* try other matches */
  strcpy(buf,rf->g_start);
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
  
  lprintf("ERROR",0,"Invalid starting gauge conf specified [%s]\n",rf->g_start);
  error(1,1,"parse_gstart " __FILE__,"invalid config name");

  return -1;
}

/* read last_conf string and fill the end parameter
 * in puregauge_flow with the id of the last conf to be generated.
 * last_conf must be one of:
 * "%d" => this will be the id of the last conf generated
 * "+%d" => if a '+' is prepended then this run will be %d config long,
 *          i.e. the end will be start+%d 
 *  rerunts:
 *  0 => string was valid
 *  -1 => invalid format
 */
static int parse_lastconf(puregauge_flow *rf) {

  int ret=0;
  int addtostart=0;

  if (rf->last_conf[0]=='+') { addtostart=1; }
  if(addtostart) {
      ret=sscanf(rf->last_conf,"+%d",&(rf->end));
  } else {
      ret=sscanf(rf->last_conf,"%d",&(rf->end));
  }
  if (ret==1) {
    if (addtostart) rf->end+=rf->start;
    else rf->end++;
    return 0;
  }

  lprintf("ERROR",0,"Invalid last conf specified [%s]\n",rf->last_conf);
  error(1,1,"parse_lastconf " __FILE__,"invalid last config name");

  return -1;
}

/* Initialize the Monte Carlo.
 * This performs the following operations:
 * 1) read from the specified input file the flow variables 
 *    and the puregauge parameters;
 * 2) set the starting gauge field
 * 3) init the puregauge update
 */
int init_mc(puregauge_flow *rf, char *ifile) {

  int start_t;

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  /* flow defaults */
  strcpy(rf->g_start,"invalid");
  rf->r_start[0]='\0';
  strcpy(rf->run_name,"run_name");
  strcpy(rf->last_conf,"invalid");
  strcpy(rf->conf_dir,"./");
  rf->save_freq=0;
  rf->meas_freq=0;
  rf->puregauge_v=&puregauge_var;

  read_input(puregauge_var.read,ifile);
  read_input(rf->read,ifile);

  /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.,0.,0.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 1.,
    .SF_BCs = 0
  };
  init_BCs(&BCs_pars);

  /* fix conf_dir name: put a / at the end of string */
  start_t=strlen(rf->conf_dir);
  if (rf->conf_dir[start_t-1]!='/') {
    rf->conf_dir[start_t]='/';
    rf->conf_dir[start_t+1]='\0';
  }

  /* set initial configuration */
  start_t=parse_gstart(rf);
  /* set last conf id */
  parse_lastconf(rf);

  /* init gauge field */
  switch(start_t) {
    case 0:
      read_gauge_field(add_dirname(rf->conf_dir,rf->g_start)); 
      break;
    case 1:
      unit_u(u_gauge);
      break;
    case 2:
      random_u(u_gauge);
      break;
  }
  apply_BCs_on_fundamental_gauge_field();
  complete_gf_sendrecv(u_gauge);
  

  /* init GAUGE */
  init_puregauge(&puregauge_var.puregauge_p);

  return 0;

}

/* save the gauge config with the specified id */
int save_conf(puregauge_flow *rf, int id) {
  char buf[256];
  
  mk_gconf_name(buf, rf, id);
  write_gauge_field(add_dirname(rf->conf_dir,buf));

  return 0;
}

/* clean up memory */
int end_mc() {
  free_puregauge();
  free_BCs();

  /* free memory */
  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  return 0;
}

#undef repr_name

