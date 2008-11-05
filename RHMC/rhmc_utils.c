/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "logger.h"
#include "rhmc_utils.h"
#include "random.h"
#include "io.h"
#include "representation.h"
#include "update.h"
#include "memory.h"

#include <stdio.h>
#include <string.h>


input_rhmc rhmc_var = init_input_rhmc(rhmc_var);

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

static void mk_gconf_name(char *name, rhmc_flow *rf, int id) {
  sprintf(name,"%s_%dx%dx%dx%dnc%dr%snf%db%.6fm%.6fn%d",
           rf->run_name,GLB_T,GLB_X,GLB_Y,GLB_Z,NG,repr_name,
           rf->rhmc_v->rhmc_p.nf,rf->rhmc_v->rhmc_p.beta,-rf->rhmc_v->rhmc_p.mass,
           id);
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

static int parse_gstart(rhmc_flow *rf) {

  int t, x, y, z, ng, nf;
  double beta, mass;
  int ret=0;
  char buf[256];

  ret=sscanf(rf->g_start,"%[^_]_%dx%dx%dx%dnc%dr" repr_name "nf%db%lfm%lfn%d",
      buf,&t,&x,&y,&z,&ng,&nf,&beta,&mass,&rf->start);

  if(ret==10) { /* we have a correct file name */
    /* increase rf->start: this will be the first conf id */
    rf->start++;

    /* do some check */
    if(t!=GLB_T || x!=GLB_X || y!=GLB_Y || z!=GLB_X) {
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

  /* try if it match a unit or random start */
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

static int parse_lastconf(rhmc_flow *rf) {

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

int init_mc(rhmc_flow *rf, char *ifile) {

  int start_t;

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
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
  rf->rhmc_v=&rhmc_var;

  read_input(rhmc_var.read,ifile);
  read_input(rf->read,ifile);

  /* fix conf_dir: put a / at the end of it */
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
  represent_gauge_field();
  

  /* init RHMC */
  
  rhmc_var.rhmc_p.integrator=&O2MN_multistep;
  rhmc_var.rhmc_p.mshift_solver=&cg_mshift; 
  rhmc_var.rhmc_p.MD_par=&rhmc_var.int_p;
  init_rhmc(&rhmc_var.rhmc_p);

  return 0;

}

int save_conf(rhmc_flow *rf, int id) {
  char buf[256];
  
  mk_gconf_name(buf, rf, id);
  write_gauge_field(add_dirname(rf->conf_dir,buf));

  return 0;
}

int end_mc() {
  free_rhmc();

  /* free memory */
  free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
  free_gfield_f(u_gauge_f);
#endif

  return 0;
}


#undef repr_name
