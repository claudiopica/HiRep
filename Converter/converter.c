/*******************************************************************************
*
* Converter from different formats
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
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
#include "moreio.h"


#ifdef WITH_MPI
#error Please compile without MPI!
#endif

#define true (0==0)
#define false (0==1)

typedef struct _format_type {
  char name[256];
  void (*read)(char*);
  void (*write)(char*);
} format_type;

int nformats=3;
format_type format[3] = {
  { .name="eolexi" , .read=read_gauge_field_eolexi , .write=write_gauge_field_eolexi },
  { .name="henty" ,  .read=read_gauge_field_henty , .write=NULL },
  { .name="mpieo" ,  .read=read_gauge_field , .write=write_gauge_field }
};

typedef struct _filename_type {
  char string[256];
  int size[4];    int size_f;
  int ng;         int ng_f;
  int nf;         int nf_f;
  char repr[256]; int repr_f;
  double beta;    int beta_f;
  double mass;    int mass_f;
  int n;          int n_f;
} filename_type;



filename_type input_filename;
char output_filename[256];
format_type* input_format;
format_type* output_format;
int swapendian=false;


int parse_cnfg_filename(char* filename, filename_type* fn) {
  int hm;
  char *tmp = NULL;
  char *basename;

  basename = filename;
  while ((tmp = strchr(basename, '/')) != NULL) {
    basename = tmp+1;
  }            

  fn->size_f=false;
  fn->ng_f=false;
  fn->nf_f=false;
  fn->repr_f=false;
  fn->beta_f=false;
  fn->mass_f=false;
  fn->n_f=false;

  strcpy(fn->string,basename);
  
  
  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%dr%[A-Z]%*[Nn]f%db%lfm%lfn%d",
      &(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng),fn->repr,&(fn->nf),&(fn->beta),&(fn->mass),&(fn->n));
  if(hm==10) {
    fn->size_f=true;
    fn->ng_f=true;
    fn->nf_f=true;
    fn->repr_f=true;
    fn->beta_f=true;
    fn->mass_f=true;
    fn->n_f=true;
    return hm;
  }


  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng),&(fn->beta),&(fn->n));
  if(hm==7) {
    fn->size_f=true;
    fn->ng_f=true;
    fn->beta_f=true;
    fn->n_f=true;
    return hm;
  }


  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d",
      &(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng));
  if(hm==5) {
    fn->size_f=true;
    fn->ng_f=true;
    return hm;
  }


  return 0;
}


void read_cmdline(int argc, char* argv[]) {
  int i;
  int ai=0, ao=0, aif=0, aof=0, ase=0, av=0;
  char def[256]="mpieo";

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=(i+1<argc)?i+1:0;
    else if (strcmp(argv[i],"-o")==0) ao=(i+1<argc)?i+1:0;
    else if (strcmp(argv[i],"-v")==0) av=(i+1<argc)?i+1:0;
    else if (strcmp(argv[i],"--swapendian")==0) ase=i;
  }
  if(ai+1<argc) if(argv[ai+1][0]!='-') aif=ai+1;
  if(ao+1<argc) if(argv[ao+1][0]!='-') aof=ao+1;

  error(ao==0 || ai==0 || aif==0,1,"parse_cmdline [converter.c]",
      "Syntax: converter -i <input file> <input format> -o <output file> [<output format>] [-v <volume>] [--swapendian]");

  strcpy(output_filename,argv[ao]);
  
  parse_cnfg_filename(argv[ai],&input_filename);

  error(input_filename.size_f==false && av==0,1,"parse_cmdline [converter.c]",
      "Volume required with this input file");
  
  if(input_filename.size_f==false) {
    error(sscanf(argv[av],"%dx%dx%dx%d",&GLB_T,&GLB_X,&GLB_Y,&GLB_Z)!=4,1,"parse_cmdline [converter.c]",
      "Wrong format for volume");
  } else {
    GLB_T=input_filename.size[0];
    GLB_X=input_filename.size[1];
    GLB_Y=input_filename.size[2];
    GLB_Z=input_filename.size[3];
  }

  error(input_filename.ng_f==true && NG!=input_filename.ng,1,"parse_cmdline [converter.c]",
      "Wrong number of colours in file name");
  
  input_format=NULL;
  for(i=0; i<nformats; i++) {
    if(strcmp(argv[aif],format[i].name)==0) {
      input_format = &format[i];
      break;
    }
  }
  error(input_format==NULL,1,"parse_cmdline [converter.c]",
      "Invalid input format (eolexi, henty, mpieo)");

  output_format=NULL;
  if(aof!=0) {
    for(i=0; i<nformats; i++) {
      if(strcmp(argv[aof],format[i].name)==0) {
        output_format = &format[i];
        break;
      }
    }
  } else {
    for(i=0; i<nformats; i++) {
      if(strcmp(def,format[i].name)==0) {
        output_format = &format[i];
        break;
      }
    }
  }
  error(output_format==NULL,1,"parse_cmdline [converter.c]",
      "Invalid output format (eolexi, henty, mpieo)");
  error(output_format->write==NULL,1,"parse_cmdline [converter.c]",
      "Output format not yet available");

  if(ase!=0) swapendian=true;
}

int main(int argc,char *argv[]) {


  /* setup process id and communications */
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,10);

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);

  input_format->read(input_filename.string);
  if(swapendian) gaugefield_swapendian();
  output_format->write(output_filename);

  free_gfield(u_gauge);

  return 0;
}

