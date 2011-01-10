#define MAIN_PROGRAM

#include <qdp-config.h>
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

#define true (0==0)
#define false (0==1)

#ifdef WITH_MPI
#error Please compile without MPI!
#endif

int lattice_size[NDIM] = { 16, 16, 16, 32 };
char output_filename[256];
char input_filename[256];

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

filename_type input_filetype;


void read_cmdline(int argc, char* argv[]) {
  int i;
  int ai=0, ao=0;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=i+1;
    else if (strcmp(argv[i],"-o")==0) ao=i+1;
  }

  error(ao==0 || ai==0 ,1,"parse_cmdline [chroma_converter.c]",
      "Syntax: chroma_converter -i <input file> -o <output file>");

  strcpy(output_filename,argv[ao]);
  strcpy(input_filename,argv[ai]);
  
 
  input_filetype.size[0]=lattice_size[3];
  input_filetype.size[1]=lattice_size[0];
  input_filetype.size[2]=lattice_size[1];
  input_filetype.size[3]=lattice_size[2];
  input_filetype.ng=QDP_Nc;
  input_filetype.size_f=true;
  input_filetype.ng_f=true;
  
  
  GLB_T=input_filetype.size[0];
  GLB_X=input_filetype.size[1];
  GLB_Y=input_filetype.size[2];
  GLB_Z=input_filetype.size[3];

  if(params==NULL) {
    params = malloc(sizeof(params_struct));
  }

  params->T=lattice_size[3];
  params->X=lattice_size[0];
  params->Y=lattice_size[1];
  params->Z=lattice_size[2];


}


int
main(int argc, char *argv[])
{
  int status = 1;
  int mu;
  QDP_ColorMatrix *U[NDIM];
  QLA_Real plaq;

  /* setup process id and communications */
  read_cmdline(argc, argv);

  /* logger setup */
  logger_setlevel(0,1000);

  /* setup communication geometry HiRep*/
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry HiRep*/
  geometry_mpi_eo();
  /* test_geometry_mpi_eo();*/

  /* alloc global gauge fields HiRep*/
  u_gauge=alloc_gfield(&glattice);

  /* start QDP */
  QDP_initialize(&argc, &argv);

  /* set lattice size and create layout */
  QDP_set_latsize(NDIM, lattice_size);
  QDP_create_layout();

  /* allocate the gauge field */
  for (mu = 0; mu < NDIM; mu++) {
    U[mu] = QDP_create_M();
  }

  /* read gauge field */
  if (read_gauge(U,input_filename) != 0) {
    printf0("ERROR: read_gauge(%s)\n", argv[1]);
    goto end;
  }

  input_filetype.beta=params->beta;
  input_filetype.beta_f=true;

  /* translate the gauge field in HiRep notation */
  
  chroma_to_HiRep(U);

  /* Compute plaquette with QDP */
  plaq = plaquette(U);

 /* Display the value */
  printf0("Chroma plaquette{%s} = %g\n\n", argv[1],
           plaq / (QDP_volume() * QDP_Nc * NDIM * (NDIM - 1) / 2 ));


  /* Compute plaquette with HiRep */
  lprintf("Converter",0,"HiRep plaq = %8.6f \n", avr_plaquette()); 
  full_plaquette();


  write_gauge_field(output_filename);

  /* delete the gauge field */
  for (mu = 0; mu < NDIM; mu++) QDP_destroy_M(U[mu]);  

   status = 0;
end:
  /* shutdown QDP */
  QDP_finalize();

  free_gfield(u_gauge);

  return status;
}
