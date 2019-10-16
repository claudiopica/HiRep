// This code will contain all the contractions necessary for rho to pi pi scattering
#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "scattering.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "clover_tools.h"
#include "setup.h"

#include "cinfo.c"
#include "IOroutines.c"
#include "scatter_functions.h"

#define PI M_PI

int main(int argc,char *argv[])
{
  //int px,py,pz, px2, py2, pz2, px3, py3, pz3;

  FILE* list;
  int tau=0;
  double m[256];

  //Copy I/O from another file
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  setup(&list, m);
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  int numsources = mes_var.nhits;
  char *path=mes_var.outdir;
  int Nmom;
  lprintf("MAIN",0,"Boundary conditions: %s\n",mes_var.bc);
  lprintf("MAIN",0,"The momenta are: %s\n",mes_var.p);
  int **p = getmomlist(mes_var.p,&Nmom);
  lprintf("MAIN",0,"Number of momenta: %d\n",Nmom);
  lprintf("MAIN",0,"The momenta are:\n");
  for(int i=0; i<Nmom; i++){
    lprintf("MAIN",0,"p%d = (%d, %d, %d)\n", i+1, p[i][0], p[i][1], p[i][2]);
  }

while(1){
    struct timeval start, end, etime;
    gettimeofday(&start,0);
    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    struct mo_0 *mo_p0[numsources]; 
    struct mo_p *mo_p[Nmom][numsources];
    for(int i=0; i<numsources; i++){
        mo_p0[i] = (struct mo_0*) malloc(sizeof(struct mo_0));
        for(int j=0; j<Nmom; j++){
            mo_p[j][i] = (struct mo_p*) malloc(sizeof(struct mo_p));
        }
        lprintf("MAIN",0,"Initiating mo, source = %d\n",i);
        init_mo_0(mo_p0[i]);
        for (int j=0;j<Nmom;j++) init_mo_p(mo_p[j][i],p[j][0],p[j][1],p[j][2]);
    }

    for (int src=0;src<numsources;++src)
    {
	    struct src_common src0;
	    struct src_p *src_pn = (struct src_p*) malloc(Nmom*sizeof(struct src_p));
	    struct prop_common prop0;
	    struct prop_p *p_p = (struct prop_p*) malloc(Nmom*sizeof(struct prop_p));


	    init_src_common(&src0,tau);
	    make_prop_common(&prop0, &src0, 4, tau,mes_var.bc);
	    gen_mo_0(mo_p0[src], &prop0, &src0, tau);

        for(int i=0; i<Nmom; i++){
            init_src_p(src_pn + i, &src0, p[i][0], p[i][1], p[i][2]);
            make_prop_p(p_p + i, src_pn + i, &src0, 4, tau, mes_var.bc);
            gen_mo_p(mo_p[i][src], &prop0, p_p + i, &src0, tau);
        }

	    free_src_common(&src0);
	    free_prop_common(&prop0);
        for(int i=0; i<Nmom; i++){
            free_src_p(src_pn + i);
            free_prop_p(p_p + i);
        }
    }
    lprintf("MAIN",0,"num sources: %d, path: %s\n",numsources,path);
    //IOold_0(mo_p0, numsources, path);
    IO_json_0(mo_p0, numsources, path);
    for(int i=0; i<Nmom; i++){
        //IOold_p(mo_p[i], numsources, path);
        IO_json_p(mo_p[i], numsources, path);
    }

    for(int src=0;src<numsources;src++){
	    free_mo_0(mo_p0[src]);
        for(int i=0; i<Nmom;i++){
            free_mo_p(mo_p[i][src]);
        }
    }

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration : analysed in [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);

    if(list==NULL) break;
}
  lprintf("DEBUG",0,"ALL done, deallocating\n");

  if(list!=NULL) fclose(list);
  finalize_process();
  free_BCs();
  free_gfield(u_gauge);
  free_propagator_eo();

  return 0;
}
