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

#include "cinfo.c"


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

#ifdef WITH_QUATERNIONS
#define nformats 10
#else
#define nformats 11
#endif

format_type format[nformats] = {
  { .name="ascii" ,     .read=read_gauge_field_ascii ,     .write=NULL },
  { .name="milc" ,     .read=read_gauge_field_milc ,     .write=NULL },
  { .name="milcn3r" ,     .read=read_gauge_field_milc_no3row ,     .write=NULL },
  { .name="mpieo" ,     .read=read_gauge_field ,           .write=write_gauge_field },
  { .name="eolexi:be" , .read=read_gauge_field_eolexi_BE , .write=write_gauge_field_eolexi_BE },
  { .name="mpieo:be" ,  .read=read_gauge_field_mpieo_BE ,  .write=write_gauge_field_mpieo_BE },
  { .name="eolexi:le" , .read=read_gauge_field_eolexi_LE , .write=write_gauge_field_eolexi_LE },
  { .name="mpieo:le" ,  .read=read_gauge_field_mpieo_LE ,  .write=write_gauge_field_mpieo_LE },
  { .name="fortran" ,  .read=read_gauge_field_fortran ,  .write=NULL },
#ifdef WITH_QUATERNIONS
  { .name="su2q" ,  .read=read_gauge_field_su2q ,  .write=write_gauge_field_su2q },
#endif
  { .name="openQCD" ,  .read=read_gauge_field_openQCD ,  .write=write_gauge_field_openQCD }
};

enum {QUENCHED_CNFG, DYNAMICAL_CNFG, UNKNOWN_CNFG};

typedef struct _filename_type {
  char string[1024];
  char label[256];int label_f;
  int size[4];    int size_f;
  int ng;         int ng_f;
  int nf;         int nf_f;
  char repr[256]; int repr_f;
  double beta;    int beta_f;
  double mass;    int mass_f;
  double kappa;   int kappa_f;
  int n;          int n_f;
  int cnfg_type;
} filename_type;



filename_type input_filename;
char output_filename[1024];
format_type* input_format;
format_type* output_format;
int check=false;


int parse_cnfg_filename(char* filename, filename_type* fn) {
  int hm;
  char *tmp = NULL;
  char *basename;

  basename = filename;
  while ((tmp = strchr(basename, '/')) != NULL) {
    basename = tmp+1;
  }            

  fn->label_f=false;
  fn->size_f=false;
  fn->ng_f=false;
  fn->nf_f=false;
  fn->repr_f=false;
  fn->beta_f=false;
  fn->mass_f=false;
  fn->n_f=false;

  strcpy(fn->string,filename);
  
  
  hm=sscanf(basename,"%[^_]_%dx%dx%dx%d%*[Nn]c%dr%[A-Z]%*[Nn]f%db%lfm%lfn%d",
      fn->label,&(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng),fn->repr,&(fn->nf),&(fn->beta),&(fn->mass),&(fn->n));
  if(hm==11) {
    fn->label_f=true;
    fn->size_f=true;
    fn->ng_f=true;
    fn->nf_f=true;
    fn->repr_f=true;
    fn->beta_f=true;
    fn->mass_f=true;
    fn->n_f=true;
    fn->cnfg_type=DYNAMICAL_CNFG;
    return hm;
  }

  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d",
      &(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng),&(fn->nf),&(fn->beta),&(fn->kappa),&(fn->n));
  if(hm==9) {
    fn->size_f=true;
    fn->ng_f=true;
    fn->nf_f=true;
    fn->beta_f=true;
    fn->kappa_f=true;
    fn->n_f=true;
    fn->cnfg_type=DYNAMICAL_CNFG;
    return hm;
  }


  hm=sscanf(basename,"%[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      fn->label,&(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng),&(fn->beta),&(fn->n));
  if(hm==8) {
    fn->label_f=true;
    fn->size_f=true;
    fn->ng_f=true;
    fn->beta_f=true;
    fn->n_f=true;
    fn->cnfg_type=QUENCHED_CNFG;
    return hm;
  }

  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng),&(fn->beta),&(fn->n));
  if(hm==7) {
    fn->size_f=true;
    fn->ng_f=true;
    fn->beta_f=true;
    fn->n_f=true;
    fn->cnfg_type=QUENCHED_CNFG;
    return hm;
  }


  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d",
      &(fn->size[0]),&(fn->size[1]),&(fn->size[2]),&(fn->size[3]),&(fn->ng));
  if(hm==5) {
    fn->size_f=true;
    fn->ng_f=true;
    fn->cnfg_type=UNKNOWN_CNFG;
    return hm;
  }

  fn->cnfg_type=UNKNOWN_CNFG;
  return 0;
}


void print_cmdline_info() {
  error(1,1,"parse_cmdline [converter.c]",
"\n\
Syntax (1): converter -m\n\
* Show compilation information.\n\n\
Syntax (2): converter -i <input file> [<input format>] -d <output directory> [<output format>] [-l <label>] [-n <n>] [-r <repr>] [-M mass]\n\
* Convert the input file from the input format to the output format. The output file is saved in the output directory, its name is generated automatically. Label, number and representation can be overridden.\n\n\
Syntax (3): converter -i <input file> [<input format>] -o <output file> [<output format>] [-v <volume>]\n\
* Convert the input file from the input format to the output format. The volume must be provided if it cannot be extracted from the input file name.\n\n\
Syntax (4): converter -i <input file> [<input format>] [-v <volume>] --check\n\
* Open the input file assuming the given format and print, if possible, plaquettes and average distance from unitarity for the link variables.\n\n\
Input formats = mpieo (be,default) | mpieo:be | mpieo:le | eolexi:be | eolexi:le | milc | milcn3r | ascii | fortran | su2q\n\
Output formats = mpieo (be,default) | mpieo:be | mpieo:le | eolexi:be | eolexi:le | su2q\n\
");
}


void read_cmdline(int argc, char* argv[]) {
  int i;
  int ai=0, ao=0, ad=0, ade=0, av=0, al=0, an=0, ar=0, aM=0, am=0;
  char *str;
  char def[256]="mpieo";

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=i;
    else if (strcmp(argv[i],"-o")==0) ao=i;
    else if (strcmp(argv[i],"-d")==0) ad=i;
    else if (strcmp(argv[i],"-v")==0) av=i;
    else if (strcmp(argv[i],"-l")==0) al=i;
    else if (strcmp(argv[i],"-r")==0) ar=i;
    else if (strcmp(argv[i],"-n")==0) an=i;
    else if (strcmp(argv[i],"-M")==0) aM=i;
    else if (strcmp(argv[i],"--check")==0) ade=i;
    else if (strcmp(argv[i],"-m")==0) am=i;
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  if(ai==0 || ai+1>=argc) {
    lprintf("ERROR",0,"Input file missing.\n");
    print_cmdline_info();
  }
  parse_cnfg_filename(argv[ai+1],&input_filename);

  if(input_filename.ng_f==true && NG!=input_filename.ng) {
    lprintf("ERROR",0,"You must compile with NG=%d!!!",input_filename.ng);
    print_cmdline_info();
  }
  
  if(input_filename.size_f==false && av==0) {
    lprintf("ERROR",0,"Since I cannot read the volume from the input file name, you must use the -v option.\n");
    print_cmdline_info();
  }
  if(input_filename.size_f==true && av!=0) {
    lprintf("ERROR",0,"You cannot override the volume read from the input file name.\n");
    print_cmdline_info();
  }
  if(input_filename.size_f==false) {
    error(sscanf(argv[av+1],"%dx%dx%dx%d",&GLB_T,&GLB_X,&GLB_Y,&GLB_Z)!=4,1,"parse_cmdline [converter.c]",
      "Wrong format for volume");
  } else {
    GLB_T=input_filename.size[0];
    GLB_X=input_filename.size[1];
    GLB_Y=input_filename.size[2];
    GLB_Z=input_filename.size[3];
  }

  input_format=NULL;
  str=def;
  if(ai+2<argc)
  if(argv[ai+2][0]!='-')
    str=argv[ai+2];
  for(i=0; i<nformats; i++) {
    if(strcmp(str,format[i].name)==0) {
      input_format = &format[i];
      break;
    }
  }
  if(input_format==NULL) {
    lprintf("ERROR",0,"Invalid input format.\n");
    print_cmdline_info();
  }

  if(ao!=0 && ad==0 && ade==0) {

    if(ao+1>=argc) {
      lprintf("ERROR",0,"Output file missing.\n");
      print_cmdline_info();
    }
    strcpy(output_filename,argv[ao+1]);
  
    output_format=NULL;
    str=def;
    if(ao+2<argc)
    if(argv[ao+2][0]!='-')
      str=argv[ao+2];
    for(i=0; i<nformats; i++) {
      if(strcmp(str,format[i].name)==0) {
        output_format = &format[i];
        break;
      }
    }
    if(output_format==NULL) {
      lprintf("ERROR",0,"Invalid output format.\n");
      print_cmdline_info();
    }
    
    if(ar!=0) lprintf("WARNING",0,"-r option ignored.\n");
    if(al!=0) lprintf("WARNING",0,"-l option ignored.\n");
    if(an!=0) lprintf("WARNING",0,"-n option ignored.\n");
    if(aM!=0) lprintf("WARNING",0,"-M option ignored.\n");

  } else if(ao==0 && ad!=0 && ade==0) {

    if(input_filename.cnfg_type==UNKNOWN_CNFG) {
      lprintf("ERROR",0,"Since the input file name is not in one of the standard formats, you cannot use the -d option.\n");
      print_cmdline_info();
    }

    if(input_filename.repr_f==true && ar!=0) {
      lprintf("ERROR",0,"You cannot override the representation read from the input file name.\n");
      print_cmdline_info();
    }

    char tmp[1024];

    if(ad+1>=argc) {
      lprintf("ERROR",0,"Output directory missing.\n");
      print_cmdline_info();
    }
    sprintf(tmp,"%s/",argv[ad+1]);
    
    if(al!=0) {
      if(al+1>=argc) {
        lprintf("ERROR",0,"Label missing.\n");
        print_cmdline_info();
      }
      sprintf(output_filename,"%s%s_",tmp,argv[al+1]);
    } else {
      if(input_filename.label_f==false) {
        lprintf("ERROR",0,"Since I cannot read the label from the input file name, you must use the -l option.\n");
        print_cmdline_info();
      }
      sprintf(output_filename,"%s%s_",tmp,input_filename.label);
    }
    
    sprintf(tmp,"%s%dx%dx%dx%dnc%d",output_filename,GLB_T,GLB_X,GLB_Y,GLB_Z,NG);
    
    if(input_filename.cnfg_type==QUENCHED_CNFG) {

      if(input_filename.beta_f==false) {
        lprintf("ERROR",0,"beta missing in input file. This is strange!!!\n");
        print_cmdline_info();
      }
      sprintf(output_filename,"%sb%.6f",tmp,input_filename.beta);
      strcpy(tmp,output_filename);

      if(aM!=0) lprintf("WARNING",0,"-M option ignored.\n");
      
    } else if(input_filename.cnfg_type==DYNAMICAL_CNFG) {
    
      if(ar!=0) {
        if(ar+1>=argc) {
          lprintf("ERROR",0,"Representation missing.\n");
          print_cmdline_info();
        }
        sprintf(output_filename,"%sr%s",tmp,argv[ar+1]);
      } else {
        if(input_filename.repr_f==false) {
          lprintf("ERROR",0,"Since I cannot read the representation from the input file name, you must use the -r option.\n");
          print_cmdline_info();
        }
        sprintf(output_filename,"%sr%s",tmp,input_filename.repr);
      }
      
      if(input_filename.nf_f==false) {
        lprintf("ERROR",0,"nf missing in input file. This is strange!!!\n");
        print_cmdline_info();
      }
      sprintf(tmp,"%snf%d",output_filename,input_filename.nf);

      if(input_filename.beta_f==false) {
        lprintf("ERROR",0,"beta missing in input file. This is strange!!!\n");
        print_cmdline_info();
      }
      sprintf(output_filename,"%sb%.6f",tmp,input_filename.beta);

      if(aM!=0) {
        if(aM+1>=argc) {
          lprintf("ERROR",0,"Mass missing.\n");
          print_cmdline_info();
        }
        sprintf(tmp,"%sm%s",output_filename,argv[aM+1]);
      } else {
        if(input_filename.mass_f==false && input_filename.kappa_f==false) {
          lprintf("ERROR",0,"mass and kappa missing in input file. This is strange!!!\n");
          print_cmdline_info();
        }
        if(input_filename.mass_f==true)
          sprintf(tmp,"%sm%.6f",output_filename,input_filename.mass);
        else if(input_filename.kappa_f==true)
          sprintf(tmp,"%sm%.6f",output_filename,-.5/input_filename.kappa+4.);
      }
      
    }
    
    if(an!=0) {
      if(an+1>=argc) {
        lprintf("ERROR",0,"Number missing.\n");
        print_cmdline_info();
      }
      sprintf(output_filename,"%sn%d",tmp,atoi(argv[an+1]));
    } else {
      if(input_filename.n_f==false) {
        lprintf("ERROR",0,"Since I cannot read the number from the input file name, you must use the -n option.\n");
        print_cmdline_info();
      }
      sprintf(output_filename,"%sn%d",tmp,input_filename.n);
    }

    output_format=NULL;
    str=def;
    if(ad+2<argc)
    if(argv[ad+2][0]!='-')
      str=argv[ad+2];
    for(i=0; i<nformats; i++) {
      if(strcmp(str,format[i].name)==0) {
        output_format = &format[i];
        break;
      }
    }
    if(output_format==NULL) {
      lprintf("ERROR",0,"Invalid output format.\n");
      print_cmdline_info();
    }
  
  } else if(ao==0 && ad==0 && ade!=0) {
    
    if(ar!=0) lprintf("WARNING",0,"-r option ignored.\n");
    if(al!=0) lprintf("WARNING",0,"-l option ignored.\n");
    if(an!=0) lprintf("WARNING",0,"-n option ignored.\n");
    if(aM!=0) lprintf("WARNING",0,"-M option ignored.\n");
    
    check=true;
  
  } else {

    lprintf("ERROR",0,"One option among -o -d --check must be used.\n");
    print_cmdline_info();

  }

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
  
  if(check) {
    full_plaquette();
    
    int mu;
    suNg A;
    double norm=0., tmp;
    
    _MASTER_FOR(&glattice,ix) {
      for(mu=0; mu<4; mu++) {
        _suNg_times_suNg_dagger(A,*pu_gauge(ix,mu),*pu_gauge(ix,mu));
        _suNg_sqnorm_m1(tmp,A);
        norm+=tmp;
      }
    }
    norm=sqrt(norm/(8.*NG*NG)/GLB_VOLUME);
    
    lprintf("MAIN",0,"Average distance from unitarity = %e\n",norm);
  } else {
    output_format->write(output_filename);
  }

  free_gfield(u_gauge);

  return 0;
}

