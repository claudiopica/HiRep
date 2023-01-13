/*******************************************************************************
 *
 * Computation of Renomalization constants (Z_a,Z_q,Z_s,Z_ps,Z_t,Z_m,Z_v)  
 * factors with gauge fixed momentum sources. 
 *
 * Rudy Arthur
 *
 * Modified by Vincent Drach and Ari Hietanen
 *
 * NOCOMPILE= ROTATED_SF || BASIC_SF || WITH_QUATERNIONS
 * NOCOMPILE= REPR_ADJOINT
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 * 
 *******************************************************************************/

#include "libhr.h"
#include <string.h>

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

#define PI 3.141592653589793238462643383279502884197
static void twist_XYZ_bc(double theta_x, double theta_y, double theta_z) {

  int index;
  int ix,iy,iz,it;
  suNf *u;
  suNf utmp;
  hr_complex eith_x, eith_y, eith_z;
  eith_x = cexp(I*PI*theta_x/(double)GLB_X);
  eith_y = cexp(I*PI*theta_y/(double)GLB_Y);
  eith_z = cexp(I*PI*theta_z/(double)GLB_Z);
  //  eith_x.re = cos(PI*theta_x/(double)GLB_X); eith_x.im = sin(PI*theta_x/(double)GLB_X);
  //eith_y.re = cos(PI*theta_y/(double)GLB_Y); eith_y.im = sin(PI*theta_y/(double)GLB_Y);
  //eith_z.re = cos(PI*theta_z/(double)GLB_Z); eith_z.im = sin(PI*theta_z/(double)GLB_Z);

  for (it=0;it<T_EXT;++it) for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
          index=ipt_ext(it,ix,iy,iz);
          u=pu_gauge_f(index,1);
          _suNf_mulc(utmp, eith_x, *u);
          *u=utmp;
          u=pu_gauge_f(index,2);
          _suNf_mulc(utmp, eith_y, *u);
          *u=utmp;
          u=pu_gauge_f(index,3);
          _suNf_mulc(utmp, eith_z, *u);
          *u=utmp;      
        }
    
}

/* Renormalization parameters */
typedef struct input_renormalization {
  char mstring[256];
  char configlist[256]; /* list of configuration */
  double precision;
  int ne;
  int n_mom;
  int n_twist;
  int pt_out;
  int px_out;
  int py_out;
  int pz_out;
  int pt_in;
  int px_in;
  int py_in;
  int pz_in;

  /* for the reading function */
  input_record_t read[15];
} input_renormalization;

#define init_input_renormalization(varname)                             \
  {                                                                     \
    .read={                                                             \
      {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring}, \
      {"Configuration list:", "mes:configlist = %s", STRING_T, &(varname).configlist},\
      {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision}, \
      {"non-excepional configuration or not?", "mes:ne = %d", INT_T, &(varname).ne}, \
      {"number of momenta in each direction", "mes:n_mom = %d", INT_T, &(varname).n_mom}, \
      {"number of twists at each momentum", "mes:n_twist = %d", INT_T, &(varname).n_twist}, \
      {"mom1 t component", "mes:pt_out = %d", INT_T, &(varname).pt_out}, \
      {"mom1 t component", "mes:px_out = %d", INT_T, &(varname).px_out}, \
      {"mom1 t component", "mes:py_out = %d", INT_T, &(varname).py_out}, \
      {"mom1 t component", "mes:pz_out = %d", INT_T, &(varname).pz_out}, \
      {"mom2 t component", "mes:pt_in = %d", INT_T, &(varname).pt_in},  \
      {"mom2 t component", "mes:px_in = %d", INT_T, &(varname).px_in},  \
      {"mom2 t component", "mes:py_in = %d", INT_T, &(varname).py_in},  \
      {"mom2 t component", "mes:pz_in = %d", INT_T, &(varname).pz_in},  \
      {NULL, NULL, INT_T, NULL}                                         \
    }                                                                   \
  }


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "renormalization.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_renormalization mes_var = init_input_renormalization(mes_var);


typedef struct {
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;


int parse_cnfg_filename(char* filename, filename_t* fn) {
  int hm;
  char *tmp = NULL;
  char *basename;

  basename = filename;
  while ((tmp = strchr(basename, '/')) != NULL) {
    basename = tmp+1;
  }            

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif
  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%dr" repr_name "%*[Nn]f%db%lfm%lfn%d",
            &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&(fn->m),&(fn->n));
  if(hm==9) {
    fn->m=-fn->m; /* invert sign of mass */
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }
#undef repr_name

  double kappa;
  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d",
            &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&kappa,&(fn->n));
  if(hm==9) {
    fn->m = .5/kappa-4.;
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }

  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
            &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
            &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  fn->type=UNKNOWN_CNFG;
  return UNKNOWN_CNFG;
}




int main(int argc,char *argv[]) {
  int i,k;
  FILE* list;
  filename_t fpars;
  int nm;
  double m[256];
  
  spinor_field* source;
  spinor_field* prop_in;
  spinor_field* prop_out;

  /* setup process communications */
  setup_process(&argc,&argv);

  setup_gauge_fields();

  read_input(mes_var.read,get_input_filename());
  strcpy(list_filename,mes_var.configlist);
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (strcmp(list_filename,"")!=0) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 

#ifdef GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#endif
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  nm=1;

  lprintf("MAIN",0,"Inverter precision = %e\n",mes_var.precision);
  m[0] = -atof(mes_var.mstring);
  for(k=0;k<nm;k++)
    lprintf("MAIN",0,"Mass[%d] = %f\n",k,m[k]);

  lprintf("MAIN",0,"Number of maximum monentum component\n",mes_var.n_mom);
  lprintf("MAIN",0,"Number of twists per momentum\n",mes_var.n_twist);
  if(mes_var.ne) lprintf("MAIN",0,"Doing non-exceptional configuration\n");
  lprintf("MAIN",0,"Momentum 1 (%d, %d, %d, %d)\n",mes_var.pt_out,mes_var.px_out,mes_var.py_out,mes_var.pz_out);
  if(mes_var.ne) lprintf("MAIN",0,"Momentum 2 (%d, %d, %d, %d)\n",mes_var.pt_in,mes_var.px_in,mes_var.py_in,mes_var.pz_in);

  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [measure_Z_mom.c]" ,
          "Failed to open list file\n");
  }

  i=0;
  
  source = alloc_spinor_field_f(4,&glattice);
  prop_in = alloc_spinor_field_f(4*nm*NF,&glattice);
  prop_out = alloc_spinor_field_f(4*nm*NF,&glattice);

  suNf_field* u_gauge_old_f=alloc_gfield_f(&glattice);
  suNg_field* u_gauge_old=alloc_gfield(&glattice);

  while(1) {
    struct timeval start, end, etime;

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

  
    parse_cnfg_filename(cnfg_filename,&fpars);

    GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
    error(fpars.type==UNKNOWN_CNFG,1,"measure_Z_mom.c","Bad name for a configuration file");
    error(fpars.nc!=NG,1,"measure_Z_mom.c","Bad NG");
    
    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
   
    suNg_field_copy(u_gauge_old,u_gauge);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());
    full_plaquette();

    //Fix Gauge
    double p2 = calc_plaq(u_gauge);
    lprintf("MAIN",0,"initial plaq %1.6f\n",p2);

    gettimeofday(&start,0);
    double act = gaugefix(10, //= 0, 1, 2, 3 for Coulomb guage else Landau
                          1.8,	//overrelax
                          10000,	//maxit
                          1e-10, //tolerance
                          u_gauge //gauge
                          );
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration #%d: Gauge Fixed in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    lprintf("MAIN",0,"action  %1.6f\n",act);
    p2 = calc_plaq(u_gauge);
    lprintf("MAIN",0,"fixed gauge plaq %1.6f\n",p2);

    represent_gauge_field();
    gettimeofday(&start,0);
    
    suNf_field_copy(u_gauge_old_f,u_gauge_f);

    init_propagator_eo(nm, m, mes_var.precision);

    int l,j,tw,num;
    num = 0;
    double mom_in[4], mom_out[4];
    mom_out[0]=mes_var.pt_out ; mom_out[1]=mes_var.px_out ; mom_out[2]=mes_var.py_out ; mom_out[3]=mes_var.pz_out ;
    mom_in[0] =mes_var.pt_in ;  mom_in[1] =mes_var.px_in ;  mom_in[2] =mes_var.py_in ;  mom_in[3] =mes_var.pz_in ;

    for(l=1; l<=mes_var.n_mom; ++l){
      double p_in[4], p_out[4];
      for(tw=-mes_var.n_twist;tw<mes_var.n_twist+1;tw++){
        suNf_field_copy(u_gauge_f,u_gauge_old_f);
        double twist = (double)tw*0.5;//(double)tw/(double)(mes_var.n_twist+1);
        lprintf("TEST",0,"<p> before twist %1.6f\n",avr_plaquette());
        twist_XYZ_bc(twist * mes_var.px_in, twist * mes_var.py_in, twist * mes_var.pz_in);
        lprintf("TEST",0,"<p> after twist %1.6f\n",avr_plaquette());

        p_in[0] = mom_in[0]*l; 
        p_in[1] = mom_in[1]*l; 
        p_in[2] = mom_in[2]*l; 
        p_in[3] = mom_in[3]*l; 


        for (k=0;k<NF;++k){
          create_gauge_fixed_momentum_source(source,p_in[0],p_in[1],p_in[2],p_in[3],k);
          calc_propagator(prop_in + 4*k,source,4);//4 for spin components
        } 

        if(mes_var.ne){
          suNf_field_copy(u_gauge_f,u_gauge_old_f);
          twist_XYZ_bc(twist * mes_var.px_out, twist * mes_var.py_out, twist * mes_var.pz_out);

          p_out[0] = mom_out[0]*l; 
          p_out[1] = mom_out[1]*l;
          p_out[2] = mom_out[2]*l; 
          p_out[3] = mom_out[3]*l;
          for (k=0;k<NF;++k){
            create_gauge_fixed_momentum_source(source,p_out[0],p_out[1],p_out[2],p_out[3],k);
            calc_propagator(prop_out + 4*k,source,4);//4 for spin components
          } 
        } else {
          p_out[0] = p_in[0]; 
          p_out[1] = p_in[1];
          p_out[2] = p_in[2]; 
          p_out[3] = p_in[3];
          for(j=0;j<4*NF;j++) spinor_field_copy_f(&prop_out[j],&prop_in[j]);
        }
        lprintf("LOOK",10,"%g%g%g%g %g%g%g%g twist %g",p_in[0],p_in[1],p_in[2],p_in[3], p_out[0],p_out[1],p_out[2],p_out[3], twist);
        measure_renormalization(prop_in, prop_out, nm, p_in[0],p_in[1],p_in[2],p_in[3], p_out[0],p_out[1],p_out[2],p_out[3]);
        char label[256];
        sprintf(label,"NPR mom_idx %d twist %d ",num, tw);
        print_renormalization(i, nm, m, label, p_in[0],p_in[1],p_in[2],p_in[3], p_out[0],p_out[1],p_out[2],p_out[3]);
        num++;
      }
    }
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration #%d: analysed in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    if(list==NULL) break;
  }
  free_spinor_field_f(source);
  free_spinor_field_f(prop_in);
  free_spinor_field_f(prop_out);

  if(list!=NULL) fclose(list);

  free_propagator_eo(); 

  free_BCs();

  free_gfield(u_gauge);
  free_gfield(u_gauge_old);
  free_gfield_f(u_gauge_old_f);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

