
// Test code for meson scattering
// Obectives:
// 1. Compare the new scattering code (measure_scattering.c) with the old code (scattering_tools.c) at zero momentum. The results should agree up to numerical precision.
// 2. Check the implementation of momenta in the new code by trying two different total momenta at the sink: one which matches the total input at the source and one that doesn't. The former should give a non-zero result while the latter should be consistent with 0 up to numerical precision.

//For some reason the linker requires the following line
#define MAIN_PROGRAM

#include "scattering_tools.c"
#include "IOroutines.c"
#include "observables.h"
#include "communications.h"

#define INDEX(px,py,pz,n_mom,tc) ((px + n_mom)*(2*n_mom+1)*(2*n_mom+1)*(GLB_T)+(py + n_mom)*(2*n_mom+1)*(GLB_T)+(pz + n_mom)*(GLB_T)+ (tc))
// There don't seem to be any noise sources with momentum, so I'll just add it a posteriori using the following function (I hope it's correct!)
// If this works, perhaps it should be moved to the sources.c file?
/*void addMomentum(spinor_field* out, spinor_field* in, int px, int py, int pz)
{
  int c[4];
  int beta, color;
  double pdotx;

  for (beta=0;beta<4;++beta){
    spinor_field_zero_f(&out[beta]);
  }
  lprintf("Adding momentum to the source",0,"mom = (%d,%d,%d)",px,py,pz);

  for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++) for(c[3]=0; c[3]<Z; c[3]++) {
	  pdotx = 2.*PI*((double)(c[1]+zerocoord[1])*(double)px/(double)GLB_X +
                         (double)(c[2]+zerocoord[2])*(double)py/(double)GLB_Y +
                         (double)(c[3]+zerocoord[3])*(double)pz/(double)GLB_Z );
	  for (beta=0;beta<4;++beta) for (color=0; color<NF; ++color){
	     _FIELD_AT(&out[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re = _FIELD_AT(&in[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re * cos(pdotx) - _FIELD_AT(&in[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].im * sin(pdotx);
	     _FIELD_AT(&out[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].im = _FIELD_AT(&in[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re * sin(pdotx) + _FIELD_AT(&in[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].im * cos(pdotx) ;
	  }
  }

  //Not sure what this loop does, but it's in every source definition, so it must be important?
  for (beta=0;beta<4;++beta){
     start_sf_sendrecv(out + beta);
     complete_sf_sendrecv(out + beta);
  }
}*/
static void do_global_sum(meson_observable* mo, double norm){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
      global_sum(motmp->corr_re,motmp->corr_size);
      global_sum(motmp->corr_im,motmp->corr_size);
      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] *= norm;
	motmp->corr_im[i] *= norm;
      }
    motmp=motmp->next;
  }
}

int compareSources(spinor_field* s1, spinor_field* s2, double tol)
{
  int c[4], beta, color;
  int res=1;
  double diffre, diffim = 0.0;

  for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++) for(c[3]=0; c[3]<Z; c[3]++) {
	  for (beta=0;beta<4;++beta) for (color=0; color<NF; ++color){
	    diffre = _FIELD_AT(&s1[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re - _FIELD_AT(&s2[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re;
	    diffim = _FIELD_AT(&s1[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].im - _FIELD_AT(&s2[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].im;
	    res= (res && (diffre*diffre < tol) && (diffim*diffim <tol));
	    if(res==0) return 0;
	  }
  }
  return 1;
}

void unit_gauge(suNg_field *gauge){
  int mu;
    int x, y, z, t, ix; for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ ix=ipt(t,x,y,z);
	    for (mu=0; mu<4; ++mu) {
		_suNg_unit(*_4FIELD_AT(gauge,ix,mu));
	    }
    }
  start_gf_sendrecv(gauge);
  complete_gf_sendrecv(gauge);
}

void resetTimes(spinor_field *out, spinor_field *in, int tau)
{
  int c[4];
  int beta, color;
  for (beta=0;beta<4;++beta){
    spinor_field_zero_f(&out[beta]);
  }
  if(COORD[0]==tau/T){
  for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++) for(c[3]=0; c[3]<Z; c[3]++) {
    for (beta=0;beta<4;++beta) for (color=0; color<NF; ++color){
      _FIELD_AT(&out[beta], ipt(tau,c[1],c[2],c[3]) )->c[beta].c[color].re = _FIELD_AT(&in[beta], ipt(tau,c[1],c[2],c[3]) )->c[beta].c[color].re;
      _FIELD_AT(&out[beta], ipt(tau,c[1],c[2],c[3]) )->c[beta].c[color].im = _FIELD_AT(&in[beta], ipt(tau,c[1],c[2],c[3]) )->c[beta].c[color].im;
    }
  }
  }
  //Not sure what this loop does, but it's in every source definition, so it must be important?
  for (beta=0;beta<4;++beta){
     start_sf_sendrecv(out + beta);
     complete_sf_sendrecv(out + beta);
  }
}
// Function for initiating meson observable (used to store the correlation function)
void init_mo(meson_observable* mo, char* name, int size)
{
  int i;
  //ind1 and ind2 don't do anything for the moment
  mo->ind1 = _g5;
  mo->ind2 = _g5;
  strcpy(mo->channel_name,name);
  strcpy(mo->channel_type, "Pi Pi scattering");
  mo->sign=1.0;
  mo->corr_size=size;
  mo->corr_re = (double * ) malloc(size * sizeof(double));
  mo->corr_im = (double * ) malloc(size * sizeof(double));
  if (mo->corr_re == NULL || mo->corr_im == NULL)
  {
    fprintf(stderr, "malloc failed in init_mo \n");
    return;
  }
  mo->next=NULL;
  for (i=0; i<size; ++i)
  {
    mo->corr_re[i]=0.0;
    mo->corr_im[i]=0.0;
  }
}

void reset_mo(meson_observable* mo)
{
  int i;
  for (i=0; i< mo->corr_size; ++i)
  {
    mo->corr_re[i]=0.0;
    mo->corr_im[i]=0.0;
  }
}

void free_mo(meson_observable* mo)
{
  free(mo->corr_re);
  free(mo->corr_im);
  free(mo);
}

int main(int argc,char *argv[])
{
  meson_observable* mo=NULL;
  int t, px, py, pz;
	int i,k;
	char tmp[256], *cptr;
	FILE* list;
	filename_t fpars;
	int nm;
	double m[256];

  //Copy I/O from another file
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  read_input(glb_var.read,input_filename);
  setup_replicas();

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0,0);
  if (PID!=0) { logger_disable(); }
  if (PID==0) { 
    sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (list_filename!=NULL) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 


  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);

  read_input(mes_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;

  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
  // geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */ 
  // test_geometry_mpi_eo();  
  /* setup random numbers */
  read_input(rlx_var.read,input_filename);
  //slower(rlx_var.rlxd_start); //convert start variable to lowercase
  if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0') {
    /*load saved state*/
    lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
    read_ranlxd_state(rlx_var.rlxd_state);
  } else {
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  }

  //                                    
  //lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed);
  //rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  //srand(glb_var.rlxd_seed+PID);

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  nm=0;
  if(fpars.type==DYNAMICAL_CNFG) {
    nm=1;
    m[0] = fpars.m;
  } else if(fpars.type==QUENCHED_CNFG) {
    strcpy(tmp,mes_var.mstring);
    cptr = strtok(tmp, ";");
    nm=0;
    while(cptr != NULL) {
      m[nm]=atof(cptr);
      nm++;
      cptr = strtok(NULL, ";");
      printf(" %3.3e \n",m[nm]);
    }            
  }



  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  init_BCs(NULL);



  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  lprintf("MAIN",0,"Inverter precision = %e\n",mes_var.precision);
  for(k=0;k<nm;k++)
  {
    lprintf("MAIN",0,"Mass[%d] = %f\n",k,m[k]);
    lprintf("CORR",0,"Mass[%d] = %f\n",k,m[k]);
  }

  read_gauge_field(cnfg_filename);
  represent_gauge_field();
  //End of the I/O block
  spinor_field* source_ts = alloc_spinor_field_f(4*NF,&glattice);
  spinor_field* source_tsp1 = alloc_spinor_field_f(4,&glattice);
  spinor_field* source_ts_mom = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop_ts =  alloc_spinor_field_f(4*NF ,&glattice);
  spinor_field* prop_tsp1 =  alloc_spinor_field_f(4 ,&glattice);
  spinor_field* prop_ts_mom =  alloc_spinor_field_f(4 ,&glattice);

  lprintf("Global vars", 0, "T=%i, T_BORDER=%i, T_EXT=%i\n",T, T_BORDER, T_EXT);

  init_propagator_eo(nm,m,mes_var.precision);

  //Create propagators
  spinor_field_zero_f(source_ts);
  create_diluted_source_equal_atau(source_ts, 0);
  calc_propagator(prop_ts,source_ts,1);

  spinor_field_zero_f(source_tsp1);
  create_diluted_source_equal_atau(source_tsp1, 1);
  calc_propagator(prop_tsp1,source_tsp1,1);

  //Add 1 unit of momentum in the z direction
/*  addMomentum(source_ts_mom, source_ts, 0, 0, 1);
  calc_propagator(prop_ts_mom,source_ts_mom,1);*/

  //Vincent's code (IO included)
  contract_pion_scatt_1spinorfield(prop_ts, prop_tsp1, 0 , 0);

  //Create a meson observable
  mo =  malloc(sizeof(meson_observable));
  init_mo(mo,"Contraction A",GLB_T);
  lprintf("DEBUG",0,"Initialise meson observable\n");
  //New code
  //Contraction A
  measure_scattering_AD_core(mo, prop_ts, prop_ts, prop_tsp1, prop_tsp1, 0, 1, 0, 0, 0, 0);
  lprintf("DEBUG",0,"New contraction\n");
  lprintf("From a",0,"%3.10e\n", mo->corr_re[0]);
  for(t=0;t<GLB_T;++t) lprintf("A",0,"%i %3.10e %3.10e \n",t,mo->corr_re[t]/GLB_VOL3/GLB_VOL3, mo->corr_im[t]/GLB_VOL3/GLB_VOL3);
  lprintf("DEBUG",0,"Print output\n");
  //Contraction B
  reset_mo(mo);
  measure_scattering_BC_core(mo, prop_ts, prop_ts, prop_tsp1, prop_tsp1, 0, 1, 0, 0, 0, 0);
  for(t=0;t<GLB_T;++t) lprintf("B",0,"%i %3.10e %3.10e \n",t,mo->corr_re[t]/GLB_VOL3/GLB_VOL3, mo->corr_im[t]/GLB_VOL3/GLB_VOL3);
  //Contraction C
  lprintf("DEBUG",0,"Starting contraction C\n");
  reset_mo(mo);
  measure_scattering_BC_core(mo, prop_ts, prop_ts, prop_tsp1, prop_tsp1, 0, -1, 0, 0, 0, 0);
  for(t=0;t<GLB_T;++t) lprintf("C",0,"%i %3.10e %3.10e \n",t,mo->corr_re[t]/GLB_VOL3/GLB_VOL3, mo->corr_im[t]/GLB_VOL3/GLB_VOL3);
  //Contraction D
  reset_mo(mo);
  measure_scattering_AD_core(mo, prop_ts, prop_ts, prop_tsp1, prop_tsp1, 0, -1, 0, 0, 0, 0);
  for(t=0;t<GLB_T;++t) lprintf("D",0,"%i %3.10e %3.10e \n",t,mo->corr_re[t]/GLB_VOL3/GLB_VOL3, mo->corr_im[t]/GLB_VOL3/GLB_VOL3);

  // Test the addMomentum function - compare existing momentum source with wall source+addMomentum; the momentum of choice here is (0,0,1)
  // The original code (Ari & Rudy)
  spinor_field_zero_f(source_ts);
  spinor_field_zero_f(source_ts_mom);
  create_gauge_fixed_momentum_source(source_ts, 0, 0, 0, 1, 0);
  //Setting all the time slices different from 0 to zero - why is it written like this anyway?!
  resetTimes(source_ts_mom, source_ts, 0);
  calc_propagator(prop_ts_mom,source_ts_mom,1);
  //Creating a wall source
  spinor_field_zero_f(source_ts);
  create_gauge_fixed_wall_source(source_ts,0,0);
  calc_propagator(prop_ts,source_ts,1);
  reset_mo(mo);
  measure_mesons_core(prop_ts,prop_ts_mom, source_ts, mo, 1, 0, 1, 0, GLB_T);
#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*GLB_T*(nm)+(py)*(n_mom)*GLB_T*(nm)+(pz)*GLB_T*(nm)+ ((cm)*GLB_T) +(tc))
  do_global_sum(mo,1.0);
  for(t=0;t<GLB_T;++t)
  {
    lprintf("Testing add_momentum - existing code",0,"%i %3.10e %3.10e \n", t, mo->corr_re[corr_ind(0,0,0,1,t,1,0)], mo->corr_im[corr_ind(0,0,0,1,t,1,0)]);
  }
  //Now try addMomentum,
  spinor_field_zero_f(source_ts_mom);
  spinor_field_zero_f(prop_ts_mom);
  add_momentum(source_ts_mom, source_ts, 0, 0, 1);
  calc_propagator(prop_ts_mom,source_ts_mom,1);
  reset_mo(mo);
  measure_mesons_core(prop_ts,prop_ts_mom, source_ts, mo, 1, 0, 1, 0, GLB_T);
  do_global_sum(mo,1.0);
  for(t=0;t<GLB_T;++t)
  {
    lprintf("Testing add_momentum - new code",0,"%i %3.10e %3.10e \n",t, mo->corr_re[corr_ind(0,0,0,1,t,1,0)], mo->corr_im[corr_ind(0,0,0,1,t,1,0)]);
  }


  //Now prepare for the correlation functions with momentum
  free_mo(mo);
  mo =  malloc(sizeof(meson_observable));
  init_mo(mo,"Contraction A",GLB_T*27);

  //Reset the gauge to unit gauge
  unit_gauge(u_gauge);
  //Create point source propagators (easy to compare with free theory results)
  spinor_field_zero_f(source_ts);
  create_point_source(source_ts, 0, 0);
  calc_propagator(prop_ts,source_ts,4);

/*  spinor_field_zero_f(source_tsp1);
  create_point_source(source_ts, 1, 0);
  calc_propagator(prop_tsp1,source_tsp1,1);*/

  reset_mo(mo);
  measure_scattering_AD_core(mo, prop_ts, prop_ts, prop_ts, prop_ts, 0, 0, 1, 0, 0, 1);
  for(px=-1;px<=1;++px) for(py=-1;py<=1;++py) for(pz=-1;pz<=1;++pz) for(t=0;t<GLB_T;++t) lprintf("A - total momentum 0 0 1",0,"%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,mo->corr_re[INDEX(px,py,pz,1,t)], mo->corr_im[INDEX(px,py,pz,1,t)]);
  reset_mo(mo);
  measure_scattering_BC_core(mo, prop_ts, prop_ts, prop_ts, prop_ts, 0, 0, 1, 0, 0, 0);
  for(px=-1;px<=1;++px) for(py=-1;py<=1;++py) for(pz=-1;pz<=1;++pz) for(t=0;t<GLB_T;++t) lprintf("B - total momentum 0 0 0",0,"%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,mo->corr_re[INDEX(px,py,pz,1,t)], mo->corr_im[INDEX(px,py,pz,1,t)]);
  reset_mo(mo);
  measure_scattering_BC_core(mo, prop_ts, prop_ts, prop_ts, prop_ts, 0, 0, 1, 0, 0, 1);
  for(px=-1;px<=1;++px) for(py=-1;py<=1;++py) for(pz=-1;pz<=1;++pz) for(t=0;t<GLB_T;++t) lprintf("B - total momentum 0 0 1",0,"%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,mo->corr_re[INDEX(px,py,pz,1,t)], mo->corr_im[INDEX(px,py,pz,1,t)]);

  //Two-point functions for comparison
  free_mo(mo);
  mo = malloc(sizeof(meson_observable));
  init_mo(mo, "Two-point", GLB_T*27);

  double c1re, c1im, c2re, c2im = 0.0;

  measure_mesons_core(prop_ts,prop_ts, source_ts, mo, 1, 0, 2, 0, GLB_T);
//  print_mesons(mo, 1.0, 0, 1, m,GLB_T,2,"TEST");
  do_global_sum(mo,1.0);
  for(t=0;t<GLB_T;++t)
  {
    c1re = mo->corr_re[corr_ind(0,0,0,2,t,1,0)];
    c1im = mo->corr_im[corr_ind(0,0,0,2,t,1,0)];
    c2re = mo->corr_re[corr_ind(0,0,1,2,t,1,0)];
    c1im = mo->corr_im[corr_ind(0,0,1,2,t,1,0)];
    lprintf("Two-point",0,"%i %3.10e %3.10e %3.10e %3.10e %3.10e %3.10e\n",t, c1re , c1im , c2re , c2im, c1re*c2re-c1im*c2im, c1re*c2im+c2re*c1im);
  }
/*  measure_point_mesons_momenta(mo, prop_ts, source_ts, 1, 0, 1);
  for(t=0;t<GLB_T;++t)
  {
    lprintf("Two-point",0,"%i %3.10e %3.10e %3.10e %3.10e\n",t, mo->corr_re[corr_ind(0,0,0,1,t,1,0)]/GLB_VOL3, mo->corr_im[corr_ind(0,0,0,1,t,1,0)]/GLB_VOL3, mo->corr_re[corr_ind(0,0,1,1,t,1,0)]/GLB_VOL3, mo->corr_im[corr_ind(0,0,1,1,t,1,0)]/GLB_VOL3);
  }*/

  spinor_field* source_seq0 = alloc_spinor_field_f(4*NF,&glattice);
  spinor_field* source_seqt = alloc_spinor_field_f(4*NF,&glattice);
  spinor_field* source_seqtmom = alloc_spinor_field_f(4*NF,&glattice);
  spinor_field* prop_seq0 =  alloc_spinor_field_f(4*NF ,&glattice);
  spinor_field* prop_seqt =  alloc_spinor_field_f(4*NF ,&glattice);
  create_sequential_source(source_seq0, 0, prop_ts);
  calc_propagator(prop_seq0,source_seq0,4);
  free_mo(mo);
  mo = malloc(sizeof(meson_observable));
  init_mo(mo, "Two-point", 27*GLB_T);
  mo->ind2=_g3;
  measure_mesons_core(prop_ts,prop_seq0, source_ts, mo, 1, 0, 2, 0, GLB_T);
  do_global_sum(mo,1.0);
  mo->ind2=_g5;
  for(t=0;t<GLB_T;++t)
  {
    lprintf("Triangle", 0, "%i  %3.10e %3.10e %3.10e %3.10e \n", t,  mo->corr_re[corr_ind(0,0,0,2,t,1,0)], mo->corr_im[corr_ind(0,0,0,2,t,1,0)],  mo->corr_re[corr_ind(0,0,1,2,t,1,0)], mo->corr_im[corr_ind(0,0,1,2,t,1,0)]);
  }
  for(t=0;t<GLB_T;++t)
  {
    reset_mo(mo);
    create_sequential_source(source_seqt, t, prop_ts);
    calc_propagator(prop_seqt,source_seqt,4);
    measure_mesons_core(prop_seq0,prop_seqt, source_ts, mo, 1, 0, 2, 0, GLB_T);
    do_global_sum(mo,1.0);
    lprintf("Rectangle", 0, "%i  %3.10e %3.10e %3.10e %3.10e \n", t,  mo->corr_re[corr_ind(0,0,0,2,t,1,0)], mo->corr_im[corr_ind(0,0,0,2,t,1,0)],  mo->corr_re[corr_ind(0,0,1,2,t,1,0)], mo->corr_im[corr_ind(0,0,1,2,t,1,0)]);
  }

  for(t=0;t<GLB_T;++t)
  {
    reset_mo(mo);
    create_sequential_source(source_seqt, t, prop_ts);
    add_momentum(source_seqtmom, source_seqt, 0,0,-1);
    calc_propagator(prop_seqt,source_seqtmom,4);
    measure_mesons_core(prop_seq0,prop_seqt, source_ts, mo, 1, 0, 2, 0, GLB_T);
    do_global_sum(mo,1.0);
    lprintf("Rectangle - momentum on sequential prop", 0, "%i  %3.10e %3.10e %3.10e %3.10e \n", t,  mo->corr_re[corr_ind(0,0,0,2,t,1,0)], mo->corr_im[corr_ind(0,0,0,2,t,1,0)],  mo->corr_re[corr_ind(0,0,1,2,t,1,0)], mo->corr_im[corr_ind(0,0,1,2,t,1,0)]);
  }

  free_mo(mo);
  finalize_process();
  free_BCs();
  free_gfield(u_gauge);
  free_propagator_eo();

  return 0;
}

