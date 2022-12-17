#include "global.h"

#ifdef ROTATED_SF
  #include "update.h"
  extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif 

static int init=1;
static spinor_field *gtmp=NULL;
static spinor_field *etmp=NULL;
static spinor_field *otmp=NULL;

/**
 * @brief Initializes the boundary conditions for fermion twisting
 */
static void init_bc_gpu(){
  #ifdef FERMION_THETA
  static int initialized=0;
  if (!initialized){
    cudaMemcpyToSymbol(eitheta_gpu, eitheta, 4*sizeof(hr_complex), 0, cudaMemcpyHostToDevice);
    CudaCheckError();
    initialized=1;
  }
  #endif
}

/**
 * @brief the following variable is used to keep trace of
 *        matrix-vector multiplication.
 *        we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

/**
 * @brief Getter for number of applications of the Dirac operator
 */
unsigned long int getMVM_gpu() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	return res;
}

/**
 * @brief Reset counter for number of applications of the Dirac operator
 */
void resetMVM_gpu() {
  MVMcounter = 0;
}

/**
 * @brief Free fields allocated for intermediate storage of field data.
 */
static void free_mem() {
  if (gtmp!=NULL) {
    free_spinor_field_f(gtmp);
    etmp=NULL;
  }
  if (etmp!=NULL) {
    free_spinor_field_f(etmp);
    etmp=NULL;
  }
  if (otmp!=NULL) {
    free_spinor_field_f(otmp);
    otmp=NULL;
  }
  init=1;
}

/**
 * @brief Allocate fields intended for storage of field data in intermediate
 *        steps
 */
static void init_Dirac() {
  if (init) {
    alloc_mem_t=GPU_MEM;

    gtmp=alloc_spinor_field_f(1, &glattice);
    etmp=alloc_spinor_field_f(1, &glat_even);
    otmp=alloc_spinor_field_f(1, &glat_odd);

    alloc_mem_t=std_mem_t;

    atexit(&free_mem);
    init=0;
  }
}
