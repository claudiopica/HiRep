/*******************************************************************************
*
* Wrapper functions for different type of measurements 
*
*******************************************************************************/

#include "observables.h"
#include "libhr_core.h"
#include "io.h"
#include "memory.h"
#include "utils.h"

#if NG == 3
void measure_baryons(double *m, int conf_num, double precision, storage_switch swc, data_storage_array **ret) {
    // declare point sources and props
    spinor_field *source =
        alloc_spinor_field(4 * NF, &glattice); //This isn't glat_even so that the odd sites will be set to zero explicitly
    spinor_field *prop = alloc_spinor_field(4 * NF, &glattice);
    int nm = 1;
    int tau = 0;

    init_propagator_eo(nm, m, precision);

    // create point source
    create_full_point_source(source, tau);

    // to test gauge invariance of the correlator  :
    //	random_gauge_transform(u_gauge);
    //  represent_gauge_field();

    // calc invert
    calc_propagator(prop, source, 4 * NF); //4x3 for QCD

    // perform contraction
    contract_baryons(prop, tau, swc, ret);

    // one should run the meson contraction as well here
    //for (k=0;k<NF;++k){
    //	for (beta=0;beta<4;++beta){
    //		zero_spinor_field(&source2[beta]);
    //		zero_spinor_field(&prop2[beta]);
    //	}

    //		measure_mesons(prop2, source2, nm, tau);
    //}
    //print_mesons(1.,conf_num,0,m,"DEFAULT_TEST_POINT");

    // free
    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(prop);
}
#endif