/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of modules
*
*******************************************************************************/

#include "libhr.h"

/// Generates an array of gaussian spinor fields
/// and copy the results in the gpu memory
void setup_random_fields(int n, spinor_field s[n]) {
    for (int i=0; i<n; i++){
        gaussian_spinor_field(&s[i]);
        copy_to_gpu_spinor_field_f(&s[i]);
    }
}

static inline double spinor_max(suNf_spinor *s) {
    double *a = (double*)s;
    double max=0.;
    for (int i=0; i<sizeof(suNf_spinor)/sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max<v) max=v;
    }
    return max;
}

static double spinor_field_findmax_f(spinor_field *in) {
    double max=0.;
    _ONE_SPINOR_FOR(in) {
         suNf_spinor c = *_SPINOR_PTR(in);
         double v = spinor_max(&c);
         if (max<v) max=v;
    }
#ifdef WITH_MPI
    global_max(&max,1);
#endif
    return max;
}

int errors = 0; // count the number of errors during this test unit
static double EPSILON=1.e-14;

/// @brief  Check if the two inputs are the same within a given relative precision of EPSILON
/// @param abs1 
/// @param abs2 
static void compare_diff(double abs1, double abs2) {
    double rel = fabs(abs1-abs2)/fabs(abs1);
    const char *msg = (rel>EPSILON) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST",2,"%s rel=%.10e abs=%.10e diff=%.10e\n", msg, rel, fabs(abs1), fabs(abs1-abs2));
}

/// @brief Compare two spinor fields in the cpu and gpu parts of out by comparing the MAX and L2 norm
/// @param out Input spinor_field. Th function compare its cpu and gpu parts
/// @param diff Additional spinor_field used for scratch work space
static void compare_cpu_gpu(spinor_field *out, spinor_field *diff) {
    spinor_field_copy_f_gpu(diff, out);
    copy_from_gpu_spinor_field_f(diff);
    spinor_field_sub_assign_f_cpu(diff, out);
    double res = spinor_field_findmax_f(diff);
    double norm2 = spinor_field_sqnorm_f_cpu(diff);
    const char *msg = (res>EPSILON) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST",2,"%s MAX norm=%.10e L2 norm=%.10e\n", msg, res, sqrt(norm2));
}

#define _TEST_LIN_ALG(_name, _ninputs, _in, _out, _test) \
    do {                                      \
        setup_random_fields(_ninputs, _in);   \
        spinor_field *out=_out;               \
        spinor_field *diff=out+1;             \
        _test                                 \
        lprintf("GPU TEST",2,"%15s: ",_name); \
        compare_cpu_gpu(out, diff);           \
    } while (0)

#define _TEST_RED_OP(_name, _ninputs, _in, _test) \
    do {                                      \
        setup_random_fields(_ninputs, _in);   \
        _test                                 \
        lprintf("GPU TEST",2,"%15s: ",_name); \
        compare_diff(abs1, abs2);             \
    } while (0)

#if 0
void unit_array(double *a, int len)
{
      for(int i=0;i<len;i++)
      {
            a[i]=1.;
      }
}

void unit_spinor_field_cpu(spinor_field *s)
{
        geometry_descriptor *type = s->type;
        _PIECE_FOR(type, ixp){
             int start = type->master_start[ixp];
             int N = type->master_end[ixp] - type->master_start[ixp]+1;
             unit_array((double*)(_FIELD_AT(s, start)), N * sizeof(suNf_spinor) / sizeof(double));
        }
}

void print_spinor_field_cpu(spinor_field *s)
{
        geometry_descriptor *type = s->type;
        _PIECE_FOR(type, ixp){

             int start = type->master_start[ixp];
             int N = type->master_end[ixp] - type->master_start[ixp]+1;
             hr_complex *r = (hr_complex*)(_FIELD_AT(s, start));

             for(int i=0;i<N*sizeof(suNf_spinor)/sizeof(hr_complex);i++) {

                 printf("spinor[%d,%d] = %f, %f\n", ixp, i, creal(r[i]), cimag(r[i]));

             }
        }

        printf("\n\nDone.\n\n");
}
#endif

int main(int argc, char *argv[]) {

    /* setup process id and communications */
    //logger_setlevel(0,10000);
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    const int niter = 1;
    // Allocate memory for CPU and GPU spinor fields
    // add 2 for the output results used in the macro TEST
    int ninputs = 3; //max number of inputs
    spinor_field *in;
    in=alloc_spinor_field_f(ninputs+2, &glattice);

    for (int k=0; k<niter; k++) {
        lprintf("TEST",0,"Loop #%d\n=====================================================\n",k);
        _TEST_LIN_ALG("s1=s2", 1, in, in+1, 
            spinor_field_copy_f(out,&in[0]);
            spinor_field_copy_f_cpu(out,&in[0]);
        );

        _TEST_RED_OP("Re<s1,s2>", 2, in,
            double abs1 = spinor_field_prod_re_f(&in[0], &in[1]);
            double abs2 = spinor_field_prod_re_f_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("Im<s1,s2>", 2, in,
            double abs1 = spinor_field_prod_im_f(&in[0], &in[1]);
            double abs2 = spinor_field_prod_im_f_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("<s1,s2>", 2, in,
            hr_complex c1 = spinor_field_prod_f(&in[0], &in[1]);
            hr_complex c2 = spinor_field_prod_f_cpu(&in[0], &in[1]);
            double abs1 = _complex_prod_re(c1,c1);
            double abs2 = _complex_prod_re(c2,c2);
        );

        _TEST_RED_OP("Re<g5*s1,s2>", 2, in,
            double abs1 = spinor_field_g5_prod_re_f(&in[0], &in[1]);
            double abs2 = spinor_field_g5_prod_re_f_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("Im<g5*s1,s2>", 2, in,
            double abs1 = spinor_field_g5_prod_im_f(&in[0], &in[1]);
            double abs2 = spinor_field_g5_prod_im_f_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("|s1|^2", 1, in,
            double abs1 = spinor_field_sqnorm_f(&in[0]);
            double abs2 = spinor_field_sqnorm_f_cpu(&in[0]);
        );

        _TEST_LIN_ALG("s2+=r*s1", 2, in, in+1,
            double r = 10.0;
            spinor_field_mul_add_assign_f(out,r,&in[0]);
            spinor_field_mul_add_assign_f_cpu(out,r,&in[0]);
        );

        _TEST_LIN_ALG("s2+=c*s1", 2, in, in+1,      
            hr_complex c = 2.0+1.0*I;
            spinor_field_mulc_add_assign_f(out,c,&in[0]);
            spinor_field_mulc_add_assign_f_cpu(out,c,&in[0]);
        );

        _TEST_LIN_ALG("s2+=c*g5*s1", 2, in, in+1,
            hr_complex c = 2.0 + 3.0*I;
            spinor_field_g5_mulc_add_assign_f(out,c,&in[0]);
            spinor_field_g5_mulc_add_assign_f_cpu(out,c,&in[0]);
        );

        _TEST_LIN_ALG("s2=r*s1", 1, in, in+1,      
            double r = 10.0;
            spinor_field_mul_f(out,r,&in[0]);
            spinor_field_mul_f_cpu(out,r,&in[0]);
        );

        _TEST_LIN_ALG("s2=c*s1", 1, in, in+1,      
            hr_complex c = 2.0+1.0*I;
            spinor_field_mulc_f(out,c,&in[0]);
            spinor_field_mulc_f_cpu(out,c,&in[0]);
        );

        _TEST_LIN_ALG("s3=s1+s2", 2, in, in+2,      
            spinor_field_add_f(out,&in[0],&in[1]);
            spinor_field_add_f_cpu(out,&in[0],&in[1]);
        );

        _TEST_LIN_ALG("s3=s1-s2", 2, in, in+2,      
            spinor_field_sub_f(out,&in[0],&in[1]);
            spinor_field_sub_f_cpu(out,&in[0],&in[1]);
        );

        _TEST_LIN_ALG("s2+=s1", 2, in, in+1,      
            spinor_field_add_assign_f(out,&in[0]);
            spinor_field_add_assign_f_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s2-=s1", 2, in, in+1,      
            spinor_field_sub_assign_f(out,&in[0]);
            spinor_field_sub_assign_f_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s1=0", 0, in, in,      
            spinor_field_zero_f(out);
            spinor_field_zero_f_cpu(out);
        );

        _TEST_LIN_ALG("s2=-s1", 1, in, in+1,
            spinor_field_minus_f(out,&in[0]);
            spinor_field_minus_f_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s3=r*s1+s*s2", 2, in, in+2,
            double r = 2.0; double s = 3.0;      
            spinor_field_lc_f(out,r,&in[0],s,&in[1]);
            spinor_field_lc_f_cpu(out,r,&in[0],s,&in[1]);
        );

        _TEST_LIN_ALG("s3+=r*s1+s*s2", 3, in, in+2,
            double r = 2.0; double s = 3.0;      
            spinor_field_lc_add_assign_f(out,r,&in[0],s,&in[1]);
            spinor_field_lc_add_assign_f_cpu(out,r,&in[0],s,&in[1]);
        );

        _TEST_LIN_ALG("s3=c*s1+d*s2", 2, in, in+2,
            hr_complex c = 2.0+3.0*I; hr_complex d = 3.0+4.0*I;      
            spinor_field_clc_f(out,c,&in[0],d,&in[1]);
            spinor_field_clc_f_cpu(out,c,&in[0],d,&in[1]);
        );

        _TEST_LIN_ALG("s3+=c*s1+d*s2", 3, in, in+2,
            hr_complex c = 2.0+3.0*I; hr_complex d = 3.0+4.0*I;      
            spinor_field_clc_add_assign_f(out,c,&in[0],d,&in[1]);
            spinor_field_clc_add_assign_f_cpu(out,c,&in[0],d,&in[1]);
        );

        _TEST_LIN_ALG("s2=g5*s1", 1, in, in+1,
            spinor_field_g5_f(out,&in[0]);
            spinor_field_g5_f_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s1=g5*s1", 1, in, in, 
            spinor_field_g5_assign_f(&in[0]);
            spinor_field_g5_assign_f_cpu(&in[0]);
        );

        _TEST_LIN_ALG("s1=g0*s2", 1, in, in+1, 
            spinor_field_g0_f(out, &in[0]);
            spinor_field_g0_f_cpu(out, &in[0]);
        );

        _TEST_LIN_ALG("s1=g1*s2", 1, in, in+1, 
            spinor_field_g1_f(out, &in[0]);
            spinor_field_g1_f_cpu(out, &in[0]);
        );

        _TEST_LIN_ALG("s1=g2*s2", 1, in, in+1, 
            spinor_field_g2_f(out, &in[0]);
            spinor_field_g2_f_cpu(out, &in[0]);
        );

        _TEST_LIN_ALG("s1=g3*s2", 1, in, in+1, 
            spinor_field_g3_f(out, &in[0]);
            spinor_field_g3_f_cpu(out, &in[0]);
        );

        _TEST_LIN_ALG("lc1", 1, in, in+1, 
            double c1 = 2.1;
            spinor_field_lc1_f(c1, out, &in[0]);
            spinor_field_lc1_f_cpu(c1, out, &in[0]);
        );

        _TEST_LIN_ALG("lc2", 1, in, in+1, 
            double c1 = 2.4;
            double c2 = 4.3;
            spinor_field_lc2_f(c1, c2, out, &in[0]);
            spinor_field_lc2_f_cpu(c1, c2, out, &in[0]);
        );

        _TEST_LIN_ALG("lc3", 3, in, in+2, 
            double c1 = 2.4;
            double c2 = 4.3;
            spinor_field_lc3_f(c1, c2, &in[0], &in[1], &in[2]);
            spinor_field_lc3_f_cpu(c1, c2, &in[0], &in[1], &in[2]);
        );
    }

    free_spinor_field_f(in);
    finalize_process();

    return errors;
}
