/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of modules (Single precision)
*
*******************************************************************************/

#include "libhr.h"

/// Generates an array of gaussian spinor fields
/// and copy the results in the gpu memory
void setup_random_fields_flt(int n, spinor_field_flt s[n]) {
    for (int i=0; i<n; i++){
        gaussian_spinor_field_flt(&s[i]);
        copy_to_gpu_spinor_field_f_flt(&s[i]);
    }
}

static inline float spinor_max_flt(suNf_spinor_flt *s) {
    float *a = (float*)s;
    float max=0.;
    for (int i=0; i<sizeof(suNf_spinor_flt)/sizeof(*a); i++) {
        float v = fabs(a[i]);
        if (max<v) max=v;
    }
    return max;
}

static float spinor_field_findmax_f_flt(spinor_field_flt *in) {
    float max=0.;
    _ONE_SPINOR_FOR(in) {
         suNf_spinor_flt c = *_SPINOR_PTR(in);
         float v = spinor_max_flt(&c);
         if (max<v) max=v;
    }
#ifdef WITH_MPI
    global_max_flt(&max,1);
#endif
    return max;
}

int errors = 0; // count the number of errors during this test unit
static float EPSILON=1.e-4;

/// @brief  Check if the two inputs are the same within a given relative precision of EPSILON
/// @param abs1 
/// @param abs2 
static void compare_diff_flt(float abs1, float abs2) {
    float rel = fabs(abs1-abs2)/fabs(abs1);
    const char *msg = (rel>EPSILON) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST",2,"%s rel=%.10e abs=%.10e diff=%.10e\n", msg, rel, fabs(abs1), fabs(abs1-abs2));
}

/// @brief Compare two spinor fields in the cpu and gpu parts of out by comparing the MAX and L2 norm
/// @param out Input spinor_field. Th function compare its cpu and gpu parts
/// @param diff Additional spinor_field used for scratch work space
static void compare_cpu_gpu_flt(spinor_field_flt *out, spinor_field_flt *diff) {
    spinor_field_copy_f_flt_gpu(diff, out);
    copy_from_gpu_spinor_field_f_flt(diff);
    spinor_field_sub_assign_f_flt_cpu(diff, out);
    float res = spinor_field_findmax_f_flt(diff);
    float norm2 = spinor_field_sqnorm_f_flt_cpu(diff);
    const char *msg = (res>EPSILON) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST",2,"%s MAX norm=%.10e L2 norm=%.10e\n", msg, res, sqrt(norm2));
}

#define _TEST_LIN_ALG(_name, _ninputs, _in, _out, _test) \
    do {                                      \
        setup_random_fields_flt(_ninputs, _in);   \
        spinor_field_flt *out=_out;               \
        spinor_field_flt *diff=out+1;             \
        _test                                 \
        lprintf("GPU TEST",2,"%15s: ",_name); \
        compare_cpu_gpu_flt(out, diff);           \
    } while (0)

#define _TEST_RED_OP(_name, _ninputs, _in, _test) \
    do {                                      \
        setup_random_fields_flt(_ninputs, _in);   \
        _test                                 \
        lprintf("GPU TEST",2,"%15s: ",_name); \
        compare_diff_flt(abs1, abs2);             \
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
    spinor_field_flt *in;
    in=alloc_spinor_field_f_flt(ninputs+2, &glattice);

    for (int k=0; k<niter; k++) {
        lprintf("TEST",0,"Loop #%d\n=====================================================\n",k);
        _TEST_LIN_ALG("s2=r*s1", 1, in, in+1,      
            float r = 10.0;
            spinor_field_mul_f_flt(out,r,&in[0]);
            spinor_field_mul_f_flt_cpu(out,r,&in[0]);
        );

        _TEST_LIN_ALG("s2=c*s1", 1, in, in+1,      
            hr_complex_flt c = 2.0+1.0*I;
            spinor_field_mulc_f_flt(out,c,&in[0]);
            spinor_field_mulc_f_flt_cpu(out,c,&in[0]);
        );

        _TEST_LIN_ALG("s2+=r*s1", 2, in, in+1,
            float r = 10.0;
            spinor_field_mul_add_assign_f_flt(out,r,&in[0]);
            spinor_field_mul_add_assign_f_flt_cpu(out,r,&in[0]);
        );

        _TEST_LIN_ALG("s2+=c*s1", 2, in, in+1,      
            hr_complex_flt c = 2.0+1.0*I;
            spinor_field_mulc_add_assign_f_flt(out,c,&in[0]);
            spinor_field_mulc_add_assign_f_flt_cpu(out,c,&in[0]);
        );

        _TEST_LIN_ALG("s3=s1+s2", 2, in, in+2,      
            spinor_field_add_f_flt(out,&in[0],&in[1]);
            spinor_field_add_f_flt_cpu(out,&in[0],&in[1]);
        );

        _TEST_LIN_ALG("s3=s1-s2", 2, in, in+2,      
            spinor_field_sub_f_flt(out,&in[0],&in[1]);
            spinor_field_sub_f_flt_cpu(out,&in[0],&in[1]);
        );

        _TEST_LIN_ALG("s2+=s1", 2, in, in+1,      
            spinor_field_add_assign_f_flt(out,&in[0]);
            spinor_field_add_assign_f_flt_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s2-=s1", 2, in, in+1,      
            spinor_field_sub_assign_f_flt(out,&in[0]);
            spinor_field_sub_assign_f_flt_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s1=0", 0, in, in,      
            spinor_field_zero_f_flt(out);
            spinor_field_zero_f_flt_cpu(out);
        );

        _TEST_LIN_ALG("s2=-s1", 1, in, in+1,
            spinor_field_minus_f_flt(out,&in[0]);
            spinor_field_minus_f_flt_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s3=r*s1+s*s2", 2, in, in+2,
            float r = 2.f; float s = 3.f;      
            spinor_field_lc_f_flt(out,r,&in[0],s,&in[1]);
            spinor_field_lc_f_flt_cpu(out,r,&in[0],s,&in[1]);
        );

        _TEST_LIN_ALG("s3+=r*s1+s*s2", 3, in, in+2,
            float r = 2.f; float s = 3.f;      
            spinor_field_lc_add_assign_f_flt(out,r,&in[0],s,&in[1]);
            spinor_field_lc_add_assign_f_flt_cpu(out,r,&in[0],s,&in[1]);
        );

        _TEST_LIN_ALG("s3=c*s1+d*s2", 2, in, in+2,
            hr_complex_flt c = 2.f+3.f*I; hr_complex_flt d = 3.f+4.f*I;      
            spinor_field_clc_f_flt(out,c,&in[0],d,&in[1]);
            spinor_field_clc_f_flt_cpu(out,c,&in[0],d,&in[1]);
        );

        _TEST_LIN_ALG("s3+=c*s1+d*s2", 3, in, in+2,
            hr_complex_flt c = 2.f+3.f*I; hr_complex_flt d = 3.f+4.f*I;      
            spinor_field_clc_add_assign_f_flt(out,c,&in[0],d,&in[1]);
            spinor_field_clc_add_assign_f_flt_cpu(out,c,&in[0],d,&in[1]);
        );

        _TEST_LIN_ALG("s2=g5*s1", 1, in, in+1,
            spinor_field_g5_f_flt(out,&in[0]);
            spinor_field_g5_f_flt_cpu(out,&in[0]);
        );

        _TEST_LIN_ALG("s2+=c*g5*s1", 2, in, in+1,
            hr_complex_flt c = 2.f + 3.f*I;
            spinor_field_g5_mulc_add_assign_f_flt(out,c,&in[0]);
            spinor_field_g5_mulc_add_assign_f_flt_cpu(out,c,&in[0]);
        );

        _TEST_RED_OP("|s1|^2", 1, in,
            float abs1 = spinor_field_sqnorm_f_flt(&in[0]);
            float abs2 = spinor_field_sqnorm_f_flt_cpu(&in[0]);
        );

        _TEST_RED_OP("Re<s1,s2>", 2, in,
            float abs1 = spinor_field_prod_re_f_flt(&in[0], &in[1]);
            float abs2 = spinor_field_prod_re_f_flt_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("Im<s1,s2>", 2, in,
            float abs1 = spinor_field_prod_im_f_flt(&in[0], &in[1]);
            float abs2 = spinor_field_prod_im_f_flt_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("<s1,s2>", 2, in,
            hr_complex c1 = spinor_field_prod_f_flt(&in[0], &in[1]);
            hr_complex c2 = spinor_field_prod_f_flt_cpu(&in[0], &in[1]);
            float abs1 = _complex_prod_re(c1,c1);
            float abs2 = _complex_prod_re(c2,c2);
        );

        _TEST_RED_OP("Re<g5*s1,s2>", 2, in,
            float abs1 = spinor_field_g5_prod_re_f_flt(&in[0], &in[1]);
            float abs2 = spinor_field_g5_prod_re_f_flt_cpu(&in[0], &in[1]);
        );

        _TEST_RED_OP("Im<g5*s1,s2>", 2, in,
            float abs1 = spinor_field_g5_prod_im_f_flt(&in[0], &in[1]);
            float abs2 = spinor_field_g5_prod_im_f_flt_cpu(&in[0], &in[1]);
        );
    }

    free_spinor_field_f_flt(in);
    finalize_process();

    return errors;
}
