#ifndef CLOVER_EXP_H
#define CLOVER_EXP_H

#ifdef WITH_EXPCLOVER

#include "libhr_core.h"
#include "Utils/factorial.h"

#ifdef __cplusplus
extern "C" {
#endif

void evaluate_sw_order(double *mass);
int get_NNexp();
int get_NN();

//C = B*A when C hermitian!
visible void forceinline _su2Nfc_times_su2Nfc_herm(suNfc *C, suNfc *B, suNfc *A) {
    // new zero component
    _suNfc_times_suNfc(C[0], B[0], A[0]);
    _suNfc_times_suNfc_assign(C[0], B[1], A[2]);

    // new one component
    _suNfc_times_suNfc(C[1], B[0], A[1]);
    _suNfc_times_suNfc_assign(C[1], B[1], A[3]);

    _suNfc_dagger(C[2], C[1]);

    // new three component
    _suNfc_times_suNfc(C[3], B[2], A[1]);
    _suNfc_times_suNfc_assign(C[3], B[3], A[3]);
}

//C = B*A
visible void forceinline _su2Nfc_times_su2Nfc(suNfc *C, suNfc *B, suNfc *A) {
    // new zero component
    _suNfc_times_suNfc(C[0], B[0], A[0]);
    _suNfc_times_suNfc_assign(C[0], B[1], A[2]);

    // new one component
    _suNfc_times_suNfc(C[1], B[0], A[1]);
    _suNfc_times_suNfc_assign(C[1], B[1], A[3]);

    // new two component
    _suNfc_times_suNfc(C[2], B[2], A[0]);
    _suNfc_times_suNfc_assign(C[2], B[3], A[2]);

    // new three component
    _suNfc_times_suNfc(C[3], B[2], A[1]);
    _suNfc_times_suNfc_assign(C[3], B[3], A[3]);
}

// C += A*B
visible void forceinline _su2Nfc_times_su2Nfc_assign(suNfc *C, suNfc *B, suNfc *A) {
    // new zero component
    _suNfc_times_suNfc_assign(C[0], B[0], A[0]);
    _suNfc_times_suNfc_assign(C[0], B[1], A[2]);

    // new one component
    _suNfc_times_suNfc_assign(C[1], B[0], A[1]);
    _suNfc_times_suNfc_assign(C[1], B[1], A[3]);

    // new two component
    _suNfc_times_suNfc_assign(C[2], B[2], A[0]);
    _suNfc_times_suNfc_assign(C[2], B[3], A[2]);

    // new three component
    _suNfc_times_suNfc_assign(C[3], B[2], A[1]);
    _suNfc_times_suNfc_assign(C[3], B[3], A[3]);
}

//C += B*A when C hermitian!
visible void forceinline _su2Nfc_times_su2Nfc_assign_herm(suNfc *C, suNfc *B, suNfc *A) {
    // new zero component
    _suNfc_times_suNfc_assign(C[0], B[0], A[0]);
    _suNfc_times_suNfc_assign(C[0], B[1], A[2]);

    // new one component
    _suNfc_times_suNfc_assign(C[1], B[0], A[1]);
    _suNfc_times_suNfc_assign(C[1], B[1], A[3]);

    _suNfc_dagger(C[2], C[1]);

    // new three component
    _suNfc_times_suNfc_assign(C[3], B[2], A[1]);
    _suNfc_times_suNfc_assign(C[3], B[3], A[3]);
}

//trace of B*A
visible void forceinline _su2Nfc_times_su2Nfc_trace(hr_complex *trace, suNfc *B, suNfc *A) {
    suNfc aux;

    _suNfc_times_suNfc(aux, B[0], A[0]);
    _suNfc_times_suNfc_assign(aux, B[1], A[2]);
    _suNfc_times_suNfc_assign(aux, B[2], A[1]);
    _suNfc_times_suNfc_assign(aux, B[3], A[3]);

    _suNfc_trace(*trace, aux);
}

//trace of the square of a 2NF hermitian matrix
visible void forceinline _su2Nfc_times_su2Nfc_trace_herm_sq(hr_complex *trace, suNfc *B) {
    suNfc aux;
    hr_complex auxtrace;

    _suNfc_times_suNfc(aux, B[1], B[2]);
    _suNfc_trace(auxtrace, aux);
    _suNfc_times_suNfc(aux, B[0], B[0]);
    _suNfc_times_suNfc_assign(aux, B[3], B[3]);

    _suNfc_trace(*trace, aux);

    *trace = *trace + 2 * auxtrace;
}

visible void forceinline _su2Nfc_unit(suNfc *A) {
    _suNfc_unit(A[0]);
    _suNfc_unit(A[3]);
    _suNfc_zero(A[1]);
    _suNfc_zero(A[2]);
}

visible void forceinline _su2Nfc_trace(hr_complex *p, suNfc *A) {
    hr_complex aux = 0.;
    _suNfc_trace(aux, A[0]);
    _suNfc_trace(*p, A[3]);
    *p = *p + aux;
}

#if (NF == 3)

visible void forceinline clover_exp_NF3(suNfc *Aplus, suNfc *expAplus, int NN) {
    suNfc A0[3], A2[4], A3[4], tmp1[4], tmp2[3];

    int i = 0, j = 0;
    hr_complex p[2 * NF - 1];
    _su2Nfc_times_su2Nfc_herm(A2, Aplus, Aplus);
    _su2Nfc_times_su2Nfc_herm(A3, A2, Aplus);
    _suNfc_unit(A0[0]);
    _suNfc_unit(A0[2]);
    _suNfc_zero(A0[1]);

    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[0], A3);
    _su2Nfc_times_su2Nfc_trace(&p[1], A3, A2);
    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[2], A2);
    _su2Nfc_trace(&p[3], A3);
    _su2Nfc_trace(&p[4], A2);

    /*
  p[0] = -p[0]/6 + p[4]*p[2]/8  + p[3]*p[3]/18 -p[4]*p[4]*p[4]/48;
  p[1] = -p[1]/5 + p[4]*p[3]/6;
  p[2] = -p[2]/4 + p[4]*p[4]/8;
  p[3] = -p[3]/3;
  p[4] = -p[4]/2;  
  */

    p[4] = -p[4] / 2;
    p[3] = -p[3] / 3;
    p[2] = -p[2] / 4;

    p[0] = p[3] * p[3] / 2 + p[4] * p[2] + (p[4] * p[4] * p[4] - p[0]) / 6;
    p[1] = -p[1] / 5 + p[4] * p[3];

    p[2] += p[4] * p[4] / 2;

    double q[2 * NF];
    for (i = 0; i < 2 * NF; i++) {
        q[i] = 0.;
    }

    double qlast;
    q[0] = inverse_fact(NN);

    for (i = NN - 1; i >= 0; i--) {
        qlast = q[2 * NF - 1];
        q[2 * NF - 1] = q[2 * NF - 2];
        for (j = 2 * NF - 2; j > 0; j--) {
            q[j] = q[j - 1] - creal(p[j]) * qlast;
        }
        q[0] = inverse_fact(i) - creal(p[0]) * qlast;
    }

    //Optimized to reduce operations!

    _suNfc_mul_add(expAplus[0], q[0], A0[0], q[1], Aplus[0]);
    _suNfc_mul(tmp1[0], q[2], A2[0]);
    _suNfc_add_assign(expAplus[0], tmp1[0]);

    _suNfc_mul_add(tmp1[0], q[3], A0[0], q[4], Aplus[0]);
    _suNfc_mul(tmp2[0], q[5], A2[0]);
    _suNfc_add_assign(tmp1[0], tmp2[0]);

    _suNfc_mul_add(expAplus[1], q[0], A0[1], q[1], Aplus[1]);
    _suNfc_mul(tmp1[1], q[2], A2[1]);
    _suNfc_add_assign(expAplus[1], tmp1[1]);

    _suNfc_mul_add(tmp1[1], q[3], A0[1], q[4], Aplus[1]);
    _suNfc_mul(tmp2[1], q[5], A2[1]);
    _suNfc_add_assign(tmp1[1], tmp2[1]);

    _suNfc_dagger(tmp1[2], tmp1[1]);

    _suNfc_mul_add(expAplus[3], q[0], A0[2], q[1], Aplus[3]);
    _suNfc_mul(tmp1[3], q[2], A2[3]);
    _suNfc_add_assign(expAplus[3], tmp1[3]);

    _suNfc_mul_add(tmp1[3], q[3], A0[2], q[4], Aplus[3]);
    _suNfc_mul(tmp2[2], q[5], A2[3]);
    _suNfc_add_assign(tmp1[3], tmp2[2]);

    _su2Nfc_times_su2Nfc_assign_herm(expAplus, A3, tmp1);
}
#endif

#if (NF == 2)

visible void forceinline clover_exp_NF2(suNfc *Aplus, suNfc *expAplus, int lNN) {
    suNfc A0[3], A2[4], tmp1[4];

    int i = 0, j = 0;
    hr_complex p[2 * NF - 1];
    _su2Nfc_times_su2Nfc_herm(A2, Aplus, Aplus);
    _suNfc_unit(A0[0]);
    _suNfc_unit(A0[2]);
    _suNfc_zero(A0[1]);

    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[0], A2);
    _su2Nfc_times_su2Nfc_trace(&p[1], Aplus, A2);
    _su2Nfc_trace(&p[2], A2);

    p[0] = -p[0] / 4 + p[2] * p[2] / 8;
    p[1] = -p[1] / 3;
    p[2] = -p[2] / 2;

    double q[2 * NF];
    for (i = 0; i < 2 * NF; i++) {
        q[i] = 0.;
    }
    double qlast;
    q[0] = inverse_fact(lNN);
    for (i = lNN - 1; i >= 0; i--) {
        qlast = q[2 * NF - 1];
        q[2 * NF - 1] = q[2 * NF - 2];
        for (j = 2 * NF - 2; j > 0; j--) {
            q[j] = q[j - 1] - creal(p[j]) * qlast;
        }
        q[0] = inverse_fact(i) - creal(p[0]) * qlast;
    }

    //Optimized to reduce operations!

    _suNfc_mul_add(expAplus[0], q[0], A0[0], q[1], Aplus[0]);
    _suNfc_mul_add(tmp1[0], q[2], A0[0], q[3], Aplus[0]);

    _suNfc_mul_add(expAplus[1], q[0], A0[1], q[1], Aplus[1]);
    _suNfc_mul_add(tmp1[1], q[2], A0[1], q[3], Aplus[1]);

    _suNfc_dagger(tmp1[2], tmp1[1]);

    _suNfc_mul_add(expAplus[3], q[0], A0[2], q[1], Aplus[3]);
    _suNfc_mul_add(tmp1[3], q[2], A0[2], q[3], Aplus[3]);

    _su2Nfc_times_su2Nfc_assign_herm(expAplus, A2, tmp1);
}

#endif

visible void forceinline clover_exp_taylor(suNfc *Xin, suNfc *u) {
    suNfc Xk[4], tmp[4];
    _su2Nfc_unit(u);
    _su2Nfc_unit(Xk);

    int k = 1;
    int i = 0;
    double error = 0., erroraux;
    while (1) {
        _su2Nfc_times_su2Nfc(tmp, Xk, Xin);

        for (i = 0; i < 4; i++) {
            _suNfc_mul(Xk[i], 1. / k, tmp[i]);
            _suNfc_add_assign(u[i], Xk[i]);
        }

        k++;

        error = 0.;

        for (i = 0; i < 4; i++) {
            _suNfc_sqnorm(erroraux, Xk[i]);
            error += erroraux;
        }

        if (sqrt(error) < 1e-28) { break; }
    }
}

visible void forceinline clover_exp(suNfc *Aplus, suNfc *expAplus, int lNN) {
#if (NF == 2)
    clover_exp_NF2(Aplus, expAplus, lNN);
#elif (NF == 3)
    clover_exp_NF3(Aplus, expAplus, lNN);
#else
    clover_exp_taylor(Aplus, expAplus);
#endif
}

#if (NF == 3)

visible void forceinline doublehornerNF3(double *C, suNfc *A, int lNNexp) {
    suNfc A2[4], A3[4];
    hr_complex p[2 * NF - 1];

    _su2Nfc_times_su2Nfc_herm(A2, A, A);
    _su2Nfc_times_su2Nfc_herm(A3, A2, A);

    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[0], A3);
    _su2Nfc_times_su2Nfc_trace(&p[1], A3, A2);
    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[2], A2);
    _su2Nfc_trace(&p[3], A3);
    _su2Nfc_trace(&p[4], A2);

    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[0], A3);
    _su2Nfc_times_su2Nfc_trace(&p[1], A3, A2);
    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[2], A2);
    _su2Nfc_trace(&p[3], A3);
    _su2Nfc_trace(&p[4], A2);

    p[0] = -p[0] / 6 + p[4] * p[2] / 8 + p[3] * p[3] / 18 - p[4] * p[4] * p[4] / 48;
    p[1] = -p[1] / 5 + p[4] * p[3] / 6;
    p[2] = -p[2] / 4 + p[4] * p[4] / 8;
    p[3] = -p[3] / 3;
    p[4] = -p[4] / 2;

    int i, j, k;
    double q[2 * NF], qlast;

    double **q2 = (double **)malloc((lNNexp + 1) * sizeof(double *));
    for (int l = 0; l < lNNexp + 1; ++l) {
        q2[l] = (double *)malloc((2 * NF) * sizeof(double));
    }

    for (j = 0; j <= lNNexp; j++) {
        q[0] = inverse_fact(lNNexp + 2);
        for (k = 1; k < 2 * NF; k++) {
            q[k] = 0.;
        }

        for (i = lNNexp - j; i > -1; i--) {
            qlast = q[2 * NF - 1];
            q[2 * NF - 1] = q[2 * NF - 2];
            for (k = 2 * NF - 2; k > 0; k--) {
                q[k] = q[k - 1] - creal(p[k]) * qlast;
            }
            q[0] = -creal(p[0]) * qlast + inverse_fact(i + j + 1);
        }

        for (i = 0; i < 2 * NF; i++) {
            q2[j][i] = q[i];
        }
    }

    for (i = 0; i < 2 * NF; i++) {
        q[0] = q2[lNNexp][i];
        for (k = 1; k < 2 * NF; k++) {
            q[k] = 0.;
        }

        for (j = lNNexp - 1; j > -1; j--) {
            qlast = q[2 * NF - 1];
            q[2 * NF - 1] = q[2 * NF - 2];
            for (k = 2 * NF - 2; k > 0; k--) {
                q[k] = q[k - 1] - creal(p[k]) * qlast;
            }
            q[0] = -creal(p[0]) * qlast + q2[j][i];
        }

        for (j = 0; j < 2 * NF; j++) {
            //  INDX = 2*NF*j + i
            C[2 * NF * j + i] = q[j];
        }
    }

    for (int l = 0; l < lNNexp + 1; ++l) {
        free(q2[l]);
    }
    free(q2);
}

#endif

#if (NF == 2)

visible void forceinline doublehornerNF2(double *C, suNfc *A, int lNNexp) {
    suNfc A2[4];
    hr_complex p[2 * NF - 1];
    _su2Nfc_times_su2Nfc_herm(A2, A, A);

    _su2Nfc_times_su2Nfc_trace_herm_sq(&p[0], A2);
    _su2Nfc_times_su2Nfc_trace(&p[1], A, A2);
    _su2Nfc_trace(&p[2], A2);

    p[0] = -p[0] / 4 + p[2] * p[2] / 8;
    p[1] = -p[1] / 3;
    p[2] = -p[2] / 2;

    int i, j, k;
    double q[2 * NF], qlast;

    double *q2 = (double *)malloc((lNNexp + 1) * 2 * NF * sizeof(double *));

    for (j = 0; j <= lNNexp; j++) {
        q[0] = inverse_fact(lNNexp + 2);
        for (k = 1; k < 2 * NF; k++) {
            q[k] = 0.;
        }

        for (i = lNNexp - j; i > -1; i--) {
            qlast = q[2 * NF - 1];
            q[2 * NF - 1] = q[2 * NF - 2];
            for (k = 2 * NF - 2; k > 0; k--) {
                q[k] = q[k - 1] - creal(p[k]) * qlast;
            }
            q[0] = -creal(p[0]) * qlast + inverse_fact(i + j + 1);
        }

        for (i = 0; i < 2 * NF; i++) {
            q2[2 * NF * j + i] = q[i];
        }
    }

    for (i = 0; i < 2 * NF; i++) {
        q[0] = q2[2 * NF * lNNexp + i];
        for (k = 1; k < 2 * NF; k++) {
            q[k] = 0.;
        }

        for (j = lNNexp - 1; j > -1; j--) {
            qlast = q[2 * NF - 1];
            q[2 * NF - 1] = q[2 * NF - 2];
            for (k = 2 * NF - 2; k > 0; k--) {
                q[k] = q[k - 1] - creal(p[k]) * qlast;
            }
            q[0] = -creal(p[0]) * qlast + q2[2 * NF * j + i];
        }

        for (j = 0; j < 2 * NF; j++) {
            //  INDX = 2*NF*j + i
            C[2 * NF * j + i] = q[j];
        }
    }

    free(q2);
}
#endif

visible void forceinline doublehorner(double *C, suNfc *A, int lNNexp) {
#if (NF == 3)
    doublehornerNF3(C, A, lNNexp);
#elif (NF == 2)
    doublehornerNF2(C, A, lNNexp);
#else
    // TODO: this does not work because error is not a host
    // device function. There is now a compile time
    // check that forbids NF > 3 to compile WITH_EXPCLOVER
    //error(0, 1, "doublehorner " __FILE__, "Force only implemented for NF=2 and NF=3");
#endif
}

visible void forceinline factorialCoef(double *C, int lNNexp) {
    int i, j;

    for (j = 0; j < lNNexp; j++) {
        for (i = 0; i < lNNexp; i++) {
            if (i + j <= lNNexp) {
                C[(lNNexp)*i + j] = inverse_fact(i + j + 1);
            } else {
                C[(lNNexp)*i + j] = 0.;
            }
        }
    }
}

#ifdef __cplusplus
}
#endif
#endif
#endif
