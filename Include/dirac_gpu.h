/**
 * @file dirac_gpu.h
 * @brief Implementation of the Dirac operator to run on GPUs
 */

#ifdef WITH_GPU
#ifndef DIRAC_GPU_H
#define DIRAC_GPU_H

#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

/**
 * @brief Number of times the dirac operator was applied to the GPU field data copy
 *        of a spinor field. Note, that this is internally implemented to count every
 *        application to a half spinor (either odd part or even part) and then 
 *        divided by two. An application to a spinor defined on the full lattice is 
 *        internally counted as two operations and then divided by two.
 *
 * @return unsigned long int  Number of applications of the Dirac operator to the full
 *                            lattice
 */
unsigned long int getMVM_gpu();


/**
 * @brief Massless Dirac operator (GPU version)
 *
 * @param out 			Out spinor field that the massless Dirac operator will
 * 				be applied to
 * @param in 			Input spinor field before operation
 */
void Dphi_gpu_(spinor_field *out, spinor_field *in);

/**
 * @brief Dirac operator (GPU version)
 *
 * @param out 			Out spinor field that the Dirac operator will
 * 				be applied to
 * @param in			Input spinor field before operation
 */
void Dphi_gpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Hermitian Dirac operator (GPU version)
 *
 * @param out 			Out spinor field that the Hermitian Dirac operator will
 *  				be applied to
 * @param in 			Input spinor field before operation
 */
void g5Dphi_gpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Squared Hermitian Dirac operator (GPU version)
 *
 * @param out 			Out spinor field that the squared Hermitian Dirac operator
 * 				will be applied to
 * @param in 			Input spinor field before operation
 */
void g5Dphi_sq_gpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Single precision massless Dirac operator (GPU version)
 *
 * @param out			Single precision output spinor field that the operator will
 * 				be applied to
 * @param in			Single precision input spinor field before Dirac operator
 * 				operation
 */
void Dphi_flt_gpu_(spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision Dirac operator (GPU version)
 *
 * @param m0			Mass parameter
 * @param out			Single precision output spinor field that the operator
 * 				will be applied to
 * @param in			Single precision input spinor field before Dirac operator
 * 				operation
 */
void Dphi_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision Hermitian Dirac operator (GPU version)
 *
 * @param m0 			Mass parameter
 * @param out			Single precision output spinor field that the operator will
 * 				be applied to
 * @param in			Single precision input spinor field before Dirac operator
 * 				operation
 */
void g5Dphi_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision squared Hermitian Dirac operator (GPU version)
 *
 * @param m0			Mass parameter
 * @param out			Single precision output spinor field that the operator will
 * 				be applied to
 * @param in			Single precision input spinor field before Dirac operator 
 * 				operation
 */
void g5Dphi_sq_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Even-odd preconditioned application of Dirac operator, odd to even (GPU version)
 *
 * @param m0			Mass paramter
 * @param out			Even output spinor field that will contain the result of the
 * 				operation
 * @param in			Odd input spinor field before Dirac operation application
 */
void Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-odd preconditioned application of Dirac operator, even to odd (GPU version)
 *
 * @param m0 			Mass parameter
 * @param out			Odd output spinor field that will contain the result of
 * 				the operation
 * @param in			Even input spinor field before Dirac operation application
 */
void Dphi_oepre_gpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-odd preconditioned application of Hermitian Dirac operator, odd to even
 *        (GPU version)
 *
 * @param m0			Mass parameter
 * @param out			Even output spinor field that will contain the result of the
 * 				operation
 * @param in			Odd input spinor field before Dirac operation application
 */
void g5Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-odd preconditioned application of Hermitian Dirac operator, even to odd
 *        (GPU version)
 *
 * @param m0			Mass paramter
 * @param out			Odd output spinor field that will contain the result of the 
 * 				operation
 * @param in			Even input spinor field before Dirac operation application
 */
void g5Dphi_eopre_sq_gpu(double m0, spinor_field *out, spinor_field *in);

#ifdef __cplusplus
    }
#endif
#endif
#endif


