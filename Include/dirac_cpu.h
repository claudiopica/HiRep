/**
 * @file dirac_cpu.h
 * @brief Implementation of the Dirac operator for calculation on CPUs
 */

#ifndef DIRAC_CPU_H
#define	DIRAC_CPU_H

#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

/**
 * @brief Number of times the dirac operator was applied to the CPU field data copy of a
 *        spinor field. Note, that this is internally implemented to count every application
 *        to a half spinor (either odd part or even part) and then divided by two. An
 *        application to a spinor defined on the full lattice is internally counted as two
 *        operations and then divided by two.
 *
 * @return unsigned long int    Number of applications of the Dirac operator to the full
 *                              lattice
 */
unsigned long int getMVM_cpu();
/**
* @brief Massless Dirac oprator (CPU version)
*
* @param out                    Spinor field that the operator will be applied to
* @param in                     Input spinor field before Dirac operator operation
*/
void Dphi_cpu_(spinor_field *out, spinor_field *in);

/**
 * @brief Dirac operator (CPU version)
 *
 * @param m0                    Mass parameter
 * @param out                   Out spinor field that the Dirac operator will be applied to
 * @param in                    Input spinor field before Dirac operator operation
 */
void Dphi_cpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Hermitian Dirac operator (CPU version)
 *
 * @param m0                    Mass parameter
 * @param out                   Out spinor field that the Hermitian Dirac operator will
 *                              be applied to
 * @param in                    Input spinor field before Dirac operator operation
 */
void g5Dphi_cpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Squared Hermitian Dirac operator (CPU version)
 *
 * @param m0                    Mass parameter
 * @param out                   Out spinor field that the squared Hermitian Dirac operator
 *                              will be applied to
 * @param in                    Input spinor field before operation
 */
void g5Dphi_sq_cpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Single precision massless Dirac operator (CPU version)
 *
 * @param out                   Single precision output spinor field that the operator will
 *                              be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void Dphi_flt_cpu_(spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision Dirac operator (CPU version)
 *
 * @param m0                    Mass parameter
 * @param out                   Single precision output spinor field that the operator will
 *                              be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void Dphi_flt_cpu(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision Hermitian Dirac operator (CPU version)
 *
 * @param m0                    Mass parameter
 * @param out                   Single precision output spinor field that the operator
 *                              will be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void g5Dphi_flt_cpu(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision squared Hermitian Dirac operator (CPU version)
 *
 * @param m0                    Mass parameter
 * @param out                   Single precision output spinor field that the operator
 *                              will be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void g5Dphi_sq_flt_cpu(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Even-Odd Preconditioned application of Dirac Operator, odd to even (CPU version)
 *
 * @param m0			Mass parameter
 * @param out			Even output spinor field that will contain the result of the 
 * 				operation
 * @param in			Odd input spinor field before Dirac operation application
 */
void Dphi_eopre_cpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-Odd Preconditioned application of Dirac Operator, even to odd (CPU version)
 *
 * @param m0			Mass parameter
 * @param out			Odd output spinor field that will contain the result of 
 * 				the operation
 * @param in			Even input spinor field before Dirac operation application
 */
void Dphi_oepre_cpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-Odd Preconditioned application of Hermitian Dirac Operator, odd to even
 *        (CPU version)
 *
 * @param m0			Mass parameter
 * @param out			Even output spinor field that will contain the result of the
 * 				operation
 * @param in			Odd input spinor field before Dirac operation application
 */
void g5Dphi_eopre_cpu(double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-Odd Preconditioned application of Hermitian Dirac Operator, even to odd
 *        (CPU version)
 *
 * @param m0			Mass parameter
 * @param out			Odd output spinor field that will contain the result of the
 * 				operation
 * @param in			Even input spinor field before Dirac operation application
 */
void g5Dphi_eopre_sq_cpu(double m0, spinor_field *out, spinor_field *in);

#ifdef __cplusplus
    }
#endif 
#endif
