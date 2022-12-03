/**
 * @file dirac_default.h
 * @brief Implementation of the Dirac operator (Default Functions)
 */

#ifndef DIRAC_DEFAULT_H
#define DIRAC_DEFAULT_H

#include "spinor_field.h"
#include "suN_types.h"
#include "utils.h"

/**
 * @brief Number of times the Dirac operator was applied.
 *
 * @return unsigned long int    Number of applications of the Dirac operator to the full
 *                              lattice
 */
extern unsigned long int (*getMVM) ();

/**
 * @brief nnumber of times the single precision Dirac operator was applied.
 *
 * @return unsigned long in     Number of applications of the single precision Dirac
 *                              operator to the full lattice
 */
unsigned long int getMVM_flt();

/**
 * @brief Massless Dirac operator. This will run on the GPU when compiled WITH_GPU,
 *        otherwise these correspond to calling the standard CPU version.
 *
 * @param out                   Output spinor field that the operator will be applied to
 * @param in                    Input spinor field before Dirac operator operation
 */
extern void (*Dphi_) (spinor_field *out, spinor_field *in);

/**
 * @brief Dirac operator. This will run on the GPU when compiled WITH_GPU, otherwise these
 *        correspond to calling the standard CPU version.
 *
 * @param out                   Output spinor field that the operator will be applied to
 * @param in                    Input spinor field before Dirac operator operation
 */
extern void (*Dphi) (double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Hermitian Dirac operator. This will run on the GPU when compiled WITH_GPU,
 *        otherwise these correspond to calling the standard CPU version.
 *
 * @param out                   Output spinor field that will contain the result of the
 *                              operation
 * @param in                    Input spinor field before Dirac operator operation
 */
extern void (*g5Dphi) (double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Square Hermitian Dirac operator. This will run on the GPU when compiled WITH_GPU,
 *        otherwise these correspond to calling the standard CPU version.
 *
 * @param out                   Output spinor field that will contain the result of the
 *                              operation
 * @param in                    Input spinor field before Dirac operator operation
 */
extern void (*g5Dphi_sq) (double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Single precision massless Dirac operator. This will run on the GPU when compiled
 *        WITH_GPU, otherwise these correspond to calling the standard CPU version.
 *
 * @param out                   Single precision output spinor field that the operator
 *                              will be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void Dphi_flt_(spinor_field_flt *out, spinor_field_flt *in);
/**
 * @brief Single precision Dirac operator. This will run on the GPU when compiled WITH_GPU,
 *        otherwise these correspond to calling the standard CPU version.
 *
 * @param out                   Single precision output spinor field that the operator
 *                              will be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision Hermitian Dirac operator. This will run on the GPU when compiled
 *        WITH_GPU, otherwise these correspond to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Single precision output spinor field that the operator
 *                              will be applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void g5Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Single precision squared Hermitian Dirac operator. This will run on the GPU when
 *        compiled WITH_GPU, otherwise these correspond to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Single precision output spinor field that operator will be
 *                              applied to
 * @param in                    Single precision input spinor field before Dirac operator
 *                              operation
 */
void g5Dphi_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
/**
 * @brief Even-odd preconditioned application of Dirac operator, odd to even. This will
 *        run on the GPU when compiled WITH_GPU, otherwise these correspond to calling the
 *        standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Even output spinor field that will contain the result of the
 *                              operation
 * @param in                    Odd input spinor field before Dirac operation application
 */
extern void (*Dphi_eopre) (double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-odd preconditioned application of Dirac operator, even to odd. This will
 *        run on the GPU when compiled WITH_GPU, otherwise these correspond to calling the
 *        standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Odd output spinor field that will contain the result of the
 *                              operation
 * @param in                    Even input spinor field before Dirac operation application
 */
extern void (*Dphi_oepre) (double m0, spinor_field *out, spinor_field *in);

/**
 * @brief Even-odd preconditioned application of Hermitian Dirac operator, odd to even.
 *        This will run on the GPU when compiled WITH_GPU, otherwise these correspond
 *        to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Even output spinor field that will contain the result of the
 *                              operation
 * @param in                    Odd input spinor field before Dirac operation application
 */
extern void (*g5Dphi_eopre) (double m0, spinor_field *out, spinor_field *in);
/**
 * @brief Even-odd preconditioned application of squared Hermitian Dirac operator, odd to
 *        even. This will run on the GPU when compiled WITH_GPU, otherwise these correspond
 *        to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Even output spinor field that will contain the result of the
 *                              operation
 * @param in                    Odd input spinor field before Dirac operation application
 */
extern void (*g5Dphi_eopre_sq) (double m0, spinor_field *out, spinor_field *in);


/**
 * @brief Even-odd preconditioned application of the single precision Dirac operator,
 *        odd to even. This will run on the GPU when compiled WITH_GPU, otherwise these
 *        correspond to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Even output single precision spinor field that will contain
 *                              the result of the operation
 * @param in                    Odd input single precision spinor field before Dirac
 *                              operation application
 */
void Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Even-odd preconditioned application of the single precision Dirac operator,
 *        even to odd. This will run on the GPU when compiled WITH_GPU, otherwise this
 *        corresponds to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Odd output single precision spinor field that will contain
 *                              the result of the operation
 * @param in                    Even input single precision spinor field before Dirac
 *                              operation application
 */
void Dphi_oepre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
/**
 * @brief Even-odd preconditioned application of the single precision Dirac operator,
 *        odd to even. This will run on the GPU when compiled WITH_GPU, otherwise this
 *        corresponds to calling the standard CPU version.
 *
 * @param m0                    Mass parameter
 * @param out                   Even output single precision spinor field that will contain
 *                              the result of the operation
 * @param in                    Odd input single precision spinor field before Dirac
 *                              operation application
 */
void g5Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/**
 * @brief Even-odd preconditioned application of the single precision squared
 *        Hermitian Dirac operator, odd to even. This will run on the GPU when
 *        compiled WITH_GPU, otherwise this corresponds to calling the standard CPU version.
 *
 * @param m0                    Mass paramter
 * @param out                   Even output single precision spinor field that will contain
 *                              the result of the operation
 * @param in                    Odd input single precision spinor field before Dirac
 *                              operation application
 */
void g5Dphi_eopre_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/* Dirac operators used in the Update */
//TODO: put this somewhere else (SAM)
void set_dirac_mass(double mass); // this is the mass used in the following operators
double get_dirac_mass();


#endif
