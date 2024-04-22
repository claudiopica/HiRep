/**
 * @file check_options.h
 * @brief Validation checks of the input parameters specified in the input file
 *        and error messages that give instructions how to fix the issue.
 */

#ifndef CHECK_OPTIONS_H
#define CHECK_OPTIONS_H

/* Boundary conditions */

#ifdef GAUGE_SPATIAL_TWIST

#ifndef REPR_ADJOINT
#error Twisted boundary conditions can be used only with the adjoint representation!!!
#endif

#undef BC_X_PERIODIC
#undef BC_X_ANTIPERIODIC
#undef BC_X_THETA
#define BC_X_PERIODIC
#define BC_X_ALREADY

#undef BC_Y_PERIODIC
#undef BC_Y_ANTIPERIODIC
#undef BC_Y_THETA
#define BC_Y_PERIODIC
#define BC_Y_ALREADY

#undef BC_Z_PERIODIC
#undef BC_Z_ANTIPERIODIC
#undef BC_Z_THETA
#define BC_Z_PERIODIC
#define BC_Z_ALREADY

#define PLAQ_WEIGHTS

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
#error(GAUGE_SPATIAL_TWIST) Twisted BCs cannot be used with Schroedinger functional!!!
#endif

#endif

#ifdef BC_T_SF

#undef BC_T_PERIODIC
#undef BC_T_ANTIPERIODIC
#undef BC_T_OPEN
#undef BC_T_THETA
#define BC_T_ALREADY

#define PLAQ_WEIGHTS

#ifdef BC_T_SF_ROTATED
#error(BC_T_SF) BC_T_SF and BC_T_SF_ROTATED cannot be used at the same time!!!
#endif

#endif

#ifdef BC_T_SF_ROTATED

#undef BC_T_PERIODIC
#undef BC_T_ANTIPERIODIC
#undef BC_T_OPEN
#undef BC_T_THETA
#define BC_T_ALREADY

#undef BC_X_PERIODIC
#undef BC_X_ANTIPERIODIC
#undef BC_X_THETA
#define BC_X_THETA

#undef BC_Y_PERIODIC
#undef BC_Y_ANTIPERIODIC
#undef BC_Y_THETA
#define BC_Y_THETA

#undef BC_Z_PERIODIC
#undef BC_Z_ANTIPERIODIC
#undef BC_Z_THETA
#define BC_Z_THETA

#define PLAQ_WEIGHTS

#ifdef BC_T_SF
#error(BC_T_SF_ROTATED) BC_T_SF and BC_T_SF_ROTATED cannot be used at the same time!!!
#endif

#endif

#if defined(HALFBG_SF) && !((NG == 2) && (defined(BC_T_SF) || defined(BC_T_SF_ROTATED)))
#error(HALFBG_SF) can be defined only if NG=2 and or BC_T_SF or BC_T_SF_ROTATED is used!!!
#endif

#if defined(BC_T_SF_ROTATED) && defined(UPDATE_EO)
#error BC_T_SF_ROTATED DOES NOT WORK WITH E/O PRECONDITIONING
#endif

#ifdef BC_T_ANTIPERIODIC

#ifdef BC_T_ALREADY
#error(BC_T_ANTIPERIODIC) BC_T already defined!!!
#endif
#define BC_T_ALREADY

#endif

#ifdef BC_T_OPEN

#ifdef BC_T_ALREADY
#error(BC_T_OPEN) BC_T already defined!!!
#endif
#define BC_T_ALREADY

#undef PLAQ_WEIGHTS
#define PLAQ_WEIGHTS

#endif

#ifdef BC_T_THETA

#ifdef BC_T_ALREADY
#error(BC_T_OPEN) BC_T already defined!!!
#endif
#define BC_T_ALREADY
#define FERMION_THETA

#endif

#if !defined(BC_T_ALREADY) && !defined(BC_T_PERIODIC)
#define BC_T_PERIODIC
#define BC_T_ALREADY
#endif

#ifdef BC_X_ANTIPERIODIC

#ifdef BC_X_ALREADY
#error(BC_X_ANTIPERIODIC) BC_X already defined!!!
#endif
#define BC_X_ALREADY

#endif

#ifdef BC_X_THETA

#ifdef BC_X_ALREADY
#error(BC_X_OPEN) BC_X already defined!!!
#endif
#define BC_X_ALREADY
#define FERMION_THETA

#endif

#if !defined(BC_X_ALREADY) && !defined(BC_X_PERIODIC)
#define BC_X_PERIODIC
#define BC_X_ALREADY
#endif

#ifdef BC_Y_ANTIPERIODIC

#ifdef BC_Y_ALREADY
#error(BC_Y_ANTIPERIODIC) BC_Y already defined!!!
#endif
#define BC_Y_ALREADY

#endif

#ifdef BC_Y_THETA

#ifdef BC_Y_ALREADY
#error(BC_Y_OPEN) BC_Y already defined!!!
#endif
#define BC_Y_ALREADY
#define FERMION_THETA

#endif

#if !defined(BC_Y_ALREADY) && !defined(BC_Y_PERIODIC)
#define BC_Y_PERIODIC
#define BC_Y_ALREADY
#endif

#ifdef BC_Z_ANTIPERIODIC

#ifdef BC_Z_ALREADY
#error(BC_Z_ANTIPERIODIC) BC_Z already defined!!!
#endif
#define BC_Z_ALREADY

#endif

#ifdef BC_Z_THETA

#ifdef BC_Z_ALREADY
#error(BC_Z_OPEN) BC_Z already defined!!!
#endif
#define BC_Z_ALREADY
#define FERMION_THETA

#endif

#if !defined(BC_Z_ALREADY) && !defined(BC_Z_PERIODIC)
#define BC_Z_PERIODIC
#define BC_Z_ALREADY
#endif

#if (defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)) && defined(WITH_QUATERNIONS)
#error Clover term has not yet been implemented with quaternions
#endif

#ifdef PURE_GAUGE_ANISOTROPY
#define PLAQ_WEIGHTS
#endif

#if defined(WITH_CLOVER) && defined(WITH_EXPCLOVER)
#error Exponential and standard clover term cannot be simultaneously defined
#endif

#if defined(WITH_SMEARING) && defined(WITH_EXPCLOVER)
#error Exponential lover term cannot be use simultaneously with the dirac smearing (not yet implemented)
#endif

// GPU checks
#if defined(WITH_GPU) && !defined(WITH_NEW_GEOMETRY)
#error Multi-GPU version does not work with old geometry. Please use new geometry.
#endif

#if defined(WITH_GPU) && defined(GAUGE_SON)
#error SO(N) gauge groups not yet implemented on GPU
#endif

#if NF > 3 && WITH_EXPCLOVER && WITH_GPU
#error "Exponential clover on GPU not implemented for NF>3"
#endif

#endif /* CHECK_OPTIONS_H */
