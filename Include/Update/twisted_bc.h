/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file twisted_bc.h
 * @brief Macros to encapsulate fermion twisting in the dirac operator.
 */

#ifndef TWISTED_BC_H
#define TWISTED_BC_H

#ifndef __cplusplus
#define _declare_vtmp(s)       \
    register suNf_vector vtmp; \
    register suNf_vector_flt vtmp_flt
// clang-format off
#define _vtmp(s) _Generic((s), suNf_vector: vtmp, suNf_vector_flt: vtmp_flt)
// clang-format on
#else
#define _declare_vtmp(s) auto vtmp = (s)
#define _vtmp(s) vtmp
#endif

#ifdef BC_T_THETA

/**
    * @brief Multiply spinor s with matrix u, apply fermion twisting in the time
    *        direction and save the result in spinor s.
    *
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, acts on spinor s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_T_multiply(r, u, s)                \
    do {                                               \
        _declare_vtmp(s);                              \
        _suNf_multiply(_vtmp(s), (u), (s));            \
        _vector_mulc_f((r), eitheta_gpu[0], _vtmp(s)); \
    } while (0)

/**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the time direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_T_inverse_multiply(r, u, s)             \
    do {                                                    \
        _declare_vtmp(s);                                   \
        _suNf_inverse_multiply(_vtmp(s), (u), (s));         \
        _vector_mulc_star_f((r), eitheta_gpu[0], _vtmp(s)); \
    } while (0)
#else

/**
    * @brief Multiply spinor s with matrix u and save the result in spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, acts on spinor s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_T_multiply(r, u, s) _suNf_multiply((r), (u), (s))

/**
    * @brief Multiply spinor s with the inverse of matrix u and save the result in
    *        spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, its inverse will be applied to s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_T_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

#ifdef BC_X_THETA

/**
    * @brief Multiply spinor s with matrix u, apply fermion twisting in the x
    *        direction and save the result in spinor s.
    *
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, acts on spinor s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_X_multiply(r, u, s)                \
    do {                                               \
        _declare_vtmp(s);                              \
        _suNf_multiply(_vtmp(s), (u), (s));            \
        _vector_mulc_f((r), eitheta_gpu[1], _vtmp(s)); \
    } while (0)
/**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the x direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_X_inverse_multiply(r, u, s)             \
    do {                                                    \
        _declare_vtmp(s);                                   \
        _suNf_inverse_multiply(_vtmp(s), (u), (s));         \
        _vector_mulc_star_f((r), eitheta_gpu[1], _vtmp(s)); \
    } while (0)
#else

/**
    * @brief Multiply spinor s with matrix u and save the result in spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, acts on spinor s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_X_multiply(r, u, s) _suNf_multiply((r), (u), (s))

/**
    * @brief Multiply spinor s with the inverse of matrix u and save the result in
    *        spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, its inverse will be applied to s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_X_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

#ifdef BC_Y_THETA

/**
    * @brief Multiply spinor s with matrix u, apply fermion twisting in the y
    *        direction and save the result in spinor s.
    *
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, acts on spinor s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_Y_multiply(r, u, s)                \
    do {                                               \
        _declare_vtmp(s);                              \
        _suNf_multiply(_vtmp(s), (u), (s));            \
        _vector_mulc_f((r), eitheta_gpu[2], _vtmp(s)); \
    } while (0)
/**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the y direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_Y_inverse_multiply(r, u, s)             \
    do {                                                    \
        _declare_vtmp(s);                                   \
        _suNf_inverse_multiply(_vtmp(s), (u), (s));         \
        _vector_mulc_star_f((r), eitheta_gpu[2], _vtmp(s)); \
    } while (0)
#else

/**
    * @brief Multiply spinor s with matrix u and save the result in spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, acts on spinor s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_Y_multiply(r, u, s) _suNf_multiply((r), (u), (s))

/**
    * @brief Multiply spinor s with the inverse of matrix u and save the result in
    *        spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, its inverse will be applied to s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_Y_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

#ifdef BC_Z_THETA

/**
    * @brief Multiply spinor s with matrix u, apply fermion twisting in the z
    *        direction and save the result in spinor s.
    *
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, acts on spinor s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_Z_multiply(r, u, s)                \
    do {                                               \
        _declare_vtmp(s);                              \
        _suNf_multiply(_vtmp(s), (u), (s));            \
        _vector_mulc_f((r), eitheta_gpu[3], _vtmp(s)); \
    } while (0)
/**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the z direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
#define _suNf_theta_Z_inverse_multiply(r, u, s)             \
    do {                                                    \
        _declare_vtmp(s);                                   \
        _suNf_inverse_multiply(_vtmp(s), (u), (s));         \
        _vector_mulc_star_f((r), eitheta_gpu[3], _vtmp(s)); \
    } while (0)
#else

/**
    * @brief Multiply spinor s with matrix u and save the result in spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, acts on spinor s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_Z_multiply(r, u, s) _suNf_multiply((r), (u), (s))

/**
    * @brief Multiply spinor s with the inverse of matrix u and save the result in
    *        spinor s.
    *
    * @param s                     input spinor, to be multiplied by matrix
    * @param u                     input matrix, its inverse will be applied to s
    * @param r                     output spinor that contains the result of the
    *                              operation.
    */
#define _suNf_theta_Z_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif
#endif
