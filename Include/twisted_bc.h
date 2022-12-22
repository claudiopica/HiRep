/**
 * @file twisted_bc.h
 * @brief Macros to encapsulate fermion twisting in the dirac operator.
 */ 

#ifndef TWISTED_BC_H
#define TWISTED_BC_H

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
    #define _suNf_theta_T_multiply(r, u, s)\
        _suNf_multiply(vtmp, (u), (s));\
        _vector_mulc_f((r), eitheta_gpu[0], vtmp)

   /**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the time direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
    #define _suNf_theta_T_inverse_multiply(r, u, s)\
        _suNf_inverse_multiply(vtmp, (u), (s));\
        _vector_mulc_star_f((r), eitheta_gpu[0], vtmp)

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
    #define _suNf_theta_X_multiply(r, u, s)\
        _suNf_multiply(vtmp, (u), (s));\
        _vector_mulc_f((r), eitheta_gpu[1], vtmp)

   /**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the x direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
    #define _suNf_theta_X_inverse_multiply(r, u, s)\
        _suNf_inverse_multiply(vtmp, (u), (s));\
        _vector_mulc_star_f((r), eitheta_gpu[1], vtmp)

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
    #define _suNf_theta_Y_multiply(r, u, s)\
        _suNf_multiply(vtmp, (u), (s));\
        _vector_mulc_f((r), eitheta_gpu[2], vtmp)

   /**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the y direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
    #define _suNf_theta_Y_inverse_multiply(r, u, s)\
        _suNf_inverse_multiply(vtmp, (u), (s));\
        _vector_mulc_star_f((r), eitheta_gpu[2], vtmp)

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
    #define _suNf_theta_Z_multiply(r, u, s)\
        _suNf_multiply(vtmp, (u), (s));\
        _vector_mulc_f((r), eitheta_gpu[3], vtmp)

   /**
    * @brief Multiply spinor s with the inverse of matrix u, apply fermion 
    *        twisting in the z direction and save the result in spinor s.
    * 
    * @param s                 input spinor, to be multiplied by matrix
    * @param u                 input matrix, its inverse will be applied to s
    * @param r                 output spinor that contains the result of the 
    *                          operation.
    */
    #define _suNf_theta_Z_inverse_multiply(r, u, s)\
        _suNf_inverse_multiply(vtmp, (u), (s));\
        _vector_mulc_star_f((r), eitheta_gpu[3], vtmp)

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
