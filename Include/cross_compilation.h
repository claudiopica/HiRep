/**
 * @file cross_compilation.h
 * @brief Macros to make cross compilation in header files simpler and more readable
 *        Use this only in header files that actually rely on multiple languages
 *        meaning: they have CUDA code and C code. Otherwise the compiler knows its
 *        C and this wrapping is not necessary
 */

/**
 * @brief For files that are cross compiled, meaning they contain C and C++ code
 *        the C functions need to be wrapped in _LANGUAGE_C and _LANGUAGE_C_END
 *        statements.
 */
 #ifdef __cplusplus
  #define _LANGUAGE_C extern "C" {
 #else
  #define _LANGUAGE_C 
 #endif

/**
 * @brief For files that are cross compiled, meaning they contain C and C++ code
 *        the C functions need to be wrapped in _LANGUAGE_C and _LANGUAGE_C_END
 *        statements.
 */
#ifdef __cplusplus
  #define _LANGUAGE_C_END }
#else
  #define _LANGUAGE_C_END
#endif

