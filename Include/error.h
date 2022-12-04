/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file error.h
 * @brief Functions and macros for error handling.
 */

#ifndef ERROR_H
#define ERROR_H

#include "cross_compilation.h"
#ifdef WITH_GPU
    #include "gpu.h"
#endif

_LANGUAGE_C

    /**
     * @brief Print message to error file defined on startup.
     *
     * @param test              Condition on whether an error should be raised.
     *                          0 for no error and continue
     *                          1 for error, stop and print error message
     * @param no                Exit Code
     *                          Value smaller than zero exits immediately with code 0.
     *                          Value larger or equal then zero exits with code given
     *                          after finalizing.
     * @param name              Function name, where the error was raised
     * @param text              Error message text
     */
    void error(int test, int no, const char *name, const char *text);

    #ifdef WITH_MPI

        /**
         * @brief Check MPI call and log error message on failure.
         *
         * @param call          Function call that should be checked.
         */
        #define CHECK_MPI(call) \
            do { \
                const int mpireturn = call; \
                if (mpireturn != MPI_SUCCESS) \
                { \
                    char message[MPI_MAX_ERROR_STRING]; \
                    int message_length; \
                    MPI_Error_string(mpireturn, message, &message_length); \
                    error(1, 1, "Error in: %s:%d, function: %s\n" \
                                "Communications call exited with code %d: %s\n", \
                                __FILE__, __LINE__, __func__, \
                                mpireturn, message); \
                } \
            } while (0)
    #endif

_LANGUAGE_C_END

#if defined(WITH_GPU) && defined(__cplusplus)

    /**
     * @brief Check CUDA call and log error message on failure 
     *        This may be more performant than CHECK_CUDA 
     *           -> TODO: replace (the other function) in the future (SAM)
     *
     * @param err           Function call that should be checked
     */
    #define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )

    /**
     * @brief Check last error after CUDA calls
     */
    #define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )

    /**
     * @brief Check CUDA call and log error message on failure
     *
     * @param cudaError_t   Error return type from CUDA call
     * @param file          File where the exception was raised
     * @param line          Line where the exception was raised
     */
    void __cudaSafeCall( cudaError_t err, const char *file, const int line );

    /**
     * @brief Check last error
     *
     * @param file          File where the last error was raised
     * @param line          Line where the last error was raised
     */
    void __cudaCheckError( const char *file, const int line );

    /**
     * @brief Check CUDA call and log error message on failure.
     *
     * @param call           Function call that should be checked.
     */
    #define CHECK_CUDA(call) \
        do {\
            const cudaError_t err1 = call; \
            if (err1 != cudaSuccess) \
            { \
                error(1, 1, "Error in: %s:%d, function: %s\n" \
                            "CUDA call exited with code %d: %s\n", \
                             __FILE__, __LINE__, __func__, \
                             err1, cudaGetErrorString(err1));\
            } \
        } while (0)

#endif
#endif
