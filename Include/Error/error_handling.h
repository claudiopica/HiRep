/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
/// Headerfile for:
/// -  error.c
/// -  print_trace.c

/**
 * @file error_handling.h
 * @brief Functions and macros for error handling.
 */

#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

//error.c
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

//print_trace.c
void print_trace(void);
void register_sighandlers(void);

//TODO: should we keep this here?
#ifdef WITH_MPI
/**
        * @brief Check MPI call and log error message on failure.
        *
        * @param call          Function call that should be checked.
        */
#define CHECK_MPI(call)                                                     \
    do {                                                                    \
        const int mpireturn = call;                                         \
        if (mpireturn != MPI_SUCCESS) {                                     \
            char message[MPI_MAX_ERROR_STRING];                             \
            int message_length;                                             \
            MPI_Error_string(mpireturn, message, &message_length);          \
            char errmsg[500];                                               \
            sprintf(errmsg, "Error in: %s:%d, function: %s\n \
                            Communications call exited with code %d: %s\n", \
                    __FILE__, __LINE__, __func__, mpireturn, message);      \
            error(1, 1, __func__, errmsg);                                  \
        }                                                                   \
    } while (0)
#endif

#ifdef __cplusplus
}
#endif
#endif