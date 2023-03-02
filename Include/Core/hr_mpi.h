/**
 * @file hr_mpi.h
 * @brief Setup for mpi
 */

#ifndef HR_MPI_H
#define HR_MPI_H

//For OpenMPI we need to disable the non-standard C++ interface
//this is enabled by defaults for C++ / CUDA code
//if OpenMPI is compiled for supporting it
//can be disabled re-compiling OpenMPI with configure --disable-mpi-cxx
#ifdef WITH_MPI
#define OMPI_SKIP_MPICXX
#include <mpi.h>

#define MPI_CHECK(stmt)                                            \
    do {                                                           \
        int mpi_errno = (stmt);                                    \
        if (MPI_SUCCESS != mpi_errno) {                            \
            char mesg[MPI_MAX_ERROR_STRING];                       \
            int mesglen;                                           \
            MPI_Error_string(mpi_eerno, mesg, &mesglen);           \
            lprintf("MPI", 0, "[%s] ERROR: %s\n", __func__, mesg); \
            error(1, 1, "File: " __FILE__, " Line: " __LINE__);    \
        }                                                          \
    } while (0)

#endif

#endif
