#ifndef HR_MPI_H
#define HR_MPI_H

//For OpenMPI we need to disable the non-standard C++ interface
//this is enabled by defaults for C++ / CUDA code 
//if OpenMPI is compiled for supporting it
//can be disabled re-compiling OpenMPI with configure --disable-mpi-cxx
#ifdef WITH_MPI
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#endif

#endif
