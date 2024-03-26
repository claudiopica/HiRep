/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file hr_sendrecv.h
 * @brief Header for custom communications, configurable over compilation flags
*/

#ifndef HR_SENDRECV_H
#define HR_SENDRECV_H

#ifdef WITH_MPI

void init_hr_comms();
void finalize_hr_comms();
void hr_sendrecv(void *sendbuffer, void *recvbuffer, geometry_descriptor *type, MPI_Datatype mpi_real_type, int field_dim,
                 int size_of_real, int mpi_chunks_per_site, int nbuffers, MPI_Request *field_reqs);
void hr_sendrecv_complete(int nreq, MPI_Request *field_reqs);

#endif

#endif