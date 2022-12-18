/**
 * @file alloc_mpi_data.c
 * @brief Allocation of sendbuffers and communication handles necessary
 *        for MPI communications for both CPU and GPU.
 */

 //TODO: sendbuf_alloc analogy for GPU
 //TODO: Deallocation of senbuffers

#ifdef WITH_MPI

    #ifdef WITH_NEW_GEOMETRY
        #ifdef WITH_GPU
            #define _SENDBUF_ALLOC(_size, _i) \
            /*TODO: GPU sendbuf not allocated correctly. (SAM)*/ \
                f[_i].sendbuf_ptr = sendbuf_alloc((_size)*sizeof(*(f[_i].ptr))); \
                int alloc_length = (_size)*sizeof(*(f[_i].ptr))*(glattice.gsize_gauge - boxVolume(geometryBoxes)); \
                cudaMalloc((void **)&(f[_i].sendbuf_gpu_ptr), alloc_length);
        #else
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = sendbuf_alloc((_size)*sizeof(*(f[_i].ptr)));
        #endif
    #else
        #ifdef WITH_GPU
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = f[_i].ptr; \
                f[_i].sendbuf_gpu_ptr = f[_i].gpu_ptr;
        #else
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = f[_i].ptr; 
        #endif
    #endif

    /**
     * @brief Free memory allocated for MPI communications
     */
    #define _FREE_MPI_CODE \
        if (u->comm_req != NULL) \
            afree(u->comm_req) /* Deallocation of sendbuffers missing */


    #define _FREE_MPI_FIELD_DATA                                                                        \
        if (f->comm_req != NULL)                                                                        \
            afree(f->comm_req)

    #define _ALLOC_MPI_FIELD_DATA(_name, _size, _geom)                                                  \
        if (type->nbuffers_##_geom > 0)                                                                 \
        {                                                                                               \
            f->comm_req = amalloc(_n * 2 * type->nbuffers_##_geom * sizeof(MPI_Request), ALIGN);        \
            error((f->comm_req) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                          \
                "Could not allocate memory space for field (MPI)");                                     \
            for (int ix = 0; ix < _n * 2 * type->nbuffers_##_geom; ++ix)                                \
                f->comm_req[ix] = MPI_REQUEST_NULL;                                                     \
            for (int i = 1; i < _n; ++i)                                                                \
                f[i].comm_req = f[i-1].comm_req + 2 * type->nbuffers_##_geom;                           \
            for (int i = 1; i < _n; ++i)                                                                \
            {                                                                                           \
                _SENDBUF_ALLOC(_size, i);                                                               \
            }                                                                                           \
                                                                                                        \
        }                                                                                               \
        else                                                                                            \
        {                                                                                               \
            f->comm_req = NULL;                                                                         \
        }
  
#else
    /**
     * @brief Empty macro if compiled without WITH_MPI compilation flag
     */
    #define _FREE_MPI_FIELD_DATA do {} while (0)
    /**
     * @brief Empty macro if compiled without WITH_MPI compilation flag
     */ 
    #define _ALLOC_MPI_FIELD_DATA(_name, _size, _geom) do {} while (0)

#endif