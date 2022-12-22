/**
 * @file alloc_gpu_field_data.c
 * @brief Allocation code snippets for the GPU/device that can be added to allocation 
 *        functions and are general to the field type.
 */
#ifdef WITH_GPU

    /**
     * @brief Code snippet to free GPU field data.
     *
     * @param _name                 Name of the field, this should
     *                              describe the function this field has in
     *                              the code given the dimension (_size)
     * @param _site_type            Elementary site type from suN.h
     */
    #define _FREE_GPU_FIELD_DATA(_name, _site_type)                                                     \
        if (f->gpu_ptr != NULL)                                                                         \
            cudaFree(f->gpu_ptr);

    /**
     * @brief Code snipped to allocate GPU field data.
     *
     * @param _name                 Name of the field, this should
     *                              describe the function this field has in
     *                              the code given the dimension (_size)
     * @param _site_type            Elementary site type from suN.h
     * @param _size                 Number of elementary site types saved
     *                              per site
     * @param _geom                 The geometries are slightly different depending 
     *                              on whether the field is located on the sites
     *                              or on the links. Fields that are located
     *                              on the sites are "spinor-like", so put 
     *                              'spinor' (without quotation marks), fields
     *                              that are located on the links are "gauge-like"
     *                              so put 'gauge'.
     */
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom)                                      \
        if (alloc_mem_t & GPU_MEM)                                                                      \
        {                                                                                               \
            cudaError_t err;                                                                            \
            /* For spinors: Allocate for all spinor array elements */                                   \
            int bytes_per_site = sizeof(*(f->gpu_ptr));                                                 \
            int number_of_sites = _n * _size * type->gsize_##_geom;                                     \
            int field_size = number_of_sites * bytes_per_site;                                          \
                                                                                                        \
            err = cudaMalloc((void **)&(f->gpu_ptr), field_size);                                       \
            error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",                             \
                            "Could not allocate GPU memory space for field");                           \
                                                                                                        \
            /* For spinors: Map the elements of the spinor arrays to the */                             \
            /* starting points in the previously allocated space. */                                    \
            for (int i = 1; i < _n; ++i)                                                                \
                f[i].gpu_ptr = f[i-1].gpu_ptr + type->gsize_##_geom * _size;                             \
        }                                                                                               \
        else                                                                                            \
            for (int i = 0; i < _n; ++i)                                                                \
                f[i].gpu_ptr = NULL;

#else

    /**
     * @brief Empty macro if code compiled without flag WITH_GPU
     */
    #define _FREE_GPU_FIELD_DATA(_name, _site_type) do {} while (0)
    /**
     * @brief Empty macro if code compiled without flag WITH_GPU
     */
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom) do {} while (0)

#endif