/**
 * @file alloc_cpu_field_data.c
 * @brief Allocation code snippets for the CPU/host that can be added to allocation 
 *        functions and are general to the field type.
 */

/**
 * @brief Allocate the poiner to the spinor field structure.
 *
 * @param _name                     Name of the field, this should describe
 *                                  the function this field has in the code
 *                                  given the dimension (_size)
 */
#define _ALLOC_FIELD_STRUCT(_name)                                                                          \
    f = amalloc(sizeof(*f), ALIGN);                                                                         \
        error(f == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                              \
                    "Could not allocate memory space for field (structure)");                               \
        f->type = type;

/**
 * @brief Allocate space for the field data of the local lattice on the CPU/host.
 *
 * @param _name                     Name of the field, this should describe
 *                                  the function this field has in the code
 *                                  given the dimension (_size)
 * @param _size                     Number of elementary site types saved
 *                                  per site
 * @param _geom                     The geometries are slightly different depending
 *                                  on whether the field is located on the sites or
 *                                  on the links. Fields that are located
 *                                  on the sites are "spinor-like", so put 
 *                                  'spinor' (without quotation marks), fields
 *                                  that are located on the links are "gauge-like"
 *                                  so put 'gauge'.
 */
#define _ALLOC_CPU_FIELD_DATA(_name, _size, _geom)                                                          \
    if (alloc_mem_t & CPU_MEM)                                                                              \
    {                                                                                                       \
        /* For spinors: Allocate for all spinor array elements */                                           \
        int bytes_per_site = sizeof(*(f->ptr));                                                             \
        int number_of_sites = _n * (_size) * type->gsize_##_geom;                                           \
        int field_size = bytes_per_site * number_of_sites;                                                  \
                                                                                                            \
        f->ptr = amalloc(field_size, ALIGN);                                                                \
        error((f->ptr) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                       \
                    "Could not allocate memory space for field (data)");                                    \
                                                                                                            \
        /* For spinors: Map the elements of the spinor arrays to the */                                     \
        /* starting points in the previously allocated space. */                                            \
        for (int i = 1; i < _n; ++i)                                                                        \
            f[i].ptr = f[i-1].ptr + type->gsize_##_geom * (_size);                                          \
    }                                                                                                       \
    else                                                                                                    \
        for (int i = 0; i < _n; ++i)                                                                        \
            f[i].ptr = NULL;
            