@page data_structures Data Structures
[TOC]
# Elementary Data Types

## Definition
The following data types are defined as C structures containing a single array of elements. These can be used to create fields in `Include/spinor_field.h` using the macro `_DECLARE_FIELD_STRUCT`. 


| Name                   | Array of                 | Size                                             |
|:-----------------------|:-------------------------|:-------------------------------------------------|
| `hr_complex`           | `double`			            | 2				                                         |
| `suNg_vector`          | `hr_complex`		          | \f$N_c\f$			                                   |
| `suNg`	               | `hr_complex`           	| \f$N_c\times N_c\f$                              |
| `suNg_spinor`          | `suNg_vector`            | 4                                                |
| `suNg_algebra_vector`  | `double`                 | \f$N_c\times N_c -1\f$                           |
| `suNf_vector`          | `hr_complex`             | \f$D_{\mathrm{R}}\f$                             |
| `suNf`                 | `hr_complex` or `double` | \f$D_{\mathrm{R}}\times D_{\mathrm{R}}\f$        |
| `suNfc`                | `hr_complex`             | \f$D_{\mathrm{R}}\times D_{\mathrm{R}}\f$        |
| `suNf_spinor`          | `suNf_vector`            | 4                                                |
| `ldl_t`                | 2 arrays, `hr_complex`   | each \fD_{\mathrm{R}}(2D_{\mathrm{R}}+1)         |

Here \f$N_c\f$ corresponds to the number of colors and \f$D_{\mathrm{R}}\f$ is the dimension of the fermion representation. The data type `suNf` can be real or complex depending on the representation being real or complex. 

Every data type has a corresponding single precision data type, the name of which is obtained adding the suffix `_flt`.

## Operations

Linear algebra operations require us to loop over the elements in the arrays of the structures. However, in order to achieve best performance, `for`-loops in bottleneck functions should be unrolled. As a result, these linear algebra functions have to be defined as macros, that expand the unrolled code. Since the number of iterations in the `for`-loop depends for many of the structures above on the number of colors and dimension of fermion representation, which have to be known at compile time. As a result, the definition of these macros depends on the compilation parameters. 

In order to have different macros available depending on compilation parameters, we autogenerate them using a `perl`-script in `Make/Utils`. 


TODO: Add more details on this somewhere else

The resulting macros can act both on single and double precision types. A list of macros can be found in the corresponding page in the function reference manual.

TODO: Since these are auto-generated, you have to first compile HiRep and then compile the function reference. This needs to be pointed out somewhere.

# Field Data Types

In ```HiRep```, field data is stored in field ```struct```s that contain an array of values on sites or links that will be allocated on the CPU and, if compiled with GPU acceleration, one that will be allocated on the GPU. The definitions of different fields are defined in ```LibHR/spinor_field.h```. Available types are

| Name               | Elementary Data Types        |
|--------------------|------------------------------|
| `suNg_field`       | `suNg`                       |
| `suNg_scalar_field`| `suNg_vector`                |
| `suNf_field`       | `suNf`                       |
| `spinor_field`     | `suNf_spinor`                |
| `suNg_av_field`    | `suNg_algebra_vector`        |
| `scalar_field`     | `double`                     |
| `ldl_field`        | `ldl_t`                      |
| `suNfc_field`      | `suNfc`                      |

... plus corresponding single precision types.

New field types can be declared by using the macro

```c
#define _DECLARE_FIELD_STRUCT(_name, _type) \
  typedef struct _name                   \
  {                                         \
    _type *ptr;                             \
    geometry_descriptor *type;              \
    _MPI_FIELD_DATA                         \
    _GPU_FIELD_DATA(_type)                  \
  } _name
```

The ```_name``` will define the field's new name, which can be anything, while the ```_type``` variable has to refer to a type that was defined in `suN_types.h`, listed in the previous section. ```_type``` defines the types of values on the lattice sites.

The field value copy of the CPU is defined by ```_type *ptr```, which is a 1D array containing the field's values at the lattice sites. The GPU copy is hidden behind the macro ```_GPU_FIELD_DATA(_type)```.

```c
#define _GPU_FIELD_DATA(_type)
#ifdef WITH_GPU
#undef _GPU_FIELD_DATA
#define _GPU_FIELD_DATA(_type) _type *gpu_ptr;
#endif //WITH_MPI
```

We need this macro instead of outright declaring the copy because we do not want to have a GPU copy in the field ```struct```s if we are only compiling for CPU. As can be seen from the macro ```_GPU_FIELD_DATA(_type)``` is defined to return nothing, but in the case of compilation with GPUs, it is overwritten to give a 1D array called ```gpu_ptr```, which can later be allocated on and accessed from the device.

Since memory access patterns have a high impact on application performance, the way that field data is stored on the GPU is different from how it is stored on the CPU in several ways that will be explained in the following. Further, in ```HiRep``` memory is managed manually instead of using a unified memory setup, which implies that from a kernel, only pointers to sites will be available but not the complete field structures. This has an impact on which functions and macros that work on the CPU are available to call from a CUDA kernel.

This means, that if we declare a spinor field

```c
spinor_field *s;
```

we may access its geometry description and sites on the CPU from a regular host function

```c
int main(void)
{
    spinor_field *s;

    // Query the value at the site with index 0
    suNf_spinor *field_value = s->ptr;

    // Check, whether the spinor field is odd
    if (s->type == &glat_odd) printf("Spinor is odd.\n")

    suNf_spinor *gpu_field_value = s->gpu_ptr;

    // The following fails, because it points to memory allocated on the GPU
    // and is therefore unavailable from the host.
    suNf_vector spinor_comp = (*gpu_field_value).c[0];
}
```

In a kernel, it is impossible to check whether the spinor is even or odd. Every call to the spinor field structure will fail.

```c
__global__ void example_kernel(spinor_field *s)
{
    // This fails because s is a host pointer, unless it was transferred
    // before being passed to the kernel.
    suNf_spinor field_value = *(s->ptr);

    // This fails because the geometry descriptor is saved on the host
    if (s->type == &glat_odd) printf("Spinor is odd.\n");

    // This fails, because s is located on the host and it is accessed in
    // order to access the field
    suNf_spinor *gpu_field_value = s->gpu_ptr;
}
```

The correct way to run a kernel that operates on the GPU field data copy is to pass the first site in the copy to the kernel and then access other sites. For example

```c
__global__ void example_kernel(suNf_spinor *start)
{
    int ix = blockIdx.x * blockDim.x  + threadIdx.x;
    // Get site with index ix
    suNf_spinor *site = start+ix;
}
```

The index in the 1D array is bijectively mapped to the coordinates in space and time.

## Special Fields

The elementary site type does not fully describe a physical field. What is missing is information on whether the field is located on the sites or the links and what its dimension is, i.e. how many elements belong to a single site or link. Operations on fields in `HiRep` refer to fields that are associated with a particular field type, dimension and geometry, which has a particular interpretation in physics. The following descriptors are used in `HiRep`

* `spinor_field_f`, `spinor field` located on sites, one element per site, corresponds to spinor field in physics
* `sfield`, `scalar field` located on sites, one element per site, corresponds to scalar field in physics
* `gfield`, `suNg_field` located on links, four elements per link, corresponds to a gauge field that has four directions, corresponds to gauge field in physics
* `gfield_f`, suNf_field, located on links, four elements per link, corresponds to a gauge field (four directions) in the fermion representation, corresponds to gauge field in physics
* `suNg_scalar_field`, located on links, one element per link
* `avfield`, located on links, four elements per link
* `gtransf`, `suNg` field located on links, one element per link, corresponds to a gauge transformation in physics
* `clover_ldl`, `ldl_field` located on links, one element per link
* `clover_term`, `suNfc_field` located on links, four elements per link
* `clover_force`, `suNf_field`, located on links, six elements per link

