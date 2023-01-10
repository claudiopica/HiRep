#ifdef WITH_GPU

#include "memory.h"
#include "libhr_core.h"

// TODO: inline comments
// TODO: doxygen docstrings

    /* Declare function to convert GPU to CPU format */
    /* This is necessary, because GPU performance is very sensitive to */
    /* memory access patterns, see documentation Doc/gpu_geometry.tex */
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)         \
        void to_gpu_format_##_name(_field_type *out, _field_type *in)                     \
        {                                                                                 \
            _CHECK_GEOMETRY_MATCHING(out, in);                                            \
            _PIECE_FOR(in->type, ixp) {                                                   \
                int stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1; \
                (void)stride;                                                             \
                _site_type *target  = _DFIELD_BLK(out, ixp, (_size));                     \
                _SITE_FOR(in->type, ixp, ix) {                                            \
                    int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                          \
                    for (int comp = 0; comp < _size; ++comp) {                            \
                        _site_type *source = _DFIELD_AT(in, ix, comp, (_size));           \
                        write_gpu_##_site_type(stride, (*source), target, ix_loc, comp);  \
                    }                                                                     \
                }                                                                         \
            }                                                                             \
        }

    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)         \
        void to_cpu_format_##_name(_field_type *out, _field_type *in)                     \
        {                                                                                 \
            _CHECK_GEOMETRY_MATCHING(out, in);                                            \
            _PIECE_FOR(in->type, ixp) {                                                   \
                int stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1; \
                (void)stride;                                                             \
                _site_type *source = _DFIELD_BLK(in, ixp, (_size));                       \
                _SITE_FOR(in->type, ixp, ix) {                                            \
                    int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                          \
                    for (int comp = 0; comp < _size; ++comp) {                            \
                        _site_type *target = _DFIELD_AT(out, ix, comp, (_size));          \
                        read_gpu_##_site_type(stride, (*target), source, ix_loc, comp);   \
                    }                                                                     \
                }                                                                         \
            }                                                                             \
        }

#define _DECLARE_CONVERT_FUNC(_name, _field_type, _site_type, _size) \
    _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size) \
    _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)

_DECLARE_CONVERT_FUNC(spinor_field_f, spinor_field, suNf_spinor, 1);
_DECLARE_CONVERT_FUNC(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1);
_DECLARE_CONVERT_FUNC(sfield, scalar_field, double, 1);
_DECLARE_CONVERT_FUNC(gfield, suNg_field, suNg, 4);
_DECLARE_CONVERT_FUNC(gfield_flt, suNg_field_flt, suNg_flt, 4);
_DECLARE_CONVERT_FUNC(gfield_f, suNf_field, suNf, 4);
_DECLARE_CONVERT_FUNC(gfield_f_flt, suNf_field_flt, suNf_flt, 4);
_DECLARE_CONVERT_FUNC(suNg_scalar_field, suNg_scalar_field, suNg_vector, 1);
_DECLARE_CONVERT_FUNC(avfield, suNg_av_field, suNg_algebra_vector, 4);
_DECLARE_CONVERT_FUNC(gtransf, suNg_field, suNg, 1);
_DECLARE_CONVERT_FUNC(clover_ldl, ldl_field, ldl_t, 1);
_DECLARE_CONVERT_FUNC(clover_term, suNfc_field, suNfc, 4);
_DECLARE_CONVERT_FUNC(clover_force, suNf_field, suNf, 6);

#undef _DECLARE_CONVERT_FUNC
#undef _DECLARE_CONVERT_TO_GPU_FORMAT
#undef _DECLARE_CONVERT_TO_CPU_FORMAT

#endif