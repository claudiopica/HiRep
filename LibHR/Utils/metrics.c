/***************************************************************************
 * Copyright (c) 2023, Sofie Martins                                        *
 * All rights reserved.                                                     *
 ***************************************************************************/

#include "utils.h"
#include "io.h"

int flops_per_site(operator_type type) {
    int flopsite = 0;

    switch (type) {
    case DPHI_CORE:
#if defined(REPR_ADJOINT)
        flopsite = 8 * NF * (7 + 8 * NF);
#else
        flopsite = 8 * NF * (7 + 16 * NF);
#endif
        break;

    case DPHI_CORE_FLT:
#if defined(REPR_ADJOINT)
        flopsite = 8 * NF * (7 + 8 * NF);
#else
        flopsite = 8 * NF * (7 + 16 * NF);
#endif
        break;

    case CPHI:
#if defined(REPR_ADJOINT)
        flopsite = 0;
#else
        flopsite = 0;
#endif
        break;

    case CPHI_FLT:
#if defined(REPR_ADJOINT)
        flopsite = 0;
#else
        flopsite = 0;
#endif
        break;

#ifdef WITH_CLOVER
    case CPHI_CORE:
        flopsite = 0;
        break;

    case CPHI_INV:
        flopsite = 0;
        break;
#endif

#ifdef WITH_EXPCLOVER
    case CPHI_CORE:
        flopsite = 0;
        break;

    case CPHI_INV:
        flopsite = 0;
        break;
#endif

#ifdef WITH_CLOVER
    case CPHI_FLT_CORE:
        flopsite = 0;
        break;

    case CPHI_INV_FLT:
        flopsite = 0;
        break;
#endif

#ifdef WITH_EXPCLOVER
    case CPHI_FLT_CORE:
        flopsite = 0;
        break;

    case CPHI_INV_FLT:
        flopsite = 0;
        break;
#endif

    case CUDA_REDUCTION:
        flopsite = 1;
        break;

    case SF_SQNORM:
        flopsite = 5 * NF + 1;
        break;

    case PLAQUETTE:
        flopsite = 6 * (2 * NG * NG * (2 * NG - 1) + NG);
        break;

    default:
        error(1, 0, __func__, "Invalid operator or FLOP count not implemented.\n");
    }

    lprintf("LA TEST", 0, "Flop per site = %d\n", flopsite);
    return flopsite;
}

int bytes_per_site(operator_type type) {
    int bytesite = 0;
    switch (type) {
    case DPHI_CORE:
#ifdef WITH_GPU
        bytesite = 36 * sizeof(suNf_vector) + 8 * sizeof(suNf);
#else
        bytesite = 40 * sizeof(suNf_vector) + 8 * sizeof(suNf);
#endif
        break;

    case DPHI_CORE_FLT:
#ifdef WITH_GPU
        bytesite = 36 * sizeof(suNf_vector_flt) + 8 * sizeof(suNf_flt);
#else
        bytesite = 40 * sizeof(suNf_vector_flt) + 8 * sizeof(suNf_flt);
#endif
        break;

    case CPHI:
#ifdef WITH_GPU
        bytesite = 44 * sizeof(suNf_vector) + 8 * sizeof(suNf) + 2 * sizeof(suNfc);
#else
        bytesite = 52 * sizeof(suNf_vector) + 8 * sizeof(suNf) + 2 * sizeof(suNfc);
#endif
        break;

    case CPHI_FLT:
#ifdef WITH_GPU
        bytesite = 44 * sizeof(suNf_vector_flt) + 8 * sizeof(suNf_flt) + 2 * sizeof(suNfc);
#else
        bytesite = 52 * sizeof(suNf_vector_flt) + 8 * sizeof(suNf_flt) + 2 * sizeof(suNfc);
#endif
        break;

    case CPHI_CORE:
#ifdef WITH_GPU
        bytesite = 8 * sizeof(suNf_vector) + sizeof(ldl_t);
#else
        bytesite = 12 * sizeof(suNf_vector) + sizeof(ldl_t);
#endif
        break;

    case CPHI_FLT_CORE:
#ifdef WITH_GPU
        bytesite = 8 * sizeof(suNf_vector_flt) + sizeof(ldl_t);
#else
        bytesite = 12 * sizeof(suNf_vector_flt) + sizeof(ldl_t);
#endif
        break;

    case CPHI_INV:
#ifdef WITH_GPU
        bytesite = 2 * sizeof(suNf_spinor) + sizeof(ldl_t);
#else
        bytesite = 3 * sizeof(suNf_spinor) + sizeof(ldl_t);
#endif
        break;

    case CPHI_INV_FLT:
#ifdef WITH_GPU
        bytesite = 2 * sizeof(suNf_spinor_flt) + sizeof(ldl_t);
#else
        bytesite = 3 * sizeof(suNf_spinor_flt) + sizeof(ldl_t);
#endif
        break;

    case CUDA_REDUCTION:
        bytesite = sizeof(double);
        break;

    case SF_SQNORM:
        bytesite = sizeof(suNf_spinor);
        break;

    case PLAQUETTE:
        bytesite = 6 * 4 * sizeof(suNg);
        break;
    }

    lprintf("LA TEST", 0, "Byte per site = %d\n", bytesite);
    lprintf("LA TEST", 0, "Dirac data movement = %e MB\n", (double)bytesite / (1024. * 1024.) * GLB_VOLUME);
    return bytesite;
}