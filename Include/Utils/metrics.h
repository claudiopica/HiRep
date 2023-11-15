/***************************************************************************
 * Copyright (c) 2023, Sofie Martins                                        *
 * All rights reserved.                                                     *
 ***************************************************************************/

typedef enum {
    DPHI_CORE = 0,
    DPHI_CORE_FLT = 1,
    CPHI = 2,
    CPHI_FLT = 3,
    CPHI_CORE = 4,
    CPHI_FLT_CORE = 5,
    CPHI_INV = 6,
    CPHI_INV_FLT = 7,
    CUDA_REDUCTION = 8,
    SF_SQNORM = 9,
    PLAQUETTE = 10
} operator_type;

int flops_per_site(operator_type);
int bytes_per_site(operator_type);