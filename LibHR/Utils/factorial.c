/***************************************************************************
 * Copyright (c) 2023, Sofie Martins                                        *
 * All rights reserved.                                                     *
 ***************************************************************************/

#include "libhr_core.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

#if MAX_FACTORIAL > 50
#error "MAX_FACTORIAL cannot be larger than 50. There is probably no reason, to calculate inverse factorials this high."
#endif

visible double inverse_fact(int i) {
    static const double inv_fact[51] = { 1.0,
                                         1.0,
                                         0.5,
                                         0.16666666666666666,
                                         0.041666666666666664,
                                         0.008333333333333333,
                                         0.001388888888888889,
                                         0.0001984126984126984,
                                         2.48015873015873e-05,
                                         2.7557319223985893e-06,
                                         2.755731922398589e-07,
                                         2.505210838544172e-08,
                                         2.08767569878681e-09,
                                         1.6059043836821613e-10,
                                         1.1470745597729725e-11,
                                         7.647163731819816e-13,
                                         4.779477332387385e-14,
                                         2.8114572543455206e-15,
                                         1.5619206968586225e-16,
                                         8.22063524662433e-18,
                                         4.110317623312165e-19,
                                         1.9572941063391263e-20,
                                         8.896791392450574e-22,
                                         3.8681701706306835e-23,
                                         1.6117375710961184e-24,
                                         6.446950284384474e-26,
                                         2.4795962632247972e-27,
                                         9.183689863795546e-29,
                                         3.279889237069838e-30,
                                         1.1309962886447718e-31,
                                         3.7699876288159054e-33,
                                         1.2161250415535181e-34,
                                         3.800390754854744e-36,
                                         1.151633562077195e-37,
                                         3.387157535521162e-39,
                                         9.67759295863189e-41,
                                         2.688220266286636e-42,
                                         7.265460179153071e-44,
                                         1.911963205040282e-45,
                                         4.902469756513544e-47,
                                         1.225617439128386e-48,
                                         2.9893108271424046e-50,
                                         7.117406731291439e-52,
                                         1.6552108677421951e-53,
                                         3.7618428812322616e-55,
                                         8.359650847182804e-57,
                                         1.817315401561479e-58,
                                         3.866628513960594e-60,
                                         8.055476070751238e-62,
                                         1.643974708316579e-63,
                                         3.2879494166331584e-65 };
    return inv_fact[i];
}

#ifdef __cplusplus
}
#endif
