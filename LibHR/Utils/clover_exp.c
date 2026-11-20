#include "Utils/clover_exp.h"
#include "Utils/factorial.h"
#include "io.h"
#include "update.h"

#include <math.h>
#include <float.h>

#ifdef WITH_EXPCLOVER

static int NN;
static int NNexp;

int get_NNexp() {
    return NNexp;
}

int get_NN() {
    return NN;
}

void evaluate_sw_order(double *mass) {
    static double m0 = 0.0;
    static double csw0 = 0.0;
    if (m0 != *mass || csw0 != get_csw()) {
        m0 = *mass;
        csw0 = get_csw();
        int n;
        double a, b, c;

        n = 0;
        c = 3.0 * csw0 / (4.0 + m0);
        a = c * exp(c);
        b = DBL_EPSILON;

        for (n = 1; n < MAX_FACTORIAL; n++) {
            a *= c;
            b *= (double)(n + 1);

            if (a < b) {
                NN = n;
                NNexp = NN + 2;
                lprintf("SWEXP", 0, "SW exp order of the taylor expansion is set to %d\n", n);
                return;
            }
        }

        error(0 == 0, 1, "set_sw_order" __FILE__, "SW parameters are out of range");
    }
}

#endif
