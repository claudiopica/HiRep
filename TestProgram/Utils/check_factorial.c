#include "libhr.h"

#define CHECK_DIFF(errors, abs1, abs2, _order)                                                                         \
    do {                                                                                                               \
        double rel = fabs(abs1 - abs2) / fabs(abs1);                                                                   \
        const char *msg = (rel > 1e-15) ? ++errors, "[FAIL]" : "[ OK ]";                                               \
        lprintf("CHECK FACTORIAL", 2, "%s, order: %d, rel=%.10e abs=%.10e diff=%.10e\n", msg, _order, rel, fabs(abs1), \
                fabs(abs1 - abs2));                                                                                    \
    } while (0)

static int errors = 0;

int main(int argc, char *argv[]) {
    setup_process(&argc, &argv);

    double inv_fact_test = 1.0;
    double inv_fact_imp;

    for (int k = 0; k <= MAX_FACTORIAL; ++k) {
        if (k > 0) { inv_fact_test = inv_fact_test / k; }
        inv_fact_imp = inverse_fact(k);
        CHECK_DIFF(errors, inv_fact_imp, inv_fact_test, k);
    }

    finalize_process();
    return errors;
}