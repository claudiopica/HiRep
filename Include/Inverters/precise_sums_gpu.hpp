
#ifndef PRECISE_SUMS_GPU_H
#define PRECISE_SUMS_GPU_H

struct __align__(16) quad_double {
    double val;
    double compensation;

    constexpr visible quad_double(const double x = 0.0, const double comp = 0.0)
        : val(x)
        , compensation(comp) {
    }

    visible inline __attribute__((always_inline)) quad_double &operator=(const quad_double &x) {
        val = x.val;
        compensation = x.compensation;
        return *this;
    }

    // Maybe we need to add the compensation term of the added variable?
    visible inline __attribute__((always_inline)) quad_double operator+(const quad_double &x) {
        double new_val = val;
        double new_compensation = compensation;
        _2Sum_acc_hirep(new_val, new_compensation, x.val);
        new_compensation += x.compensation;
        return quad_double(new_val, new_compensation);
    }

    visible inline __attribute__((always_inline)) quad_double &operator+=(const quad_double &x) {
        _2Sum_acc_hirep(val, compensation, x.val);
        compensation += x.compensation;
        return *this;
    }

    visible inline __attribute__((always_inline)) void correct() {
        val += compensation;
        compensation = 0;
    }
};

#endif