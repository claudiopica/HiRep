/*This is an automatically generated function, do not edit.*/
#include "libhr.h"
#define Complex(a, b) ((a) + I * (b))
int fullgbcheck(int rotid, hr_complex *rotated, hr_complex *unrotated) {
#define rotfun(a) rotated[(a)]
#define unrotfun(a) unrotated[(a)]
#if total_n_glue_op > 0
    hr_complex tmp[3];
#endif
    int return_value = 0;
    if (0 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(0);
        tmp[1] = 1. * unrotfun(0);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=4 Pre[ax, ay, -ax, -ay][0, 0, 0] + 4 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(1);
        tmp[1] = 1. * unrotfun(1);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                -1, 0, 0, 1, 1, 1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(2);
        tmp[1] = 1. * unrotfun(2);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, -1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(3);
        tmp[1] = 1. * unrotfun(3);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, -1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(4);
        tmp[1] = 1. * unrotfun(4);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_A1plusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = -1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(5);
        tmp[1] = 1. * unrotfun(5);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Pre[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pre[ay, ay, ax, -ay, -ay, -ax][0, 0, 0] + 8 Pre[ay, ay, az, -ay, -ay, -az][0, 0, 0] + 8 Pre[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pre[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) - 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. - 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) + 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = -0.5 * unrotfun(6) - 0.8660254037844386 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = -0.8660254037844386 * unrotfun(6) + 0.5 * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(6);
        tmp[1] = 0. + 1. * unrotfun(6);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_1_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(7);
        tmp[1] = 0. + 1. * unrotfun(7);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=000_EplusOhP+_2_001_A1Dic4+_xy-x-y_001_A1Dic4+_xy-x-y\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 3, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (1 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (2 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (3 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (4 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (5 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (6 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (7 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (8 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (9 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (10 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (11 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (12 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (13 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (14 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (24 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = -0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (31 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (32 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (33 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., 0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (34 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (35 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (36 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] =
            Complex(0., 0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (37 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] =
            Complex(0., -0.5) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (38 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = Complex(0., -0.5) * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + Complex(0., 0.5) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (39 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (40 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - Complex(0., 1.) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 0.7071067811865475 * unrotfun(8) - 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (41 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) + 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 0.7071067811865475 * unrotfun(8) + 0.7071067811865475 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (42 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0.5 * unrotfun(8) - 0.7071067811865475 * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (43 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) + Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) - 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (44 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = -0.5 * unrotfun(8) - Complex(0., 0.7071067811865475) * unrotfun(9) + 0.5 * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. + 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. + 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. + 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(8);
        tmp[1] = 0. - 1. * unrotfun(10);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] - 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] + 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(9);
        tmp[1] = 0. - 1. * unrotfun(9);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Sqrt[2] Pim[ax, ax, ay, -ax, -ax, -ay][0, 0, 0] - 8 Sqrt[2] Pim[ay, ay, ax, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 2, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(10);
        tmp[1] = 0. - 1. * unrotfun(8);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=(-8 I) Pim[ax, ax, az, -ax, -ax, -az][0, 0, 0] + 8 Pim[ay, ay, az, -ay, -ay, -az][0, 0, 0] + (8 I) Pim[az, az, ax, -az, -az, -ax][0, 0, 0] - 8 Pim[az, az, ay, -az, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 0, 4, 3, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (15 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (16 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (21 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (25 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (26 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(11);
        tmp[1] = 1. * unrotfun(11);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[ax, ay, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 0, 1, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(12);
        tmp[1] = 1. * unrotfun(12);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ax, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = 1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (17 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = 1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (18 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = 1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (22 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = 1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (29 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = -1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (30 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = -1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (47 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = -1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(13);
        tmp[1] = -1. * unrotfun(13);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=Pre[ax, ax, -ay, -ax, az, -ax, ay, -az][0, 0, 0] + Pre[ax, ay, az, -ax, -ax, -ay, ax, -az][0, -1, 0] - Pre[ax, az, ax, ay, -az, -ax, -ax, -ay][0, -1, 0] + Pre[ay, ax, az, az, -ay, -az, -ax, -az][0, -1, 0] - Pre[ay, az, ax, ax, -ay, -ax, -az, -ax][0, -1, 0] + Pre[az, ax, az, ay, -ax, -az, -az, -ay][0, -1, 0] - Pre[az, ay, ax, -az, -az, -ay, az, -ax][0, -1, 0] - Pre[az, az, -ay, -az, ax, -az, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                0, 1, 0, 2, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (0 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (20 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (19 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (23 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (27 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (28 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (45 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
    if (46 == rotid) {
        tmp[0] = rotfun(14);
        tmp[1] = 1. * unrotfun(14);
        _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
        if (sqrt(creal(tmp[2])) >= 1.e-10) {
            lprintf(
                "Error", 0,
                " Op=8 Pre[az, ay, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",
                1, 0, 0, 1, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
            return_value++;
        }
    }
#undef unrotfunreturn
#undef rotfun
#undef Complex
    return return_value;
}
