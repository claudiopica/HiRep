#include "libhr_core.h"

#define XG(m, a, b) ((m) + (a) * NG + (b))
#define XF(m, a, b) ((m) + (a) * NF + (b))

#ifdef __cplusplus
extern "C" {
#endif

visible void _group_represent2(suNf *v, suNg *u) {
#ifdef WITH_QUATERNIONS
    *v = *((suNf *)u);
#elif defined REPR_ADJOINT

    int A, C;
    int a, b, i, j, k, c, d;
    double *vf = (double *)v;
    hr_complex *uf = (hr_complex *)u;

    suNg m;
    hr_complex *mf = (hr_complex *)(&m);

    A = 0;
    for (a = 0; a < NG; a++) {
        for (b = (a == 0) ? 1 : 0; b < NG; b++) {
            if (a > b) {
                for (i = 0; i < NG; i++) {
                    for (j = i; j < NG; j++) {
                        *XG(mf, i, j) = (*XG(uf, i, a)) * conj(*XG(uf, j, b)) + (*XG(uf, i, b)) * conj(*XG(uf, j, a));
                    }
                }
            } else if (a < b) {
                for (i = 0; i < NG; i++) {
                    for (j = i; j < NG; j++) {
                        *XG(mf, i, j) =
                            I * ((*XG(uf, i, a)) * conj(*XG(uf, j, b))) - I * ((*XG(uf, i, b)) * conj(*XG(uf, j, a)));
                    }
                }
            } else if (a == b) {
                for (i = 0; i < NG; i++) {
                    for (j = i; j < NG; j++) {
                        *XG(mf, i, j) = -a * (*XG(uf, i, a)) * conj(*XG(uf, j, a));
                        for (k = 0; k < a; k++) {
                            *XG(mf, i, j) += (*XG(uf, i, k)) * conj(*XG(uf, j, k));
                        }
                        *XG(mf, i, j) *= sqrt(2. / (a * (a + 1.)));
                    }
                }
            }

            C = 0;
            for (c = 0; c < NG; c++) {
                for (d = (c == 0) ? 1 : 0; d < NG; d++) {
                    if (c > d) {
                        *(XF(vf, C, A)) = creal(*XG(mf, d, c));
                    } else if (c < d) {
                        *(XF(vf, C, A)) = cimag(*XG(mf, c, d));
                    } else if (c == d) {
                        *(XF(vf, C, A)) = -c * creal(*XG(mf, c, c));
                        for (k = 0; k < c; k++) {
                            *(XF(vf, C, A)) += creal(*XG(mf, k, k));
                        }
                        *(XF(vf, C, A)) *= sqrt(.5 / (c * (c + 1.)));
                    }

                    C++;
                }
            }

            A++;
        }
    }

#elif defined REPR_SYMMETRIC

    const double st = sqrt(2.);
    int A, C;
    int a, b, i, j, c, d;
    hr_complex *vf = (hr_complex *)v;
    hr_complex *uf = (hr_complex *)u;

    suNg m;
    hr_complex *mf = (hr_complex *)(&m);

    A = 0;
    for (a = 0; a < NG; a++) {
        for (b = 0; b < a; b++) {
            for (i = 0; i < NG; i++) {
                for (j = i; j < NG; j++) {
                    *XG(mf, i, j) = (*XG(uf, i, a)) * (*XG(uf, j, b)) + (*XG(uf, i, b)) * (*XG(uf, j, a));
                }
            }

            C = 0;
            for (c = 0; c < NG; c++) {
                for (d = 0; d < c; d++) {
                    *XF(vf, C, A) = *XG(mf, d, c);
                    C++;
                }
                *XF(vf, C, A) = (*XG(mf, c, c)) / st;
                C++;
            }

            A++;
        }

        for (i = 0; i < NG; i++) {
            for (j = i; j < NG; j++) {
                *XG(mf, i, j) = (*XG(uf, i, a)) * (*XG(uf, j, a));
            }
        }

        C = 0;
        for (c = 0; c < NG; c++) {
            for (d = 0; d < c; d++) {
                *XF(vf, C, A) = (*XG(mf, d, c)) * st;
                C++;
            }
            *XF(vf, C, A) = (*XG(mf, c, c));
            C++;
        }

        A++;
    }

#elif defined REPR_ANTISYMMETRIC

    int A, C;
    int a, b, i, j, c, d;
    hr_complex *vf = (hr_complex *)v;
    hr_complex *uf = (hr_complex *)u;

    suNg m;
    hr_complex *mf = (hr_complex *)(&m);

    A = 0;
    for (a = 1; a < NG; a++) {
        for (b = 0; b < a; b++) {
            for (i = 0; i < NG; i++) {
                for (j = i; j < NG; j++) {
                    *XG(mf, i, j) = (*XG(uf, i, a)) * (*XG(uf, j, b)) - (*XG(uf, i, b)) * (*XG(uf, j, a));
                }
            }

            C = 0;
            for (c = 1; c < NG; c++) {
                for (d = 0; d < c; d++) {
                    *XF(vf, C, A) = -(*XG(mf, d, c));
                    C++;
                }
            }

            A++;
        }
    }

#elif defined REPR_FUNDAMENTAL
    *v = *((suNf *)u);
#endif
}

visible void _group_represent_flt(suNf_flt *v, suNg_flt *u) {
#ifdef REPR_ADJOINT

    int A, C;
    int a, b, i, j, k, c, d;
    float *vf = (float *)v;
    hr_complex_flt *uf = (hr_complex_flt *)u;

    suNg_flt m;
    hr_complex_flt *mf = (hr_complex_flt *)(&m);

    A = 0;
    for (a = 0; a < NG; a++) {
        for (b = (a == 0) ? 1 : 0; b < NG; b++) {
            if (a > b) {
                for (i = 0; i < NG; i++) {
                    for (j = i; j < NG; j++) {
                        *XG(mf, i, j) = (*XG(uf, i, a)) * conj(*XG(uf, j, b)) + (*XG(uf, i, b)) * conj(*XG(uf, j, a));
                    }
                }
            } else if (a < b) {
                for (i = 0; i < NG; i++) {
                    for (j = i; j < NG; j++) {
                        *XG(mf, i, j) =
                            I * ((*XG(uf, i, a)) * conj(*XG(uf, j, b))) - I * ((*XG(uf, i, b)) * conj(*XG(uf, j, a)));
                    }
                }
            } else if (a == b) {
                for (i = 0; i < NG; i++) {
                    for (j = i; j < NG; j++) {
                        *XG(mf, i, j) = -a * (*XG(uf, i, a)) * conj(*XG(uf, j, a));
                        for (k = 0; k < a; k++) {
                            *XG(mf, i, j) += (*XG(uf, i, k)) * conj(*XG(uf, j, k));
                        }
                        *XG(mf, i, j) *= (float)sqrt(2. / (a * (a + 1.)));
                    }
                }
            }

            C = 0;
            for (c = 0; c < NG; c++) {
                for (d = (c == 0) ? 1 : 0; d < NG; d++) {
                    if (c > d) {
                        *(XF(vf, C, A)) = creal(*XG(mf, d, c));
                    } else if (c < d) {
                        *(XF(vf, C, A)) = cimag(*XG(mf, c, d));
                    } else if (c == d) {
                        *(XF(vf, C, A)) = -c * creal(*XG(mf, c, c));
                        for (k = 0; k < c; k++) {
                            *(XF(vf, C, A)) += creal(*XG(mf, k, k));
                        }
                        *(XF(vf, C, A)) *= (float)sqrt(.5 / (c * (c + 1.)));
                    }

                    C++;
                }
            }

            A++;
        }
    }

#elif defined REPR_SYMMETRIC

    const double st = sqrt(2.);
    int A, C;
    int a, b, i, j, c, d;
    hr_complex_flt *vf = (hr_complex_flt *)v;
    hr_complex_flt *uf = (hr_complex_flt *)u;

    suNg_flt m;
    hr_complex_flt *mf = (hr_complex_flt *)(&m);

    A = 0;
    for (a = 0; a < NG; a++) {
        for (b = 0; b < a; b++) {
            for (i = 0; i < NG; i++) {
                for (j = i; j < NG; j++) {
                    *XG(mf, i, j) = (*XG(uf, i, a)) * (*XG(uf, j, b)) + (*XG(uf, i, b)) * (*XG(uf, j, a));
                }
            }

            C = 0;
            for (c = 0; c < NG; c++) {
                for (d = 0; d < c; d++) {
                    *XF(vf, C, A) = *XG(mf, d, c);
                    C++;
                }
                *XF(vf, C, A) = (*XG(mf, c, c)) / st;
                C++;
            }

            A++;
        }

        for (i = 0; i < NG; i++) {
            for (j = i; j < NG; j++) {
                *XG(mf, i, j) = (*XG(uf, i, a)) * (*XG(uf, j, a));
            }
        }

        C = 0;
        for (c = 0; c < NG; c++) {
            for (d = 0; d < c; d++) {
                *XF(vf, C, A) = (*XG(mf, d, c)) * st;
                C++;
            }
            *XF(vf, C, A) = (*XG(mf, c, c));
            C++;
        }

        A++;
    }

#elif defined REPR_ANTISYMMETRIC

    int A, C;
    int a, b, i, j, c, d;
    hr_complex_flt *vf = (hr_complex_flt *)v;
    hr_complex_flt *uf = (hr_complex_flt *)u;

    suNg_flt m;
    hr_complex_flt *mf = (hr_complex_flt *)(&m);

    A = 0;
    for (a = 1; a < NG; a++) {
        for (b = 0; b < a; b++) {
            for (i = 0; i < NG; i++) {
                for (j = i; j < NG; j++) {
                    *XG(mf, i, j) = (*XG(uf, i, a)) * (*XG(uf, j, b)) - (*XG(uf, i, b)) * (*XG(uf, j, a));
                }
            }

            C = 0;
            for (c = 1; c < NG; c++) {
                for (d = 0; d < c; d++) {
                    *XF(vf, C, A) = -(*XG(mf, d, c));
                    C++;
                }
            }

            A++;
        }
    }

#elif defined REPR_FUNDAMENTAL
    *v = *((suNf_flt *)u);
#endif
}

#ifdef __cplusplus
}
#endif

#undef XG
#undef XF