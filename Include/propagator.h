/*******************************************************************************
*
* File propagator.h
*
* Type definitions and macros for propagator
*
* 2013 Rudy Arthur, Ari Hietanen
*
*******************************************************************************/

#ifndef PROPAGATOR_H
#define PROPAGATOR_H
#include "suN_types.h"
#include "suN.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"

typedef struct
{
  suNf_spin_matrix c[NF];
} suNf_propagator;

//r propagator, s spinor, i color index, j spin index
#define _propagator_assign(p, s, i, j)               \
  do                                                 \
  {                                                  \
    int ITMP;                                        \
    for (ITMP = 0; ITMP < NF; ++ITMP)                \
    {                                                \
      (p).c[ITMP].c[0].c[j].c[i] = (s).c[0].c[ITMP]; \
      (p).c[ITMP].c[1].c[j].c[i] = (s).c[1].c[ITMP]; \
      (p).c[ITMP].c[2].c[j].c[i] = (s).c[2].c[ITMP]; \
      (p).c[ITMP].c[3].c[j].c[i] = (s).c[3].c[ITMP]; \
    }                                                \
  } while (0)

//(r).c[i].c[j] = (s)
//_spinmatrix_assign_col(r.c[i], s, j)

//r propagator, s spinmatrix, i color index
#define _propagator_assign_spin_matrix(r, s, i) \
  (r).c[i] = (s)

//r propagator
#define _PROP_AT(r, a, alpha, beta, b) ((r).c[(a)].c[(alpha)].c[(beta)].c[(b)])

#define _PROP_IDX(r, i, j) ((r).c[((i) / 4)].c[((i) % 4)].c[((j) % 4)].c[((j) / 4)])

//  r=0  (r propagator)
#define _propagator_zero(r)           \
  do                                  \
  {                                   \
    int ITMP;                         \
    for (ITMP = 0; ITMP < NF; ++ITMP) \
    {                                 \
      _spinmatrix_zero((r).c[ITMP]);  \
    }                                 \
  } while (0)

//r propagator k result; Tr [ r ]
#define _propagator_one(r)                     \
  _propagator_zero(r);                         \
  do                                           \
  {                                            \
    int ITMP;                                  \
    for (ITMP = 0; ITMP < NF; ITMP++)          \
    {                                          \
      _complex_1(r.c[ITMP].c[0].c[0].c[ITMP]); \
      _complex_1(r.c[ITMP].c[1].c[1].c[ITMP]); \
      _complex_1(r.c[ITMP].c[2].c[2].c[ITMP]); \
      _complex_1(r.c[ITMP].c[3].c[3].c[ITMP]); \
    }                                          \
  } while (0)

//r propagator k result; Tr [ r ]
#define _propagator_trace(k, r)             \
  do                                        \
  {                                         \
    _complex_0(k);                          \
    int ITMP;                               \
    for (ITMP = 0; ITMP < NF; ITMP++)       \
    {                                       \
      (k) += (r).c[ITMP].c[0].c[0].c[ITMP]; \
      (k) += (r).c[ITMP].c[1].c[1].c[ITMP]; \
      (k) += (r).c[ITMP].c[2].c[2].c[ITMP]; \
      (k) += (r).c[ITMP].c[3].c[3].c[ITMP]; \
    }                                       \
  } while (0)

#define _propagator_add(p, q, r)                                             \
  do                                                                         \
  {                                                                          \
    int _a, _beta;                                                           \
    for (_a = 0; _a < NF; ++_a)                                              \
      for (_beta = 0; _beta < 4; ++_beta)                                    \
      {                                                                      \
        _spinor_add_f(p.c[_a].c[_beta], q.c[_a].c[_beta], r.c[_a].c[_beta]); \
      }                                                                      \
  } while (0)

#define _propagator_sub(p, q, r)                                             \
  do                                                                         \
  {                                                                          \
    int _a, _beta;                                                           \
    for (_a = 0; _a < NF; ++_a)                                              \
      for (_beta = 0; _beta < 4; ++_beta)                                    \
      {                                                                      \
        _spinor_sub_f(p.c[_a].c[_beta], q.c[_a].c[_beta], r.c[_a].c[_beta]); \
      }                                                                      \
  } while (0)

//S propagator k factor;
#define _propagator_mul_assign(S, k)                              \
  do                                                              \
  {                                                               \
    int _a, _beta;                                                \
    for (_a = 0; _a < NF; ++_a)                                   \
      for (_beta = 0; _beta < 4; ++_beta)                         \
      {                                                           \
        _spinor_mul_f((S).c[_a].c[_beta], k, (S).c[_a].c[_beta]); \
      }                                                           \
  } while (0)

#define _propagator_mulc_assign(S, k)                   \
  do                                                    \
  {                                                     \
    int ITMP, JTMP;                                     \
    for (ITMP = 0; ITMP < 4 * NF; ITMP++)               \
    {                                                   \
      for (JTMP = 0; JTMP < 4 * NF; JTMP++)             \
      {                                                 \
        double complex tmp;                             \
        tmp = _PROP_IDX(S, ITMP, JTMP);                 \
        _complex_mul(_PROP_IDX(S, ITMP, JTMP), k, tmp); \
      }                                                 \
    }                                                   \
  } while (0)

//r propagator = s^T
#define _propagator_transpose(r, s)                              \
  do                                                             \
  {                                                              \
    int ITMP, JTMP;                                              \
    for (ITMP = 0; ITMP < 4 * NF; ITMP++)                        \
    {                                                            \
      for (JTMP = 0; JTMP < 4 * NF; JTMP++)                      \
      {                                                          \
        _PROP_IDX((r), ITMP, JTMP) = _PROP_IDX((s), JTMP, ITMP); \
      }                                                          \
    }                                                            \
  } while (0)

//r propagator = s^dagger
#define _propagator_dagger(r, s)                                       \
  do                                                                   \
  {                                                                    \
    int ITMP, JTMP;                                                    \
    for (ITMP = 0; ITMP < 4 * NF; ITMP++)                              \
    {                                                                  \
      for (JTMP = 0; JTMP < 4 * NF; JTMP++)                            \
      {                                                                \
        _PROP_IDX((r), ITMP, JTMP) = conj(_PROP_IDX((s), JTMP, ITMP)); \
      }                                                                \
    }                                                                  \
  } while (0)

//s spinor = r propagator t spinor
#define _propagator_mul_spinor(s, r, t)                                                                                       \
  do                                                                                                                          \
  {                                                                                                                           \
    _spinor_zero_f(s);                                                                                                        \
    int _a, _b, _alpha, _beta;                                                                                                \
    for (_alpha = 0; _alpha < NF; ++_alpha)                                                                                   \
    {                                                                                                                         \
      for (_a = 0; _a < NF; ++_a)                                                                                             \
      {                                                                                                                       \
        for (_beta = 0; _beta < NF; ++_beta)                                                                                  \
        {                                                                                                                     \
          for (_b = 0; _b < NF; ++_b)                                                                                         \
          {                                                                                                                   \
            _complex_mul_assign((s).c[(_alpha)].c[(_a)], (r).c[(_a)].c[(_alpha)].c[(_beta)].c[(_b)], (t).c[(_beta)].c[(_b)]); \
          }                                                                                                                   \
        }                                                                                                                     \
      }                                                                                                                       \
    }                                                                                                                         \
  } while (0)

//s spinor =  t^dagger spinor r propagator
#define _propagator_leftmul_spinor(s, t, r)                                                                                        \
  do                                                                                                                               \
  {                                                                                                                                \
    _spinor_zero_f(s);                                                                                                             \
    int _a, _b, _alpha, _beta;                                                                                                     \
    for (_beta = 0; _beta < NF; ++_beta)                                                                                           \
    {                                                                                                                              \
      for (_b = 0; _b < NF; ++_b)                                                                                                  \
      {                                                                                                                            \
        for (_alpha = 0; _alpha < NF; ++_alpha)                                                                                    \
        {                                                                                                                          \
          for (_a = 0; _a < NF; ++_a)                                                                                              \
          {                                                                                                                        \
            _complex_mul_star_assign((s).c[(_beta)].c[(_b)], (t).c[(_alpha)].c[(_a)], (r).c[(_a)].c[(_alpha)].c[(_beta)].c[(_b)]); \
          }                                                                                                                        \
        }                                                                                                                          \
      }                                                                                                                            \
    }                                                                                                                              \
  } while (0)

//Q propagator S propagator R propagator factor; Q = SR
#define _propagator_mul(Q, S, R)                                                                                                            \
  do                                                                                                                                        \
  {                                                                                                                                         \
    int _a, _b, _c, _alpha, _beta, _gamma;                                                                                                  \
    for (_a = 0; _a < NF; ++_a)                                                                                                             \
    {                                                                                                                                       \
      for (_b = 0; _b < NF; ++_b)                                                                                                           \
      {                                                                                                                                     \
        for (_alpha = 0; _alpha < 4; ++_alpha)                                                                                              \
        {                                                                                                                                   \
          for (_beta = 0; _beta < 4; ++_beta)                                                                                               \
          {                                                                                                                                 \
            _complex_0(Q.c[_a].c[_alpha].c[_beta].c[_b]);                                                                                   \
            for (_c = 0; _c < NF; _c++)                                                                                                     \
            {                                                                                                                               \
              for (_gamma = 0; _gamma < 4; _gamma++)                                                                                        \
              {                                                                                                                             \
                _complex_mul_assign(Q.c[_a].c[_alpha].c[_beta].c[_b], S.c[_a].c[_alpha].c[_gamma].c[_c], R.c[_c].c[_gamma].c[_beta].c[_b]); \
              }                                                                                                                             \
            }                                                                                                                               \
          }                                                                                                                                 \
        }                                                                                                                                   \
      }                                                                                                                                     \
    }                                                                                                                                       \
  } while (0)

//tr complex S propagator R propagator factor; tr = Trace[ S R ]
#define _propagator_mul_trace(tr, S, R)                                                                  \
  do                                                                                                     \
  {                                                                                                      \
    int _a, _b, _alpha, _beta;                                                                           \
    _complex_0(tr);                                                                                      \
    for (_a = 0; _a < NF; ++_a)                                                                          \
    {                                                                                                    \
      for (_alpha = 0; _alpha < 4; ++_alpha)                                                             \
      {                                                                                                  \
        for (_b = 0; _b < NF; ++_b)                                                                      \
        {                                                                                                \
          for (_beta = 0; _beta < 4; ++_beta)                                                            \
          {                                                                                              \
            _complex_mul_assign(tr, S.c[_a].c[_alpha].c[_beta].c[_b], R.c[_b].c[_beta].c[_alpha].c[_a]); \
          }                                                                                              \
        }                                                                                                \
      }                                                                                                  \
    }                                                                                                    \
  } while (0)

//tr complex S propagator R propagator factor; tr = Trace[ S^dagger R ]
#define _propagator_muldag_trace(tr, S, R)                                                                \
  do                                                                                                      \
  {                                                                                                       \
    int _a, _b, _alpha, _beta;                                                                            \
    _complex_0(tr);                                                                                       \
    for (_a = 0; _a < NF; ++_a)                                                                           \
    {                                                                                                     \
      for (_alpha = 0; _alpha < 4; ++_alpha)                                                              \
      {                                                                                                   \
        for (_b = 0; _b < NF; ++_b)                                                                       \
        {                                                                                                 \
          for (_beta = 0; _beta < 4; ++_beta)                                                             \
          {                                                                                               \
            _complex_prod_assign(tr, S.c[_b].c[_beta].c[_alpha].c[_a], R.c[_b].c[_beta].c[_alpha].c[_a]); \
          }                                                                                               \
        }                                                                                                 \
      }                                                                                                   \
    }                                                                                                     \
  } while (0)

//Color matrix U x propagator
#define _suNf_prop_multiply(us, u, s)                       \
  do                                                        \
  {                                                         \
    suNf_vector v1, v2;                                     \
    int _a, _b, _alpha, _beta;                              \
    for (_beta = 0; _beta < 4; ++_beta)                     \
      for (_alpha = 0; _alpha < 4; ++_alpha)                \
      {                                                     \
        for (_b = 0; _b < NF; ++_b)                         \
        {                                                   \
          for (_a = 0; _a < NF; ++_a)                       \
          {                                                 \
            v1.c[_a] = (s).c[_a].c[_alpha].c[_beta].c[_b];  \
          }                                                 \
          _suNf_multiply(v2, (u), v1);                      \
          for (_a = 0; _a < NF; ++_a)                       \
          {                                                 \
            (us).c[_a].c[_alpha].c[_beta].c[_b] = v2.c[_a]; \
          }                                                 \
        }                                                   \
      }                                                     \
  } while (0)

//Color matrix U^dag x propagator
#define _suNf_inverse_prop_multiply(us, u, s)               \
  do                                                        \
  {                                                         \
    suNf_vector v1, v2;                                     \
    int _a, _b, _alpha, _beta;                              \
    for (_beta = 0; _beta < 4; ++_beta)                     \
      for (_alpha = 0; _alpha < 4; ++_alpha)                \
      {                                                     \
        for (_b = 0; _b < NF; ++_b)                         \
        {                                                   \
          for (_a = 0; _a < NF; ++_a)                       \
          {                                                 \
            v1.c[_a] = (s).c[_a].c[_alpha].c[_beta].c[_b];  \
          }                                                 \
          _suNf_inverse_multiply(v2, (u), v1);              \
          for (_a = 0; _a < NF; ++_a)                       \
          {                                                 \
            (us).c[_a].c[_alpha].c[_beta].c[_b] = v2.c[_a]; \
          }                                                 \
        }                                                   \
      }                                                     \
  } while (0)

#define _id_propagator(p, q)                                     \
  do                                                             \
  {                                                              \
    int ITMP, JTMP;                                              \
    for (ITMP = 0; ITMP < 4 * NF; ITMP++)                        \
    {                                                            \
      for (JTMP = 0; JTMP < 4 * NF; JTMP++)                      \
      {                                                          \
        _PROP_IDX((p), ITMP, JTMP) = _PROP_IDX((q), ITMP, JTMP); \
      }                                                          \
    }                                                            \
  } while (0)

//P = Gamma Q
#define _g0_propagator(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g0_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _g1_propagator(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g1_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _g2_propagator(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g2_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _g3_propagator(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g3_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _g5_propagator(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g5_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _g5g0_propagator(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _g5g0_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _g5g3_propagator(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _g5g3_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _g5g1_propagator(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _g5g1_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _g5g2_propagator(p, q)            \
  do                                      \
  {                                       \
    int _a;                               \
    for (_a = 0; _a < NF; ++_a)           \
    {                                     \
      _g5g2_spinmatrix(p.c[_a], q.c[_a]); \
    }                                     \
  } while (0)

//P = Q Gamma
#define _propagator_g0(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _spinmatrix_g0((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _propagator_g2(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _spinmatrix_g2((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _propagator_g5(p, q)                \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _spinmatrix_g5((p).c[_a], (q).c[_a]); \
    }                                       \
  } while (0)

#define _propagator_g5g0(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _spinmatrix_g5g0((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _g5g0g1_propagator(p, q)            \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g5g0g1_spinmatrix(p.c[_a], q.c[_a]); \
    }                                       \
  } while (0)

#define _g5g0g2_propagator(p, q)                \
  do                                            \
  {                                             \
    int _a;                                     \
    for (_a = 0; _a < NF; ++_a)                 \
    {                                           \
      _g5g0g2_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                           \
  } while (0)

#define _propagator_g5g0g2(p, q)                \
  do                                            \
  {                                             \
    int _a;                                     \
    for (_a = 0; _a < NF; ++_a)                 \
    {                                           \
      _spinmatrix_g5g0g2((p).c[_a], (q).c[_a]); \
    }                                           \
  } while (0)

#define _g5g0g3_propagator(p, q)            \
  do                                        \
  {                                         \
    int _a;                                 \
    for (_a = 0; _a < NF; ++_a)             \
    {                                       \
      _g5g0g3_spinmatrix(p.c[_a], q.c[_a]); \
    }                                       \
  } while (0)

#define _propagator_g5g3(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _spinmatrix_g5g3((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _propagator_g5g1(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _spinmatrix_g5g1((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _g0g1_propagator(p, q)            \
  do                                      \
  {                                       \
    int _a;                               \
    for (_a = 0; _a < NF; ++_a)           \
    {                                     \
      _g0g1_spinmatrix(p.c[_a], q.c[_a]); \
    }                                     \
  } while (0)

#define _g0g2_propagator(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _g0g2_spinmatrix((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _propagator_g0g2(p, q)                \
  do                                          \
  {                                           \
    int _a;                                   \
    for (_a = 0; _a < NF; ++_a)               \
    {                                         \
      _spinmatrix_g0g2((p).c[_a], (q).c[_a]); \
    }                                         \
  } while (0)

#define _g0g3_propagator(p, q)            \
  do                                      \
  {                                       \
    int _a;                               \
    for (_a = 0; _a < NF; ++_a)           \
    {                                     \
      _g0g3_spinmatrix(p.c[_a], q.c[_a]); \
    }                                     \
  } while (0)

#endif
