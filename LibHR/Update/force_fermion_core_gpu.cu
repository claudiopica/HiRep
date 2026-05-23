
#include "geometry.h"
#include "libhr_core.h"
#include "update.h"
#include "memory.h"
#include "utils.h"
#include "inverters.h"

#ifdef WITH_GPU

#ifdef BC_T_THETA
#define _T_theta_mulc(r)                       \
    _vector_mulc_f(ptmp, eitheta_gpu[0], (r)); \
    (r) = ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r)                       \
    _vector_mulc_f(ptmp, eitheta_gpu[1], (r)); \
    (r) = ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r)                       \
    _vector_mulc_f(ptmp, eitheta_gpu[2], (r)); \
    (r) = ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r)                       \
    _vector_mulc_f(ptmp, eitheta_gpu[3], (r)); \
    (r) = ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_DIR0(u, chi1, chi2, ug)               \
    _vector_add_f(ptmp, chi2.c[0], chi2.c[2]);   \
    _suNf_multiply(p.c[0], ug, ptmp);            \
    _T_theta_mulc(p.c[0]);                       \
    _vector_add_f(ptmp, chi2.c[1], chi2.c[3]);   \
    _suNf_multiply(p.c[1], ug, ptmp);            \
    _T_theta_mulc(p.c[1]);                       \
    _vector_sub_f(p.c[2], chi1.c[0], chi1.c[2]); \
    _vector_sub_f(p.c[3], chi1.c[1], chi1.c[3]); \
    _suNf_FMAT((u), p)

#define _F_DIR1(u, chi1, chi2, ug)                 \
    _vector_i_add_f(ptmp, chi2.c[0], chi2.c[3]);   \
    _suNf_multiply(p.c[0], ug, ptmp);              \
    _X_theta_mulc(p.c[0]);                         \
    _vector_i_add_f(ptmp, chi2.c[1], chi2.c[2]);   \
    _suNf_multiply(p.c[1], ug, ptmp);              \
    _X_theta_mulc(p.c[1]);                         \
    _vector_i_sub_f(p.c[2], chi1.c[0], chi1.c[3]); \
    _vector_i_sub_f(p.c[3], chi1.c[1], chi1.c[2]); \
    _suNf_FMAT((u), p)

#define _F_DIR2(u, chi1, chi2, ug)               \
    _vector_add_f(ptmp, chi2.c[0], chi2.c[3]);   \
    _suNf_multiply(p.c[0], ug, ptmp);            \
    _Y_theta_mulc(p.c[0]);                       \
    _vector_sub_f(ptmp, chi2.c[1], chi2.c[2]);   \
    _suNf_multiply(p.c[1], ug, ptmp);            \
    _Y_theta_mulc(p.c[1]);                       \
    _vector_sub_f(p.c[2], chi1.c[0], chi1.c[3]); \
    _vector_add_f(p.c[3], chi1.c[1], chi1.c[2]); \
    _suNf_FMAT((u), p)

#define _F_DIR3(u, chi1, chi2, ug)                 \
    _vector_i_add_f(ptmp, chi2.c[0], chi2.c[2]);   \
    _suNf_multiply(p.c[0], ug, ptmp);              \
    _Z_theta_mulc(p.c[0]);                         \
    _vector_i_sub_f(ptmp, chi2.c[1], chi2.c[3]);   \
    _suNf_multiply(p.c[1], ug, ptmp);              \
    _Z_theta_mulc(p.c[1]);                         \
    _vector_i_sub_f(p.c[2], chi1.c[0], chi1.c[2]); \
    _vector_i_add_f(p.c[3], chi1.c[1], chi1.c[3]); \
    _suNf_FMAT((u), p)

static suNg_av_field *force_sum = NULL;

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)

__device__ static __forceinline__ void g5_sigma(suNf_spinor *s, suNf_spinor *u, int mu, int nu) {
    if (mu == 0 && nu == 1) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = -I * u->c[1].c[i];
            s->c[1].c[i] = -I * u->c[0].c[i];
            s->c[2].c[i] = -I * u->c[3].c[i];
            s->c[3].c[i] = -I * u->c[2].c[i];
        }
    } else if (nu == 0 && mu == 1) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = I * u->c[1].c[i];
            s->c[1].c[i] = I * u->c[0].c[i];
            s->c[2].c[i] = I * u->c[3].c[i];
            s->c[3].c[i] = I * u->c[2].c[i];
        }
    } else if (mu == 0 && nu == 2) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = -u->c[1].c[i];
            s->c[1].c[i] = u->c[0].c[i];
            s->c[2].c[i] = -u->c[3].c[i];
            s->c[3].c[i] = u->c[2].c[i];
        }
    } else if (nu == 0 && mu == 2) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = u->c[1].c[i];
            s->c[1].c[i] = -u->c[0].c[i];
            s->c[2].c[i] = u->c[3].c[i];
            s->c[3].c[i] = -u->c[2].c[i];
        }
    } else if (mu == 0 && nu == 3) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = -I * u->c[0].c[i];
            s->c[1].c[i] = I * u->c[1].c[i];
            s->c[2].c[i] = -I * u->c[2].c[i];
            s->c[3].c[i] = I * u->c[3].c[i];
        }
    } else if (nu == 0 && mu == 3) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = I * u->c[0].c[i];
            s->c[1].c[i] = -I * u->c[1].c[i];
            s->c[2].c[i] = I * u->c[2].c[i];
            s->c[3].c[i] = -I * u->c[3].c[i];
        }
    } else if (mu == 1 && nu == 2) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = I * u->c[0].c[i];
            s->c[1].c[i] = -I * u->c[1].c[i];
            s->c[2].c[i] = -I * u->c[2].c[i];
            s->c[3].c[i] = I * u->c[3].c[i];
        }
    } else if (nu == 1 && mu == 2) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = -I * u->c[0].c[i];
            s->c[1].c[i] = I * u->c[1].c[i];
            s->c[2].c[i] = I * u->c[2].c[i];
            s->c[3].c[i] = -I * u->c[3].c[i];
        }
    } else if (mu == 1 && nu == 3) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = -u->c[1].c[i];
            s->c[1].c[i] = u->c[0].c[i];
            s->c[2].c[i] = u->c[3].c[i];
            s->c[3].c[i] = -u->c[2].c[i];
        }
    } else if (nu == 1 && mu == 3) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = u->c[1].c[i];
            s->c[1].c[i] = -u->c[0].c[i];
            s->c[2].c[i] = -u->c[3].c[i];
            s->c[3].c[i] = u->c[2].c[i];
        }
    } else if (mu == 2 && nu == 3) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = I * u->c[1].c[i];
            s->c[1].c[i] = I * u->c[0].c[i];
            s->c[2].c[i] = -I * u->c[3].c[i];
            s->c[3].c[i] = -I * u->c[2].c[i];
        }
    } else if (nu == 2 && mu == 3) {
        for (int i = 0; i < NF; i++) {
            s->c[0].c[i] = -I * u->c[1].c[i];
            s->c[1].c[i] = -I * u->c[0].c[i];
            s->c[2].c[i] = I * u->c[3].c[i];
            s->c[3].c[i] = I * u->c[2].c[i];
        }
    }
}

__device__ static __forceinline__ suNf fmat_create(suNf_spinor *a_lhs, suNf_spinor *a_rhs, suNf_spinor *b_lhs,
                                                   suNf_spinor *b_rhs) {
    suNf fmat;
    _suNf_zero(fmat);
    for (int i = 0; i < NF; i++) {
        for (int j = 0; j < NF; j++) {
            for (int k = 0; k < 4; k++) {
#ifdef REPR_IS_REAL
                fmat.c[i * NF + j] +=
                    creal(a_lhs->c[k].c[i] * conj(a_rhs->c[k].c[j]) + b_lhs->c[k].c[i] * conj(b_rhs->c[k].c[j]));
#else
                fmat.c[i * NF + j] += a_lhs->c[k].c[i] * conj(a_rhs->c[k].c[j]) + b_lhs->c[k].c[i] * conj(b_rhs->c[k].c[j]);
#endif
            }
        }
    }
    return fmat;
}

#endif

__global__ static void _force_fermion_core(suNf_spinor *Xs, suNf_spinor *Ys, suNg_algebra_vector *force_sum_gpu, suNf *gauge,
                                           int *iup_gpu, double coeff, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        suNg_algebra_vector f, fs;
        suNf_vector ptmp;
        suNf_spinor p;
        suNf ug;

        suNf_FMAT s1;
        int iy;
        suNf_spinor chi1, chi2;
        iy = iup_gpu[4 * ix + 0];
        _suNf_FMAT_zero(s1);
        read_gpu<double>(0, &chi1, Xs, ix, 0, 1);
        read_gpu<double>(0, &chi2, Ys, iy, 0, 1);
        read_gpu<double>(0, &ug, gauge, ix, 0, 4);
        _F_DIR0(s1, chi1, chi2, ug);
        read_gpu<double>(0, &chi1, Ys, ix, 0, 1);
        read_gpu<double>(0, &chi2, Xs, iy, 0, 1);
        _F_DIR0(s1, chi1, chi2, ug);

        read_gpu<double>(0, &fs, force_sum_gpu, ix, 0, 4);
        _algebra_project_FMAT(f, s1);
        _algebra_vector_mul_add_assign_g(fs, coeff, f);
        write_gpu<double>(0, &fs, force_sum_gpu, ix, 0, 4);

        iy = iup_gpu[4 * ix + 1];
        _suNf_FMAT_zero(s1);
        read_gpu<double>(0, &chi1, Xs, ix, 0, 1);
        read_gpu<double>(0, &chi2, Ys, iy, 0, 1);
        read_gpu<double>(0, &ug, gauge, ix, 1, 4);
        _F_DIR1(s1, chi1, chi2, ug);
        read_gpu<double>(0, &chi1, Ys, ix, 0, 1);
        read_gpu<double>(0, &chi2, Xs, iy, 0, 1);
        _F_DIR1(s1, chi1, chi2, ug);

        read_gpu<double>(0, &fs, force_sum_gpu, ix, 1, 4);
        _algebra_project_FMAT(f, s1);
        _algebra_vector_mul_add_assign_g(fs, coeff, f);
        write_gpu<double>(0, &fs, force_sum_gpu, ix, 1, 4);

        iy = iup_gpu[4 * ix + 2];
        _suNf_FMAT_zero(s1);
        read_gpu<double>(0, &chi1, Xs, ix, 0, 1);
        read_gpu<double>(0, &chi2, Ys, iy, 0, 1);
        read_gpu<double>(0, &ug, gauge, ix, 2, 4);
        _F_DIR2(s1, chi1, chi2, ug);
        read_gpu<double>(0, &chi1, Ys, ix, 0, 1);
        read_gpu<double>(0, &chi2, Xs, iy, 0, 1);
        _F_DIR2(s1, chi1, chi2, ug);

        read_gpu<double>(0, &fs, force_sum_gpu, ix, 2, 4);
        _algebra_project_FMAT(f, s1);
        _algebra_vector_mul_add_assign_g(fs, coeff, f);
        write_gpu<double>(0, &fs, force_sum_gpu, ix, 2, 4);

        iy = iup_gpu[4 * ix + 3];
        _suNf_FMAT_zero(s1);
        read_gpu<double>(0, &chi1, Xs, ix, 0, 1);
        read_gpu<double>(0, &chi2, Ys, iy, 0, 1);
        read_gpu<double>(0, &ug, gauge, ix, 3, 4);
        _F_DIR3(s1, chi1, chi2, ug);
        read_gpu<double>(0, &chi1, Ys, ix, 0, 1);
        read_gpu<double>(0, &chi2, Xs, iy, 0, 1);
        _F_DIR3(s1, chi1, chi2, ug);

        read_gpu<double>(0, &fs, force_sum_gpu, ix, 3, 4);
        _algebra_project_FMAT(f, s1);
        _algebra_vector_mul_add_assign_g(fs, coeff, f);
        write_gpu<double>(0, &fs, force_sum_gpu, ix, 3, 4);
    }
}

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
__global__ static void _force_clover_core(suNf *cl_force, suNg_algebra_vector *force_sum_gpu, suNf *gauge, int *iup_gpu,
                                          int *idn_gpu, double dt, double coeff, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;

        suNf Zl[6], W[9];
        suNf s1, s2, s3, fmat;
        suNf u;
        suNg_algebra_vector f, f2;
        int num, sign;

        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                if (mu == nu) { continue; }

                int o1 = iup_gpu[4 * ix + mu];
                int o2 = iup_gpu[4 * ix + nu];
                int o3 = idn_gpu[4 * ix + nu];
                int o4 = iup_gpu[4 * o3 + mu];
                int o5 = iup_gpu[4 * o2 + mu];

                if (mu < nu) {
                    num = nu * (nu - 1) / 2 + mu;
                    sign = +1;
                } else {
                    num = mu * (mu - 1) / 2 + nu;
                    sign = -1;
                }

                read_gpu<double>(0, &Zl[0], cl_force, ix, num, 6);
                read_gpu<double>(0, &Zl[1], cl_force, o1, num, 6);
                read_gpu<double>(0, &Zl[2], cl_force, o3, num, 6);
                read_gpu<double>(0, &Zl[3], cl_force, o4, num, 6);
                read_gpu<double>(0, &Zl[4], cl_force, o5, num, 6);
                read_gpu<double>(0, &Zl[5], cl_force, o2, num, 6);

                // Construct links
                read_gpu<double>(0, &u, gauge, o3, mu, 4);
                _suNf_dagger(W[0], u);
                read_gpu<double>(0, &W[1], gauge, o3, nu, 4);
                read_gpu<double>(0, &W[2], gauge, o1, nu, 4);

                read_gpu<double>(0, &u, gauge, o2, mu, 4);
                _suNf_dagger(W[3], u);
                read_gpu<double>(0, &u, gauge, ix, nu, 4);
                _suNf_dagger(W[4], u);
                read_gpu<double>(0, &u, gauge, o4, nu, 4);
                _suNf_dagger(W[5], u);

                _suNf_times_suNf(W[6], W[0], W[1]);
                _suNf_times_suNf(W[7], W[2], W[3]);
                _suNf_times_suNf(s1, W[5], W[6]);
                _suNf_times_suNf(W[8], W[7], W[4]);
                _suNf_sub_assign(W[8], s1);

                // Calculate sum of forces
                _suNf_times_suNf(fmat, W[8], Zl[0]);
                _suNf_times_suNf(s1, Zl[1], W[8]);
                _suNf_add_assign(fmat, s1);

                _suNf_times_suNf(s1, W[0], Zl[2]);
                _suNf_times_suNf(s2, s1, W[1]);
                _suNf_times_suNf(s3, Zl[3], W[6]);
                _suNf_add_assign(s2, s3);
                _suNf_times_suNf(s1, W[5], s2);
                _suNf_sub_assign(fmat, s1);
                _suNf_times_suNf(s1, W[2], Zl[4]);
                _suNf_times_suNf(s2, s1, W[3]);
                _suNf_times_suNf(s3, W[7], Zl[5]);
                _suNf_add_assign(s2, s3);
                _suNf_times_suNf(s1, s2, W[4]);
                _suNf_add_assign(fmat, s1);
                read_gpu<double>(0, &u, gauge, ix, mu, 4);

                _suNf_times_suNf(s1, u, fmat);

                // Project on force
                _algebra_project(f, s1);
                read_gpu<double>(0, &f2, force_sum_gpu, ix, mu, 4);
                _algebra_vector_mul_add_assign_g(f2, sign * coeff, f);
                write_gpu<double>(0, &f2, force_sum_gpu, ix, mu, 4);
            }
        }
    }
}

#ifdef WITH_CLOVER
__global__ static void _force_clover_fermion(suNf *cl_force, suNf_spinor *Xs, suNf_spinor *Ys, double residue, int N,
                                             int block_start) {
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += gridDim.x * blockDim.x) {
        suNf_spinor tmp_lhs, tmp_rhs;
        suNf_spinor rhs, lhs;
        suNf fm, fm_tmp;

        // (mu,nu) = (0,1)
        read_gpu<double>(0, &rhs, Xs, ix, 0, 1);
        read_gpu<double>(0, &lhs, Ys, ix, 0, 1);
        read_gpu<double>(0, &fm, cl_force, ix, 0, 6);

        g5_sigma(&tmp_rhs, &rhs, 0, 1);
        g5_sigma(&tmp_lhs, &lhs, 0, 1);
        fm_tmp = fmat_create(&tmp_lhs, &rhs, &tmp_rhs, &lhs);
        _suNf_mul(fm_tmp, residue, fm_tmp);
        _suNf_add_assign(fm, fm_tmp);
        write_gpu<double>(0, &fm, cl_force, ix, 0, 6);

        // (mu,nu) = (0,2)
        read_gpu<double>(0, &fm, cl_force, ix, 1, 6);
        g5_sigma(&tmp_rhs, &rhs, 0, 2);
        g5_sigma(&tmp_lhs, &lhs, 0, 2);
        fm_tmp = fmat_create(&tmp_lhs, &rhs, &tmp_rhs, &lhs);
        _suNf_mul(fm_tmp, residue, fm_tmp);
        _suNf_add_assign(fm, fm_tmp);
        write_gpu<double>(0, &fm, cl_force, ix, 1, 6);

        // (mu,nu) = (1,2)
        read_gpu<double>(0, &fm, cl_force, ix, 2, 6);
        g5_sigma(&tmp_rhs, &rhs, 1, 2);
        g5_sigma(&tmp_lhs, &lhs, 1, 2);
        fm_tmp = fmat_create(&tmp_lhs, &rhs, &tmp_rhs, &lhs);
        _suNf_mul(fm_tmp, residue, fm_tmp);
        _suNf_add_assign(fm, fm_tmp);
        write_gpu<double>(0, &fm, cl_force, ix, 2, 6);

        // (mu,nu) = (0,3)
        read_gpu<double>(0, &fm, cl_force, ix, 3, 6);
        g5_sigma(&tmp_rhs, &rhs, 0, 3);
        g5_sigma(&tmp_lhs, &lhs, 0, 3);
        fm_tmp = fmat_create(&tmp_lhs, &rhs, &tmp_rhs, &lhs);
        _suNf_mul(fm_tmp, residue, fm_tmp);
        _suNf_add_assign(fm, fm_tmp);
        write_gpu<double>(0, &fm, cl_force, ix, 3, 6);

        // (mu,nu) = (1,3)
        read_gpu<double>(0, &fm, cl_force, ix, 4, 6);
        g5_sigma(&tmp_rhs, &rhs, 1, 3);
        g5_sigma(&tmp_lhs, &lhs, 1, 3);
        fm_tmp = fmat_create(&tmp_lhs, &rhs, &tmp_rhs, &lhs);
        _suNf_mul(fm_tmp, residue, fm_tmp);
        _suNf_add_assign(fm, fm_tmp);
        write_gpu<double>(0, &fm, cl_force, ix, 4, 6);

        // (mu,nu) = (2,3)
        read_gpu<double>(0, &fm, cl_force, ix, 5, 6);
        g5_sigma(&tmp_rhs, &rhs, 2, 3);
        g5_sigma(&tmp_lhs, &lhs, 2, 3);
        fm_tmp = fmat_create(&tmp_lhs, &rhs, &tmp_rhs, &lhs);
        _suNf_mul(fm_tmp, residue, fm_tmp);
        _suNf_add_assign(fm, fm_tmp);
        write_gpu<double>(0, &fm, cl_force, ix, 5, 6);
    }
}

#endif

#ifdef WITH_EXPCLOVER

visible static __forceinline__ void A_times_spinor(suNf_spinor *out, suNfc *Aplus, suNfc *Aminus, suNf_spinor *in) {
    suNf_vector aux;

    // Comp 0 1
    _suNfc_multiply(out->c[0], Aplus[0], in->c[0]);
    _suNfc_multiply(aux, Aplus[1], in->c[1]);
    _vector_add_assign_f(out->c[0], aux);

    _suNfc_multiply(out->c[1], Aplus[2], in->c[0]);
    _suNfc_multiply(aux, Aplus[3], in->c[1]);
    _vector_add_assign_f(out->c[1], aux);
    // Comp 2 3
    _suNfc_multiply(out->c[2], Aminus[0], in->c[2]);
    _suNfc_multiply(aux, Aminus[1], in->c[3]);
    _vector_add_assign_f(out->c[2], aux);

    _suNfc_multiply(out->c[3], Aminus[2], in->c[2]);
    _suNfc_multiply(aux, Aminus[3], in->c[3]);
    _vector_add_assign_f(out->c[3], aux);
}

// EXP CSW FORCE TERM
__global__ void _force_clover_fermion(double invexpmass, suNf *cl_force, suNfc *cl_term, suNf_spinor *Xs, suNf_spinor *Ys,
                                      double residue, int NNexp, int N, int block_start) {
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += gridDim.x * blockDim.x) {
        // Construct force matrices
        suNf_spinor tmp_lhs, tmp_rhs;
        suNf_spinor lhs_k, rhs_k;
        suNf_spinor xi[2 * NF];
        suNf_vector v1, v2, v3, v4;

        suNf_spinor rhs, lhs;
        suNf fm, fm_tmp;

        suNfc Aplus[4];
        suNfc Aminus[4];
        suNfc s0, s1, s2, s3;

        double Cplus[2 * NF * 2 * NF];
        double Cminus[2 * NF * 2 * NF];

        int k = 0, i = 0;

        // Create matrix Aplus, Aminus
        read_gpu<double>(0, &s0, cl_term, ix, 0, 4);
        read_gpu<double>(0, &s1, cl_term, ix, 1, 4);
        read_gpu<double>(0, &s2, cl_term, ix, 2, 4);
        read_gpu<double>(0, &s3, cl_term, ix, 3, 4);

        _suNf_mul(Aplus[0], invexpmass, s0);
        _suNf_mul(Aplus[1], invexpmass, s1);
        _suNf_dagger(Aplus[2], Aplus[1]);
        _suNf_mul(Aplus[3], -invexpmass, s0);

        _suNf_mul(Aminus[0], invexpmass, s2);
        _suNf_mul(Aminus[1], invexpmass, s3);
        _suNf_dagger(Aminus[2], Aminus[1]);
        _suNf_mul(Aminus[3], -invexpmass, s2);

        // double horner scheme
        doublehorner(Cplus, Aplus, NNexp);
        doublehorner(Cminus, Aminus, NNexp);

        // Remember rhs = eta, lhs  = xi
        read_gpu<double>(0, &rhs, Xs, ix, 0, 1);
        read_gpu<double>(0, &lhs, Ys, ix, 0, 1);

        xi[0] = lhs;
        for (k = 0; k < 2 * NF - 1; k++) {
            A_times_spinor(&xi[k + 1], Aplus, Aminus, &xi[k]);
        }

        for (k = 0; k < 2 * NF; k++) {
            // Calculate eta_k = (Aplus^k*etaplus, Aminus^k*etaminus)
            if (k == 0) {
                rhs_k = rhs;
            } else {
                // Comp 0 1
                _suNfc_multiply(v1, Aplus[0], rhs_k.c[0]);
                _suNfc_multiply(v2, Aplus[1], rhs_k.c[1]);
                _suNfc_multiply(v3, Aplus[2], rhs_k.c[0]);
                _suNfc_multiply(v4, Aplus[3], rhs_k.c[1]);
                _vector_add_f(rhs_k.c[0], v1, v2);
                _vector_add_f(rhs_k.c[1], v3, v4);
                // Comp 2 3
                _suNfc_multiply(v1, Aminus[0], rhs_k.c[2]);
                _suNfc_multiply(v2, Aminus[1], rhs_k.c[3]);
                _suNfc_multiply(v3, Aminus[2], rhs_k.c[2]);
                _suNfc_multiply(v4, Aminus[3], rhs_k.c[3]);
                _vector_add_f(rhs_k.c[2], v1, v2);
                _vector_add_f(rhs_k.c[3], v3, v4);
            }

            // Calculate xi_k = sum_{l} C_{kl} A^l xi
            _spinor_zero_f(lhs_k);
            for (i = 0; i < 2 * NF; i++) {
                _vector_mulc_add_assign_f(lhs_k.c[0], Cplus[k * 2 * NF + i], xi[i].c[0]);
                _vector_mulc_add_assign_f(lhs_k.c[1], Cplus[k * 2 * NF + i], xi[i].c[1]);
                _vector_mulc_add_assign_f(lhs_k.c[2], Cminus[k * 2 * NF + i], xi[i].c[2]);
                _vector_mulc_add_assign_f(lhs_k.c[3], Cminus[k * 2 * NF + i], xi[i].c[3]);
            }

            // (mu,nu) = (0,1)
            read_gpu<double>(0, &fm, cl_force, ix, 0, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 0, 1);
            g5_sigma(&tmp_lhs, &lhs_k, 0, 1);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 0, 6);

            // (mu,nu) = (0,2)
            read_gpu<double>(0, &fm, cl_force, ix, 1, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 0, 2);
            g5_sigma(&tmp_lhs, &lhs_k, 0, 2);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 1, 6);

            // (mu,nu) = (1,2)
            read_gpu<double>(0, &fm, cl_force, ix, 2, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 1, 2);
            g5_sigma(&tmp_lhs, &lhs_k, 1, 2);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 2, 6);

            // (mu,nu) = (0,3)
            read_gpu<double>(0, &fm, cl_force, ix, 3, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 0, 3);
            g5_sigma(&tmp_lhs, &lhs_k, 0, 3);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 3, 6);

            // (mu,nu) = (1,3)
            read_gpu<double>(0, &fm, cl_force, ix, 4, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 1, 3);
            g5_sigma(&tmp_lhs, &lhs_k, 1, 3);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 4, 6);

            // (mu,nu) = (2,3)
            read_gpu<double>(0, &fm, cl_force, ix, 5, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 2, 3);
            g5_sigma(&tmp_lhs, &lhs_k, 2, 3);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 5, 6);
        }
    }
}

__global__ static void _force_clover_fermion_taylor(double invexpmass, suNf *cl_force, suNfc *cl_term, suNf_spinor *Xs,
                                                    suNf_spinor *Ys, double residue, int NNexp, int N, int block_start) {
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += gridDim.x * blockDim.x) {
        double *Coef = (double *)malloc(NNexp * NNexp * sizeof(double));

        int i, k;

        suNf_spinor tmp_lhs, tmp_rhs;
        suNf_spinor lhs_k, rhs_k;
        suNf_spinor *xi = (suNf_spinor *)malloc(NNexp * sizeof(suNf_spinor));
        suNf_vector v1, v2, v3, v4;

        suNf_spinor rhs, lhs;
        suNf fm, fm_tmp;

        suNfc Aplus[4];
        suNfc Aminus[4];
        suNfc s0, s1, s2, s3;
        factorialCoef(Coef, NNexp);

        // Construct force matrices
        // Create matrix Aplus, Aminus
        read_gpu<double>(0, &s0, cl_term, ix, 0, 4);
        read_gpu<double>(0, &s1, cl_term, ix, 1, 4);
        read_gpu<double>(0, &s2, cl_term, ix, 2, 4);
        read_gpu<double>(0, &s3, cl_term, ix, 3, 4);

        _suNf_mul(Aplus[0], invexpmass, s0);
        _suNf_mul(Aplus[1], invexpmass, s1);
        _suNf_dagger(Aplus[2], Aplus[1]);
        _suNf_mul(Aplus[3], -invexpmass, s0);

        _suNf_mul(Aminus[0], invexpmass, s2);
        _suNf_mul(Aminus[1], invexpmass, s3);
        _suNf_dagger(Aminus[2], Aminus[1]);
        _suNf_mul(Aminus[3], -invexpmass, s2);

        // Remember rhs = eta, lhs  = xi
        read_gpu<double>(0, &rhs, Xs, ix, 0, 1);
        read_gpu<double>(0, &lhs, Ys, ix, 0, 1);

        xi[0] = lhs;
        for (k = 0; k < NNexp - 1; k++) {
            A_times_spinor(&xi[k + 1], Aplus, Aminus, &xi[k]);
        }

        for (k = 0; k < NNexp; k++) {
            // Calculate eta_k = (Aplus^k*etaplus, Aminus^k*etaminus)
            if (k == 0) {
                rhs_k = rhs;
            } else {
                // Comp 0 1
                _suNfc_multiply(v1, Aplus[0], rhs_k.c[0]);
                _suNfc_multiply(v2, Aplus[1], rhs_k.c[1]);
                _suNfc_multiply(v3, Aplus[2], rhs_k.c[0]);
                _suNfc_multiply(v4, Aplus[3], rhs_k.c[1]);
                _vector_add_f(rhs_k.c[0], v1, v2);
                _vector_add_f(rhs_k.c[1], v3, v4);
                // Comp 2 3
                _suNfc_multiply(v1, Aminus[0], rhs_k.c[2]);
                _suNfc_multiply(v2, Aminus[1], rhs_k.c[3]);
                _suNfc_multiply(v3, Aminus[2], rhs_k.c[2]);
                _suNfc_multiply(v4, Aminus[3], rhs_k.c[3]);
                _vector_add_f(rhs_k.c[2], v1, v2);
                _vector_add_f(rhs_k.c[3], v3, v4);
            }

            // Calculate xi_k = sum_{l} C_{kl} A^l xi
            _spinor_zero_f(lhs_k);
            for (i = 0; i < NNexp; i++) {
                _vector_mulc_add_assign_f(lhs_k.c[0], Coef[k * NNexp + i], xi[i].c[0]);
                _vector_mulc_add_assign_f(lhs_k.c[1], Coef[k * NNexp + i], xi[i].c[1]);
                _vector_mulc_add_assign_f(lhs_k.c[2], Coef[k * NNexp + i], xi[i].c[2]);
                _vector_mulc_add_assign_f(lhs_k.c[3], Coef[k * NNexp + i], xi[i].c[3]);
            }

            // (mu,nu) = (0,1)
            read_gpu<double>(0, &fm, cl_force, ix, 0, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 0, 1);
            g5_sigma(&tmp_lhs, &lhs_k, 0, 1);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 0, 6);

            // (mu,nu) = (0,2)
            read_gpu<double>(0, &fm, cl_force, ix, 1, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 0, 2);
            g5_sigma(&tmp_lhs, &lhs_k, 0, 2);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 1, 6);

            // (mu,nu) = (1,2)
            read_gpu<double>(0, &fm, cl_force, ix, 2, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 1, 2);
            g5_sigma(&tmp_lhs, &lhs_k, 1, 2);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 2, 6);

            // (mu,nu) = (0,3)
            read_gpu<double>(0, &fm, cl_force, ix, 3, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 0, 3);
            g5_sigma(&tmp_lhs, &lhs_k, 0, 3);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 3, 6);

            // (mu,nu) = (1,3)
            read_gpu<double>(0, &fm, cl_force, ix, 4, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 1, 3);
            g5_sigma(&tmp_lhs, &lhs_k, 1, 3);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 4, 6);

            // (mu,nu) = (2,3)
            read_gpu<double>(0, &fm, cl_force, ix, 5, 6);
            g5_sigma(&tmp_rhs, &rhs_k, 2, 3);
            g5_sigma(&tmp_lhs, &lhs_k, 2, 3);
            fm_tmp = fmat_create(&tmp_lhs, &rhs_k, &tmp_rhs, &lhs_k);
            _suNf_mul(fm_tmp, residue, fm_tmp);
            _suNf_add_assign(fm, fm_tmp);
            write_gpu<double>(0, &fm, cl_force, ix, 5, 6);
        }
        free(Coef);
    }
}

#endif

void force_clover_fermion_gpu(spinor_field *Xs, spinor_field *Ys, double residue) {
#ifdef WITH_EXPCLOVER
    compute_clover_term();
    double invexpmass = get_dirac_mass();
    evaluate_sw_order(&invexpmass);
    invexpmass = 1.0 / (4.0 + invexpmass);

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        suNf *cl_force_gpu = cl_force->gpu_ptr + 6 * block_start;
        suNfc *cl_term_gpu = cl_term->gpu_ptr + 4 * block_start;
        _force_clover_fermion<<<grid, BLOCK_SIZE, 0, 0>>>(invexpmass, cl_force_gpu, cl_term_gpu, _GPU_FIELD_BLK(Xs, ixp),
                                                          _GPU_FIELD_BLK(Ys, ixp), residue, get_NNexp(), N, block_start);
        CudaCheckError();
    }
#endif

#ifdef WITH_CLOVER
    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        suNf *cl_force_gpu = cl_force->gpu_ptr + 6 * block_start;
        _force_clover_fermion<<<grid, BLOCK_SIZE, 0, 0>>>(cl_force_gpu, _GPU_FIELD_BLK(Xs, ixp), _GPU_FIELD_BLK(Ys, ixp),
                                                          residue, N, block_start);
        CudaCheckError();
    }
#endif
}

#ifdef WITH_EXPCLOVER
void force_clover_fermion_taylor_gpu(spinor_field *Xs, spinor_field *Ys, double residue) {
    compute_clover_term();
    double invexpmass = get_dirac_mass();
    evaluate_sw_order(&invexpmass);
    invexpmass = 1.0 / (4.0 + invexpmass);
    int NNexp = get_NNexp();

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        suNf *cl_force_gpu = cl_force->gpu_ptr + 6 * block_start;
        suNfc *cl_term_gpu = cl_term->gpu_ptr + 4 * block_start;
        _force_clover_fermion_taylor<<<grid, BLOCK_SIZE, 0, 0>>>(invexpmass, cl_force_gpu, cl_term_gpu, _GPU_FIELD_BLK(Xs, ixp),
                                                                 _GPU_FIELD_BLK(Ys, ixp), residue, get_NNexp(), N, block_start);
        CudaCheckError();
    }
}
#endif

void force_clover_core_gpu(double dt) {
    double coeff = dt * (_REPR_NORM2 / _FUND_NORM2) * (1. / 8.) * get_csw();

    start_sendrecv_clover_force(cl_force);
    complete_sendrecv_clover_force(cl_force);

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        _force_clover_core<<<grid, BLOCK_SIZE, 0, 0>>>(cl_force->gpu_ptr, force_sum->gpu_ptr, u_gauge_f->gpu_ptr, iup_gpu,
                                                       idn_gpu, dt, coeff, N, block_start);
        CudaCheckError();
    }
}

#endif

void force_fermion_core_gpu(spinor_field *Xs, spinor_field *Ys, int auto_fill_odd, double dt, double residue) {
    double coeff;
    spinor_field Xtmp, Ytmp;

    coeff = residue * dt * (_REPR_NORM2 / _FUND_NORM2);
    Xtmp = *Xs;
    Ytmp = *Ys;
    Xs->type = &glattice;
    Ys->type = &glattice;

#ifdef UPDATE_EO

    if (auto_fill_odd) {
        spinor_field Xe, Xo, Ye, Yo;

        Xe = *Xs;
        Xe.type = &glat_even;
        Xo = *Xs;
        Xo.type = &glat_odd;

        Ye = *Ys;
        Ye.type = &glat_even;
        Yo = *Ys;
        Yo.type = &glat_odd;

        Xo.gpu_ptr = Xs->gpu_ptr + glat_odd.master_shift;
        Yo.gpu_ptr = Ys->gpu_ptr + glat_odd.master_shift;

        Dphi_(&Xo, &Xe);
        Dphi_(&Yo, &Ye);
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
        Cphi_diag_inv(get_dirac_mass(), &Xo, &Xo);
        Cphi_diag_inv(get_dirac_mass(), &Yo, &Yo);
#endif
    }

    coeff = -coeff;
#endif

    // TODO: comms while calculating clover force
    start_sendrecv_spinor_field(Xs);
    complete_sendrecv_spinor_field(Xs);
    start_sendrecv_spinor_field(Ys);
    complete_sendrecv_spinor_field(Ys);

#if defined(WITH_CLOVER)
    force_clover_fermion_gpu(Xs, Ys, residue);
#endif
#if defined(WITH_EXPCLOVER)
#if (NF == 3 || NF == 2)
    //	force_clover_fermion_taylor(Xs, Ys, residue);
    force_clover_fermion_gpu(Xs, Ys, residue);
#else
    force_clover_fermion_taylor_gpu(Xs, Ys, residue);
#endif
#endif //EXP_CLOVER

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        _force_fermion_core<<<grid, BLOCK_SIZE, 0, 0>>>(Xs->gpu_ptr, Ys->gpu_ptr, force_sum->gpu_ptr, u_gauge_f->gpu_ptr,
                                                        iup_gpu, coeff, N, block_start);
        CudaCheckError();
    }

    Xs->type = Xtmp.type;
    Ys->type = Ytmp.type;
}

void fermion_force_begin_gpu() {
    if (force_sum == NULL) { force_sum = alloc_suNg_av_field(&glattice); }

    cudaMemset(force_sum->gpu_ptr, 0, force_sum->type->gsize_gauge * 4 * sizeof(suNg_algebra_vector));
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    cudaMemset(cl_force->gpu_ptr, 0, cl_force->type->gsize_gauge * 6 * sizeof(suNf));
#endif
}

void fermion_force_end_gpu(double dt, suNg_av_field *force) {
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    force_clover_core_gpu(dt);
#endif

#ifdef WITH_SMEARING
    error(1, 1, __func__, "This function is not ported to GPU for WITH_SMEARING");
#else
    add_assign(force, force_sum);
#endif
    apply_BCs_on_momentum_field(force);
}

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
void (*force_clover_fermion)(spinor_field *Xs, spinor_field *Ys, double residue) = force_clover_fermion_gpu;
#endif

void (*force_fermion_core)(spinor_field *Xs, spinor_field *Ys, int auto_fill_odd, double dt,
                           double residue) = force_fermion_core_gpu;
void (*fermion_force_begin)(void) = fermion_force_begin_gpu;
void (*fermion_force_end)(double dt, suNg_av_field *) = fermion_force_end_gpu;

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3

#undef _T_theta_mulc
#undef _X_theta_mulc
#undef _Y_theta_mulc
#undef _Z_theta_mulc

#endif