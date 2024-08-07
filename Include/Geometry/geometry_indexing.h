/**
 * @file geometry_indexing.h
 * @brief Finding indices from coordinates and the other way around
 */

#ifndef GEOMETRY_INDEXING_H
#define GEOMETRY_INDEXING_H

#define _lexi(_T, _X, _Y, _Z, _t, _x, _y, _z) ((((_t) * (_X) + (_x)) * (_Y) + (_y)) * (_Z) + (_z))
#define ipt_ext(t, x, y, z) ipt[_lexi(T_EXT, X_EXT, Y_EXT, Z_EXT, t, x, y, z)]
#define ipt_ext_gpu(t, x, y, z) ipt_gpu[_lexi(T_EXT_GPU, X_EXT_GPU, Y_EXT_GPU, Z_EXT_GPU, t, x, y, z)]
#define ipt(t, x, y, z) ipt_ext((t) + T_BORDER, (x) + X_BORDER, (y) + Y_BORDER, (z) + Z_BORDER)
#define imask(ix) imask[ix]

#define ipt_4d(t, x) ipt_4d[(t) * (VOL3) + (x)]
#define iup(site, dir) iup[(site) * 4 + (dir)]
#define idn(site, dir) idn[(site) * 4 + (dir)]

#ifdef WITH_GPU
#define _PTR(_field) (_field)->gpu_ptr
#else
#define _PTR(_field) (_field)->ptr
#endif

#endif
