

/*template<typename SITE_TYPE, typename REAL>
__device__ void read_gpu(int stride, SITE_TYPE *s, SITE_TYPE *in, int ix, int dim) {
    #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * 4*NF  + ix % THREADSIZE;
    #else
        int iz = ix;
    #endif
    REAL* in_cpx = (REAL*)in;
    for (int s_comp = 0; s_comp < 4; ++s_comp) {
        for (int vec_comp = 0; vec_comp < NF; ++vec_comp) {
            (*s).c[s_comp].c[vec_comp] = in_cpx[iz];
            iz+=THREADSIZE; 
        }
    }
}

template <typename SITE_TYPE, typename REAL>
__device__ void write_gpu(int stride, SITE_TYPE s, SITE_TYPE *out, int ix) {
   #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * 4*NF  + ix % THREADSIZE;
    #else
        int iz = ix;
    #endif
    REAL * out_cpx = (REAL*)out;
    for (int s_comp = 0; s_comp < 4; ++s_comp) {
      for (int vec_comp = 0; vec_comp < NF; ++vec_comp) {
         out_cpx[iz] = s.c[s_comp].c[vec_comp];
         iz+=THREADSIZE;
      }
    }
}*/
