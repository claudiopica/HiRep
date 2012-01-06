/*Matrix multiplcation and trace routines for four-fermion correlators*/

//z+=Tr(g5 m n)+Tr(g5 n m)
#define _smatrix_trace_g5mul_plus(z,m,n) \
  { \
	for(int _a=0;_a<2*NF;_a++){ \
		for(int _b=0;_b<4*NF;_b++){ \
			_complex_mul_assign(z,(m)[_b].c[_a/NF].c[_a%NF],(n)[_a].c[_b/NF].c[_b%NF]); \
			_complex_mul_assign(z,(n)[_b].c[_a/NF].c[_a%NF],(m)[_a].c[_b/NF].c[_b%NF]); } } \
	for(int _a=2*NF;_a<4*NF;_a++){ \
		for(int _b=0;_b<4*NF;_b++){ \
			_complex_mul_sub_assign(z,(m)[_b].c[_a/NF].c[_a%NF],(n)[_a].c[_b/NF].c[_b%NF]); \
			_complex_mul_sub_assign(z,(n)[_b].c[_a/NF].c[_a%NF],(m)[_a].c[_b/NF].c[_b%NF]); } } \
}((void)0)

//z-=Tr(g5 m n)+Tr(g5 n m)
#define _smatrix_trace_g5mul_minus(z,m,n) \
  { \
	for(int _a=0;_a<2*NF;_a++){ \
		for(int _b=0;_b<4*NF;_b++){ \
			_complex_mul_sub_assign(z,(m)[_b].c[_a/NF].c[_a%NF],(n)[_a].c[_b/NF].c[_b%NF]); \
			_complex_mul_sub_assign(z,(n)[_b].c[_a/NF].c[_a%NF],(m)[_a].c[_b/NF].c[_b%NF]); } } \
	for(int _a=2*NF;_a<4*NF;_a++){ \
		for(int _b=0;_b<4*NF;_b++){ \
			_complex_mul_assign(z,(m)[_b].c[_a/NF].c[_a%NF],(n)[_a].c[_b/NF].c[_b%NF]); \
			_complex_mul_assign(z,(n)[_b].c[_a/NF].c[_a%NF],(m)[_a].c[_b/NF].c[_b%NF]); } } \
}((void)0)

//m = g5 ndag
#define _smatrix_g5ndagp(m,n) \
  { \
	for(int _a=0;_a<2*NF;_a++){ \
		for(int _b=0;_b<4*NF;_b++){ \
			_complex_star((m)[_b].c[_a/NF].c[_a%NF],(*((n)[_a])).c[_b/NF].c[_b%NF]); } } \
	for(int _a=2*NF;_a<4*NF;_a++){ \
		for(int _b=0;_b<4*NF;_b++){ \
			_complex_star_minus((m)[_b].c[_a/NF].c[_a%NF],(*((n)[_a])).c[_b/NF].c[_b%NF]); } } \
}((void)0)

#define _smatrix_m_n(l,m,n) \
  { \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_mul_assign((l)[_a].c[_c/NF].c[_c%NF],(m)[_b].c[_c/NF].c[_c%NF],(n)[_a].c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_mp_n(l,m,n) \
  { \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_mul_assign((l)[_a].c[_c/NF].c[_c%NF],(*((m)[_b])).c[_c/NF].c[_c%NF],(n)[_a].c[_b/NF].c[_b%NF]); } } } \
}((void)0)

  #define _smatrix_m_np(l,m,n) \
  { \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_mul_assign((l)[_a].c[_c/NF].c[_c%NF],(m)[_b].c[_c/NF].c[_c%NF],(*((n)[_a])).c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_mp_np(l,m,n) \
  { \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_mul_assign((l)[_a].c[_c/NF].c[_c%NF],(*((m)[_b])).c[_c/NF].c[_c%NF],(*((n)[_a])).c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_mdag_n(l,m,n) \
  { \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(m)[_c].c[_b/NF].c[_b%NF],(n)[_a].c[_b/NF].c[_b%NF]); } } } \
}((void)0)
  
#define _smatrix_mdagp_n(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(*((m)[_c])).c[_b/NF].c[_b%NF],(n)[_a].c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_mdag_np(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(m)[_c].c[_b/NF].c[_b%NF],(*((n)[_a])).c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_mdagp_np(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(*((m)[_c])).c[_b/NF].c[_b%NF],(*((n)[_a])).c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_m_ndag(l,m,n) \
  { \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(n)[_b].c[_a/NF].c[_a%NF],(m)[_b].c[_c/NF].c[_c%NF]); } } } \
}((void)0)
  
#define _smatrix_mp_ndag(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(n)[_b].c[_a/NF].c[_a%NF],(*((m)[_b])).c[_c/NF].c[_c%NF]); } } } \
}((void)0)

#define _smatrix_m_ndagp(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(*((n)[_b])).c[_a/NF].c[_a%NF],(m)[_b].c[_c/NF].c[_c%NF]); } } } \
}((void)0)

#define _smatrix_mp_ndagp(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_prod_assign((l)[_a].c[_c/NF].c[_c%NF],(*((n)[_b])).c[_a/NF].c[_a%NF],(*((m)[_b])).c[_c/NF].c[_c%NF]); } } } \
}((void)0)

#define _smatrix_mdag_ndagp(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_mul_star_star_assign((l)[_a].c[_c/NF].c[_c%NF],(*((n)[_b])).c[_a/NF].c[_a%NF],(m)[_c].c[_b/NF].c[_b%NF]); } } } \
}((void)0)

#define _smatrix_mdagp_ndag(l,m,n) \
{ \
	for(int _a=0;_a<4*NF;_a++){ \
		_spinor_zero_f((l)[_a]); \
		for(int _b=0;_b<4*NF;_b++){ \
			for(int _c=0;_c<4*NF;_c++){ \
				_complex_mul_star_star_assign((l)[_a].c[_c/NF].c[_c%NF],(n)[_b].c[_a/NF].c[_a%NF],(*((m)[_c])).c[_b/NF].c[_b%NF]); } } } \
}((void)0)

//Returns:
//_disc-= Tr(gmu _A)Tr(g5gmu _B) + Tr(g5gmu _A)Tr(gmu _B)
//_conn-= Tr(g5gmu _A gmu _B) + Tr(gmu _A g5gmu _B)
//assumes existence of ptmp1, ptmp2, z1, z2 as temporary variables
#define _TRACE_minus(_A,_B,_disc,_conn) \
{ \
			_smatrix_g0_f(ptmp1,(_A)); \
			_smatrix_g0_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_minus((_conn),ptmp1,ptmp2); \
			_smatrix_g1_f(ptmp1,(_A)); \
			_smatrix_g1_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_minus((_conn),ptmp1,ptmp2); \
			_smatrix_g2_f(ptmp1,(_A)); \
			_smatrix_g2_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_minus((_conn),ptmp1,ptmp2); \
			_smatrix_g3_f(ptmp1,(_A)); \
			_smatrix_g3_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_minus((_conn),ptmp1,ptmp2); \
}((void)0)

//Returns:
//_disc+= Tr(gmu _A)Tr(g5gmu _B) + Tr(g5gmu _A)Tr(gmu _B)
//_conn+= Tr(g5gmu _A gmu _B) + Tr(gmu _A g5gmu _B)
//assumes existence of ptmp1, ptmp2, z1, z2 as temporary variables
#define _TRACE_plus(_A,_B,_disc,_conn) \
{ \
			_smatrix_g0_f(ptmp1,(_A)); \
			_smatrix_g0_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_plus((_conn),ptmp1,ptmp2); \
			_smatrix_g1_f(ptmp1,(_A)); \
			_smatrix_g1_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_plus((_conn),ptmp1,ptmp2); \
			_smatrix_g2_f(ptmp1,(_A)); \
			_smatrix_g2_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_plus((_conn),ptmp1,ptmp2); \
			_smatrix_g3_f(ptmp1,(_A)); \
			_smatrix_g3_f(ptmp2,(_B)); \
			_smatrix_traceg5_f(z1,ptmp1); \
			_smatrix_trace_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_f(z1,ptmp1); \
			_smatrix_traceg5_f(z2,ptmp2); \
			_complex_mul_sub_assign((_disc),z1,z2); \
			_smatrix_trace_g5mul_plus((_conn),ptmp1,ptmp2); \
}((void)0)
