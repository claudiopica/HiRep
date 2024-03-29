/******************************************************************************
 *
 *
 * File check_scattering_length.c
 * Checks of the pi pi scattering length calculations 
 * Author: Vincent Drach & Fernando Romero Lopez
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "scattering.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "clover_tools.h"
#include "setup.h"
#include "cinfo.c"
#include "clover_tools.h"


#define PI 3.1415926535
#define SQR(A) ((A) * (A))
#define CMUL(a,b) (a).re*=b;(a).im*=b



/// \cond
#define BASENAME(filename) (strrchr((filename),'/') ? strrchr((filename),'/')+1 : filename )

#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*GLB_T*(nm)+(py)*(n_mom)*GLB_T*(nm)+(pz)*GLB_T*(nm)+ ((cm)*GLB_T) +(tc))
/// \endcond



/// \cond
#define INDEX(px,py,pz,n_mom,tc) ((px + n_mom)*(2*n_mom+1)*(2*n_mom+1)*(GLB_T)+(py + n_mom)*(2*n_mom+1)*(GLB_T)+(pz + n_mom)*(GLB_T)+ (tc))
/// \endcond

/**
 * @brief Structure containing data from the input file relevant to scattering.
 */
typedef struct _input_scatt {
	char mstring[256];
    double csw;
	double precision;
	int nhits;
	int tsrc;
	char outdir[256], bc[16], p[256],configlist[256];

	/* for the reading function */
	input_record_t read[12];

} input_scatt;

#define init_input_scatt(varname) \
{ \
	.read={\
		{"Fermion masses", "mes:mass = %s", STRING_T, (varname).mstring},\
		{"csw", "mes:csw = %lf", DOUBLE_T, &(varname).csw},\
		{"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
		{"number of inversions per cnfg", "sigma_triangle:nhits = %d", INT_T, &(varname).nhits},\
		{"Source time:", "mes:tsrc = %d", INT_T, &(varname).tsrc},\
		{"Output directory:", "mes:outdir = %s", STRING_T, &(varname).outdir},\
		{"Configuration list:", "mes:configlist = %s", STRING_T, &(varname).configlist},\
		{"Boundary conditions:", "mes:bc = %s", STRING_T, &(varname).bc},\
		{NULL, NULL, INT_T, NULL}\
	}\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char prop_filename[256]="";
char source_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "meson_scattering.out";
int Nsource;
double M;

enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };





input_scatt mes_var = init_input_scatt(mes_var);

typedef struct {
	char string[256];
    char configlist[256];
	int t, x, y, z;
	int nc, nf;
	double b, m;
	int n;
	int type;
} filename_t;

/**
 * @file scatter_test.c
 *
 * Tests the scattering code
 *
 * @author Tadeusz Janowski
 */

/**
 * @brief Sets the configuration to unit gauge
 */
void unit_gauge(suNg_field *gauge){
    int mu;
    int x, y, z, t, ix; for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ ix=ipt(t,x,y,z);
        for (mu=0; mu<4; ++mu) {
            _suNg_unit(*_4FIELD_AT(gauge,ix,mu));
        }
    }
    start_gf_sendrecv(gauge);
    complete_gf_sendrecv(gauge);
}

typedef struct fourvector{
    double v[4];
} fourvec;

/**
 * @brief Adds two four-vectors together replacing the first one with the sum, v1 += v2
 */
static void iadd(fourvec *v1, fourvec *v2)
{
    for(int i=0;i<4;++i)
    {
        v1->v[i] += v2->v[i];
    }
}

/**
 * @brief Multiply four-vector by a real number
 */
static void imul(fourvec *v1, double a)
{
    for(int i=0;i<4;++i)
    {
        v1->v[i] *= a;
    }
}

/**
 * @brief Returns sum over phat^2 + m (part of the propagator)
 */
double f1(fourvec p, double m)
{
    double tmp = 0.0;
    int i;

    for(i=0;i<4;++i)
    {
        tmp += sin(p.v[i]/2)*sin(p.v[i]/2);
    }
    return m + 2*tmp;
}

/**
 * @brief Part of the propagator
 */
static double f2(fourvec v1, fourvec v2)
{
    int i;
    double result = 0.0;
    for(i=0;i<4;++i)
    {
        result += sin(v1.v[i])*sin(v2.v[i]);
    }
    return result;
}

/**
 * @brief Part of the propagator
 */
double b_mu(fourvec p1, int mu){
    return sin(p1.v[mu]);
}

/**
 * @brief Calculates analytic expression for pion 2-point function
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex twopoint(fourvec p, double m,int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    hr_complex res;
    res = 0.;
    double tmp;
    
    for (q1=0;q1<L;++q1) for (q2=0; q2<L; ++q2)  for (q3=0; q3<L; ++q3) for (q41=0; q41<LT; ++q41) for (q42=0; q42<LT; ++q42){
        mom1 = (fourvec) {{q1,q2,q3,((double) q41)*L/LT}};
        iadd(&mom1, &p);
        imul(&mom1, 2.0* PI / L);
        mom2 = (fourvec) {{q1,q2,q3,((double) q42)*L/LT}};
        imul(&mom2, 2.0* PI / L);

        tmp = (f1(mom1,m)*f1(mom2,m) + f2(mom1,mom2)) / ( (SQR(f1(mom1,m)) + f2(mom1,mom1)) * (SQR(f1(mom2,m)) + f2(mom2,mom2) ) );
        res += tmp * cexp(I*(2.0 * PI / LT) * t * (q42 - q41));
    }

    res = 4*res/L/L/L/LT/LT;
    return res;
}
/**
 * @brief Calculates analytic expression for scalar 2-point function
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex twopoint_id(fourvec p, double m,int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    hr_complex res;
    res = 0.;
    double tmp;
    
    for (q1=0;q1<L;++q1) for (q2=0; q2<L; ++q2)  for (q3=0; q3<L; ++q3) for (q41=0; q41<LT; ++q41) for (q42=0; q42<LT; ++q42){
        mom1 = (fourvec) {{q1,q2,q3,((double) q41)*L/LT}};
        iadd(&mom1, &p);
        imul(&mom1, 2.0* PI / L);
        mom2 = (fourvec) {{q1,q2,q3,((double) q42)*L/LT}};
        imul(&mom2, 2.0* PI / L);

        tmp = (f1(mom1,m)*f1(mom2,m) - f2(mom1,mom2)) / ( (SQR(f1(mom1,m)) + f2(mom1,mom1)) * (SQR(f1(mom2,m)) + f2(mom2,mom2) ) );
        res += tmp * cexp(I*(2.0 * PI / LT) * t * (q42 - q41));
    }

    res = 4*res/L/L/L/LT/LT;
    return res;
}



/**
 * @brief Calculates analytic expression for pion 2-point function
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex disc_id(fourvec p, double m,int L, int LT, int t)
{
    fourvec mom1;
    int q1, q2, q3, q41;
    hr_complex res;
    res = 0.;
    double tmp;
    
    for (q1=0;q1<L;++q1) for (q2=0; q2<L; ++q2)  for (q3=0; q3<L; ++q3) for (q41=0; q41<LT; ++q41){
        mom1 = (fourvec) {{q1,q2,q3,((double) q41)*L/LT}};
        iadd(&mom1, &p);
        imul(&mom1, 2.0* PI / L);
    
        tmp = f1(mom1,m) /  (SQR(f1(mom1,m)) + f2(mom1,mom1))  ;
        res += tmp;
    }

    res = 4*res/L/L/L/LT;
    return res;
}

/**
 * @brief Calculates analytic expression for rho 2-point function with gamma3 at source and sink
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex twopoint_rho(fourvec p, double m,int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    hr_complex res;
    res = 0.;
    double tmp;

    for (q1=0;q1<L;++q1) for (q2=0; q2<L; ++q2)  for (q3=0; q3<L; ++q3) for (q41=0; q41<LT; ++q41) for (q42=0; q42<LT; ++q42){
        mom1 = (fourvec) {{q1,q2,q3,((double) q41)*L/LT}};
        iadd(&mom1, &p);
        imul(&mom1, 2.0* PI / L);
        mom2 = (fourvec) {{q1,q2,q3,((double) q42)*L/LT}};
        imul(&mom2, 2.0* PI / L);

        tmp = (f1(mom1,m)*f1(mom2,m) + f2(mom1,mom2) - 2*sin(mom1.v[2])*sin(mom2.v[2])) / ( (SQR(f1(mom1,m)) + f2(mom1,mom1)) * (SQR(f1(mom2,m)) + f2(mom2,mom2) ) );
        res += tmp * cexp(I*(2.0 * PI / LT) * t * (q42 - q41));
    }

    res = 4*res/L/L/L/LT/LT;
    return res;
}

/**
 * @brief Calculates analytic expression for rho 2-point function with gamma1 at source and gamma2 at the sink
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex twopoint_rho12(fourvec p, double m,int L, int LT, int t)
{
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    hr_complex res;
    res=0.;
    double tmp;

    for (q1=0;q1<L;++q1) for (q2=0; q2<L; ++q2)  for (q3=0; q3<L; ++q3) for (q41=0; q41<LT; ++q41) for (q42=0; q42<LT; ++q42){
        mom1 = (fourvec) {{q1,q2,q3,((double) q41)*L/LT}};
        iadd(&mom1, &p);
        imul(&mom1, 2.0* PI / L);
        mom2 = (fourvec) {{q1,q2,q3,((double) q42)*L/LT}};
        imul(&mom2, 2.0* PI / L);

        tmp = ( -(sin(mom1.v[0])*sin(mom2.v[1]) + sin(mom1.v[1])*sin(mom2.v[0]) )) / ( (SQR(f1(mom1,m)) + f2(mom1,mom1)) * (SQR(f1(mom2,m)) + f2(mom2,mom2) ) );
        res += tmp* cexp( I*(2.0 * PI / LT) * t * (q42 - q41));
    }

    res = 4*res/L/L/L/LT/LT;
    return res;
}
#define Q(A,L) (2*PI*(A)/(L))
#define FV(A,B) (fourvec) {{Q(A##1,L), Q(A##2,L), Q(A##3,L), Q(B,LT)}}

/**
 * @brief Calculates analytic expression for the triangle graph with gamma_3 at the sink
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex Triangle(fourvec p, double m, int L, int LT, int t)
{
    fourvec mom[3];
    int q1, q2, q3, q14, q24, q34, i,j;
    hr_complex res;
    res = 0.;
    double numerator, denominator;
    double af1[3];
    double af2[3][3];

    for(q1=0; q1<L; ++q1)  for(q2=0; q2<L; ++q2) for(q3=0; q3<L; ++q3) for (q14=0; q14<LT; ++q14) for (q24=0; q24<LT; ++q24) for (q34=0; q34<LT; ++q34){
        mom[0] = FV(q,q14);
        mom[1] = FV(q,q24);
        mom[2] = (fourvec) {{q1,q2,q3,((double) q34)*L/LT}};
        iadd(&mom[2],& p);
        imul(&mom[2], 2.0* PI / L);

        denominator=1.0;
        for(i=0;i<3;++i)
        {
            af1[i] = f1(mom[i],m);
            for(j=0;j<3;++j)
            {
                af2[i][j] = f2(mom[i],mom[j]);
            }
            denominator *= (SQR(af1[i]) + af2[i][i]);
        }

        numerator \
            = (af1[0]*af1[1] + af2[0][1])*sin(mom[2].v[2]) \
            + (af1[0]*af1[2] + af2[0][2])*sin(mom[1].v[2]) \
            - (af1[1]*af1[2] + af2[1][2])*sin(mom[0].v[2]) ;

        //res += ( sin((double) (t * (q24 - q34)) * 2.0 * PI/LT)  -I*cos((double) (t * (q24 - q34)) * 2.0 * PI/LT))*numerator / denominator;
        res += cexp(I*(t * (q24 - q34) * 2.0 * PI/LT - PI/2.))*numerator / denominator;
    }

    res = 4*res/L/L/L/LT/LT/LT;
    return res;
}



/**
 * @brief Calculates analytic expression for the triangle graph with identity at the sink
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex Triangle_id(fourvec p, double m, int L, int LT, int t)
{
    fourvec mom[3];
    int q1, q2, q3, q14, q24, q34, i,j;
    hr_complex res;
    res = 0.;
    double numerator, denominator;
    double af1[3];
    double af2[3][3];

    for(q1=0; q1<L; ++q1)  for(q2=0; q2<L; ++q2) for(q3=0; q3<L; ++q3) for (q14=0; q14<LT; ++q14) for (q24=0; q24<LT; ++q24) for (q34=0; q34<LT; ++q34){
        mom[0] = FV(q,q14);
        mom[1] = FV(q,q24);
        mom[2] = (fourvec) {{q1,q2,q3,((double) q34)*L/LT}};
        iadd(&mom[2],& p);
        imul(&mom[2], 2.0* PI / L);

        denominator=1.0;
        for(i=0;i<3;++i)
        {
            af1[i] = f1(mom[i],m);
            for(j=0;j<3;++j)
            {
                af2[i][j] = f2(mom[i],mom[j]);
            }
            denominator *= (SQR(af1[i]) + af2[i][i]);
        }


	numerator = (af1[0]*af1[1]*af1[2]) - af1[0]*( af2[1][2] ) + af1[1]*( af2[0][2] ) + af1[2]*( af2[0][1] );

        res += cexp(I*(t * (q24 - q34) * 2.0 * PI/LT ))*numerator / denominator;
    }

    res = 4*res/L/L/L/LT/LT/LT;
    return res;
}




/**
 * @brief Calculates analytic expression for the rectangle pipi->pipi contraction.
 * @param px Momentum at the sink connected to the source
 * @param py Momentum at the source
 * @param pz Momentum at the sink on the vertex opposite to the source
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex R_(fourvec px, fourvec py, fourvec pz, double m, int L, int LT, int t)
{
    fourvec mom[4];
    int q11, q12, q13, q14, q24, q34, q44, i,j;
    hr_complex res;
    res = 0.;
    double numerator, denominator;
    double af1[4];
    double af2[4][4];

    for (q11=0;q11<L;++q11) for (q12=0; q12<L; ++q12)  for (q13=0; q13<L; ++q13) for (q14=0; q14<LT; ++q14) for (q24=0; q24<LT; ++q24) for (q34=0; q34<LT; ++q34) for (q44=0; q44<LT; ++q44){
        mom[0] = FV(q1,q14);
        mom[1] = (fourvec) {{q11,q12,q13,((double) q24 )*L/LT}};
        iadd(&mom[1],& px);
        imul(&mom[1], 2.0* PI / L);
        mom[2] = (fourvec) {{q11,q12,q13,((double) q34 )*L/LT}};
        iadd(&mom[2],& px);
        iadd(&mom[2],& pz);
        imul(&mom[2], 2.0* PI / L);
        mom[3] = (fourvec) {{q11,q12,q13,((double) q44 )*L/LT}};
        iadd(&mom[3],& px);
        iadd(&mom[3],& py);
        iadd(&mom[3],& pz);
        imul(&mom[3], 2.0* PI / L);

        denominator=1.0;
        for(i=0;i<4;++i)
        {
            af1[i] = f1(mom[i],m);
            for(j=0;j<4;++j)
            {
                af2[i][j] = f2(mom[i],mom[j]);
            }
            denominator *= (SQR(af1[i]) + af2[i][i]);
        }

        numerator \
            = af1[0]*af1[1]*af1[2]*af1[3] \
            + af1[0]*af1[1]*af2[2][3] \
            - af1[0]*af1[2]*af2[1][3] \
            + af1[0]*af1[3]*af2[1][2] \
            + af1[1]*af1[2]*af2[0][3] \
            - af1[1]*af1[3]*af2[0][2] \
            + af1[2]*af1[3]*af2[0][1] \
            + af2[0][1] * af2[2][3] \
            - af2[0][2] * af2[1][3] \
            + af2[0][3] * af2[1][2];

              res += cexp(I*(t * (q24-q44)) * 2.0 * PI/LT)* numerator / denominator;
    }

    res = 4*res/L/L/L/LT/LT/LT/LT;
    return res;
}

/**
 * @brief Compares numerical and analytic values of the correlation function.
 * @param mo meson_observable containing the correlation funtion
 * @param corr analytically computed correlation function
 * @param px,py,pz momentum at the sink
 * @param pmax maximum value of p used to generate mo
 * @param tol tolerance, the program returns an error if abs(numeric-analytic)>tol
 * @returns 0 if comparison successful, 1 otherwise
 */
int compare_2pt(meson_observable *mo, hr_complex *corr, int px, int py, int pz, int pmax, double tol ){
    int retval = 0;
    for(int t=0; t<GLB_T; t++){
        double num_re = mo->corr_re[corr_ind(px,py,pz,pmax,t,1,0)];
        double ana_re = creal(corr[t]);

        double num_im = mo->corr_im[corr_ind(px,py,pz,pmax,t,1,0)];
        double ana_im = cimag(corr[t]);
        if(fabs(num_re - ana_re) > tol || fabs(num_im - ana_im)>tol){
            lprintf("TEST",0,"Mismatch, t=%d, numeric = %e + I*(%e), analytic = %e + I*(%e)",t,num_re,num_im,ana_re,ana_im);
            retval = 1;
        }
        else
        {
            lprintf("TEST",0,"Match, t=%d, numeric = %e + I*(%e), analytic = %e + I*(%e)",t,num_re,num_im,ana_re,ana_im);
        }
        
    }
    return retval;
}



int main(int argc,char *argv[])
{
  int return_value=0;
  int ncorr=3;
  double m[256];
  meson_observable **mo_arr;  
  
  fourvec zero_p = (fourvec){{0,0,0,0}};
  
  double threshold = 0.03;

  error(!(GLB_X==GLB_Y && GLB_X==GLB_Z),1,"main", "This test works only for GLB_X=GLB_Y=GLB_Z");

  setup_process(&argc,&argv);

  setup_gauge_fields();
  
  read_input(glb_var.read,get_input_filename());
  read_input(mes_var.read,get_input_filename());
  read_input(rlx_var.read,get_input_filename());

  int numsources = mes_var.nhits;
   
  #if defined(WITH_CLOVER) ||  defined(WITH_EXPCLOVER)
  set_csw(&mes_var.csw);
  #endif
	
   m[0] = atof(mes_var.mstring); 
   init_propagator_eo(1,m,mes_var.precision);
  
   lprintf("MAIN",0,"Boundary conditions: %s\n",mes_var.bc);
   lprintf("MAIN",0,"mass is : %e\n",m[0]);
   lprintf("MAIN",0,"Number of hits : %d \n",numsources);
  

   struct timeval start, end, etime;
   gettimeofday(&start,0);
 
   unit_gauge(u_gauge);
   represent_gauge_field();
   gettimeofday(&start,0);
   //char path[100]="./output/";
   hr_complex Tstoch[GLB_T],Ttheo[GLB_T];
   hr_complex Dstoch[GLB_T],Dtheo[GLB_T];
   hr_complex Cstoch[GLB_T],Ctheo[GLB_T];

   for (int t=0;t < GLB_T ; t++)    {
       Tstoch[t] = 0.0;
       Dstoch[t] = 0.0;
       Cstoch[t] = 0.0;
   }

    mo_arr= (meson_observable**)malloc(sizeof(meson_observable*)*ncorr);
    for(int i=0; i<ncorr; i++) mo_arr[i] =  (meson_observable*) malloc(sizeof(meson_observable));

    init_mo(mo_arr[0],"sigmaconn",GLB_T);
    init_mo(mo_arr[1],"sigmadisc",GLB_T);
    init_mo(mo_arr[2],"T",GLB_T);
   
   for (int n=0;n<numsources;n++)
    {   
        if (n%100==0) lprintf("MAIN",0,"nhits: %d / %d \n/",n,numsources);
        measure_pion_scattering_I0_TS( m,1 ,mes_var.precision,NULL,"test",1,mo_arr);
        
        for (int t=0;t < GLB_T ; t++) 
            {
                Tstoch[t] += (mo_arr[2]->corr_re[corr_ind(0,0,0,0,t,1,0)]+I*mo_arr[2]->corr_im[corr_ind(0,0,0,0,t,1,0)])/(numsources*GLB_X*GLB_Z*GLB_Y); 
                Dstoch[t] += (mo_arr[1]->corr_re[corr_ind(0,0,0,0,t,1,0)]+I*mo_arr[1]->corr_im[corr_ind(0,0,0,0,t,1,0)])/(numsources*GLB_X*GLB_Z*GLB_Y); 
                Cstoch[t] += (mo_arr[0]->corr_re[corr_ind(0,0,0,0,t,1,0)]+I*mo_arr[0]->corr_im[corr_ind(0,0,0,0,t,1,0)])/(numsources*GLB_X*GLB_Z*GLB_Y); 
            }
        for (int i=0;i<ncorr;i++) reset_mo(mo_arr[i]);
        
    }   
 
    for (int t=0;t < GLB_T ; t++)  
    {
      Ttheo[t] =  NF*Triangle_id( zero_p, m[0],GLB_X,GLB_T,t );  // R_(zero_p,zero_p,zero_p,m[0],GLB_X,GLB_T,t);
      Dtheo[t] =  NF*disc_id(zero_p,m[0],GLB_X,GLB_T, t)*0.5;    // factor 2 is due to the eo normalisation of the stochastic sources.
      Ctheo[t] =  NF*twopoint_id(zero_p,m[0],GLB_X,GLB_T, t)*0.5; // source is not even/odd but the io2pt is called with a 0.5 normalisation
    }
  
  
    for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"T analytical T %e + 1I %e  numerical %e + 1I %e\n" ,creal(Ttheo[i]),cimag(Ttheo[i]),creal(Tstoch[i]), cimag(Tstoch[i]));
    for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"diff analytical T %e + 1I %e, and relative %e  \n" ,creal(Ttheo[i]) - creal(Tstoch[i]),cimag(Ttheo[i])- cimag(Tstoch[i]), (creal(Ttheo[i]) - creal(Tstoch[i]))/( creal(Ttheo[i])  ) );

    for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"Disc analytical Disc %e + 1I %e  numerical %e + 1I %e\n" ,creal(Dtheo[i]),cimag(Dtheo[i]),creal(Dstoch[i]), cimag(Dstoch[i]));
    for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"diff analytical Disc %e + 1I %e, and relative %e  \n" ,creal(Dtheo[i]) - creal(Dstoch[i]),cimag(Dtheo[i])- cimag(Dstoch[i]), (creal(Dtheo[i]) - creal(Dstoch[i]))/( creal(Dtheo[i])  ) );

    for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"Conn analytical C %e + 1I %e  numerical %e + 1I %e\n" ,creal(Ctheo[i]),cimag(Ctheo[i]),creal(Cstoch[i]), cimag(Cstoch[i]));
    for (int i=0;i<GLB_T;++i)    lprintf("MAIN",0,"diff analytical C %e + 1I %e, and relative %e  \n" ,creal(Ctheo[i]) - creal(Cstoch[i]),cimag(Ctheo[i])- cimag(Cstoch[i]), (creal(Ctheo[i]) - creal(Cstoch[i]))/( creal(Ctheo[i])  ) );


    for (int i=0;i<GLB_T;i++){

      lprintf("MAIN",0,"T check %e\n",  (creal(Ttheo[i]) - creal(Tstoch[i])) /creal(Ttheo[i])    ) ;
      lprintf("MAIN",0,"D check %e\n",  (creal(Dtheo[i]) - creal(Dstoch[i])) /creal(Dtheo[i])    ) ;
      lprintf("MAIN",0,"Conn check %e\n",  (creal(Ctheo[i]) - creal(Cstoch[i])) /creal(Ctheo[i])    ) ;


      if(  fabs( (creal(Ttheo[i]) - creal(Tstoch[i]))/creal(Ttheo[i]) ) > threshold  ) return_value += 1;
      if(  fabs( (creal(Dtheo[i]) - creal(Dstoch[i]))/creal(Dtheo[i]) ) > threshold  ) return_value += 1;
      if(  fabs( (creal(Ctheo[i]) - creal(Cstoch[i]))/creal(Ctheo[i]) ) > threshold  ) return_value += 1;
    }


    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration : analysed in [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);
    
  
  global_sum_int(&return_value,1);
  lprintf("MAIN", 0, "return_value= %d\n ",  return_value);

  lprintf("DEBUG",0,"ALL done, deallocating\n");
  for (int i=0;i<ncorr;i++) free_mo(mo_arr[i]);

  finalize_process();

  return return_value;
}
