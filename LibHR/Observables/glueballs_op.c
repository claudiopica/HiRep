#include <stdlib.h>
#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "utils.h"
#include "glueballs.h"

#include <string.h>
#define npaths 229
static double PI = 3.141592653589793238462643383279502884197;
static double complex *mom_def_Re_tr_paths = NULL;
static double complex *mom_def_Im_tr_paths = NULL;
static double complex *path_storage = NULL;
int **direct_spatial_rotations()
{
    int i;
    int **res = malloc(sizeof(int *) * 48);
    int *res1 = malloc(sizeof(int *) * 48 * 4);
    for (i = 0; i < 48; i++)
        res[i] = res1 + 4 * i;

    res1[0] = 0;
    res1[1] = 1;
    res1[2] = 2;
    res1[3] = 3;
    res1[4] = 0;
    res1[5] = 2;
    res1[6] = 1;
    res1[7] = -3;
    res1[8] = 0;
    res1[9] = -2;
    res1[10] = -1;
    res1[11] = -3;
    res1[12] = 0;
    res1[13] = -1;
    res1[14] = 3;
    res1[15] = 2;
    res1[16] = 0;
    res1[17] = -1;
    res1[18] = -3;
    res1[19] = -2;
    res1[20] = 0;
    res1[21] = 3;
    res1[22] = -2;
    res1[23] = 1;
    res1[24] = 0;
    res1[25] = -3;
    res1[26] = -2;
    res1[27] = -1;
    res1[28] = 0;
    res1[29] = -2;
    res1[30] = -3;
    res1[31] = 1;
    res1[32] = 0;
    res1[33] = 3;
    res1[34] = -1;
    res1[35] = -2;
    res1[36] = 0;
    res1[37] = -3;
    res1[38] = -1;
    res1[39] = 2;
    res1[40] = 0;
    res1[41] = -2;
    res1[42] = 3;
    res1[43] = -1;
    res1[44] = 0;
    res1[45] = 3;
    res1[46] = 1;
    res1[47] = 2;
    res1[48] = 0;
    res1[49] = 2;
    res1[50] = 3;
    res1[51] = 1;
    res1[52] = 0;
    res1[53] = 2;
    res1[54] = -3;
    res1[55] = -1;
    res1[56] = 0;
    res1[57] = -3;
    res1[58] = 1;
    res1[59] = -2;
    res1[60] = 0;
    res1[61] = -2;
    res1[62] = 1;
    res1[63] = 3;
    res1[64] = 0;
    res1[65] = 2;
    res1[66] = -1;
    res1[67] = 3;
    res1[68] = 0;
    res1[69] = 3;
    res1[70] = 2;
    res1[71] = -1;
    res1[72] = 0;
    res1[73] = -3;
    res1[74] = 2;
    res1[75] = 1;
    res1[76] = 0;
    res1[77] = 1;
    res1[78] = 3;
    res1[79] = -2;
    res1[80] = 0;
    res1[81] = 1;
    res1[82] = -3;
    res1[83] = 2;
    res1[84] = 0;
    res1[85] = -1;
    res1[86] = -2;
    res1[87] = 3;
    res1[88] = 0;
    res1[89] = -1;
    res1[90] = 2;
    res1[91] = -3;
    res1[92] = 0;
    res1[93] = 1;
    res1[94] = -2;
    res1[95] = -3;
    res1[96] = 0;
    res1[97] = -1;
    res1[98] = -2;
    res1[99] = -3;
    res1[100] = 0;
    res1[101] = -2;
    res1[102] = -1;
    res1[103] = 3;
    res1[104] = 0;
    res1[105] = 2;
    res1[106] = 1;
    res1[107] = 3;
    res1[108] = 0;
    res1[109] = 1;
    res1[110] = -3;
    res1[111] = -2;
    res1[112] = 0;
    res1[113] = 1;
    res1[114] = 3;
    res1[115] = 2;
    res1[116] = 0;
    res1[117] = -3;
    res1[118] = 2;
    res1[119] = -1;
    res1[120] = 0;
    res1[121] = 3;
    res1[122] = 2;
    res1[123] = 1;
    res1[124] = 0;
    res1[125] = 2;
    res1[126] = 3;
    res1[127] = -1;
    res1[128] = 0;
    res1[129] = -3;
    res1[130] = 1;
    res1[131] = 2;
    res1[132] = 0;
    res1[133] = 3;
    res1[134] = 1;
    res1[135] = -2;
    res1[136] = 0;
    res1[137] = 2;
    res1[138] = -3;
    res1[139] = 1;
    res1[140] = 0;
    res1[141] = -3;
    res1[142] = -1;
    res1[143] = -2;
    res1[144] = 0;
    res1[145] = -2;
    res1[146] = -3;
    res1[147] = -1;
    res1[148] = 0;
    res1[149] = -2;
    res1[150] = 3;
    res1[151] = 1;
    res1[152] = 0;
    res1[153] = 3;
    res1[154] = -1;
    res1[155] = 2;
    res1[156] = 0;
    res1[157] = 2;
    res1[158] = -1;
    res1[159] = -3;
    res1[160] = 0;
    res1[161] = -2;
    res1[162] = 1;
    res1[163] = -3;
    res1[164] = 0;
    res1[165] = -3;
    res1[166] = -2;
    res1[167] = 1;
    res1[168] = 0;
    res1[169] = 3;
    res1[170] = -2;
    res1[171] = -1;
    res1[172] = 0;
    res1[173] = -1;
    res1[174] = -3;
    res1[175] = 2;
    res1[176] = 0;
    res1[177] = -1;
    res1[178] = 3;
    res1[179] = -2;
    res1[180] = 0;
    res1[181] = 1;
    res1[182] = 2;
    res1[183] = -3;
    res1[184] = 0;
    res1[185] = 1;
    res1[186] = -2;
    res1[187] = 3;
    res1[188] = 0;
    res1[189] = -1;
    res1[190] = 2;
    res1[191] = 3;
    return res;
}
int **inverse_spatial_rotations()
{
    int i;
    int **res = malloc(sizeof(int *) * 48);
    int *res1 = malloc(sizeof(int *) * 48 * 4);
    for (i = 0; i < 48; i++)
        res[i] = res1 + 4 * i;
    res1[0] = 0;
    res1[1] = 1;
    res1[2] = 2;
    res1[3] = 3;
    res1[4] = 0;
    res1[5] = 2;
    res1[6] = 1;
    res1[7] = -3;
    res1[8] = 0;
    res1[9] = -2;
    res1[10] = -1;
    res1[11] = -3;
    res1[12] = 0;
    res1[13] = -1;
    res1[14] = 3;
    res1[15] = 2;
    res1[16] = 0;
    res1[17] = -1;
    res1[18] = -3;
    res1[19] = -2;
    res1[20] = 0;
    res1[21] = 3;
    res1[22] = -2;
    res1[23] = 1;
    res1[24] = 0;
    res1[25] = -3;
    res1[26] = -2;
    res1[27] = -1;
    res1[28] = 0;
    res1[29] = 3;
    res1[30] = -1;
    res1[31] = -2;
    res1[32] = 0;
    res1[33] = -2;
    res1[34] = -3;
    res1[35] = 1;
    res1[36] = 0;
    res1[37] = -2;
    res1[38] = 3;
    res1[39] = -1;
    res1[40] = 0;
    res1[41] = -3;
    res1[42] = -1;
    res1[43] = 2;
    res1[44] = 0;
    res1[45] = 2;
    res1[46] = 3;
    res1[47] = 1;
    res1[48] = 0;
    res1[49] = 3;
    res1[50] = 1;
    res1[51] = 2;
    res1[52] = 0;
    res1[53] = -3;
    res1[54] = 1;
    res1[55] = -2;
    res1[56] = 0;
    res1[57] = 2;
    res1[58] = -3;
    res1[59] = -1;
    res1[60] = 0;
    res1[61] = 2;
    res1[62] = -1;
    res1[63] = 3;
    res1[64] = 0;
    res1[65] = -2;
    res1[66] = 1;
    res1[67] = 3;
    res1[68] = 0;
    res1[69] = -3;
    res1[70] = 2;
    res1[71] = 1;
    res1[72] = 0;
    res1[73] = 3;
    res1[74] = 2;
    res1[75] = -1;
    res1[76] = 0;
    res1[77] = 1;
    res1[78] = -3;
    res1[79] = 2;
    res1[80] = 0;
    res1[81] = 1;
    res1[82] = 3;
    res1[83] = -2;
    res1[84] = 0;
    res1[85] = -1;
    res1[86] = -2;
    res1[87] = 3;
    res1[88] = 0;
    res1[89] = -1;
    res1[90] = 2;
    res1[91] = -3;
    res1[92] = 0;
    res1[93] = 1;
    res1[94] = -2;
    res1[95] = -3;
    res1[96] = 0;
    res1[97] = -1;
    res1[98] = -2;
    res1[99] = -3;
    res1[100] = 0;
    res1[101] = -2;
    res1[102] = -1;
    res1[103] = 3;
    res1[104] = 0;
    res1[105] = 2;
    res1[106] = 1;
    res1[107] = 3;
    res1[108] = 0;
    res1[109] = 1;
    res1[110] = -3;
    res1[111] = -2;
    res1[112] = 0;
    res1[113] = 1;
    res1[114] = 3;
    res1[115] = 2;
    res1[116] = 0;
    res1[117] = -3;
    res1[118] = 2;
    res1[119] = -1;
    res1[120] = 0;
    res1[121] = 3;
    res1[122] = 2;
    res1[123] = 1;
    res1[124] = 0;
    res1[125] = -3;
    res1[126] = 1;
    res1[127] = 2;
    res1[128] = 0;
    res1[129] = 2;
    res1[130] = 3;
    res1[131] = -1;
    res1[132] = 0;
    res1[133] = 2;
    res1[134] = -3;
    res1[135] = 1;
    res1[136] = 0;
    res1[137] = 3;
    res1[138] = 1;
    res1[139] = -2;
    res1[140] = 0;
    res1[141] = -2;
    res1[142] = -3;
    res1[143] = -1;
    res1[144] = 0;
    res1[145] = -3;
    res1[146] = -1;
    res1[147] = -2;
    res1[148] = 0;
    res1[149] = 3;
    res1[150] = -1;
    res1[151] = 2;
    res1[152] = 0;
    res1[153] = -2;
    res1[154] = 3;
    res1[155] = 1;
    res1[156] = 0;
    res1[157] = -2;
    res1[158] = 1;
    res1[159] = -3;
    res1[160] = 0;
    res1[161] = 2;
    res1[162] = -1;
    res1[163] = -3;
    res1[164] = 0;
    res1[165] = 3;
    res1[166] = -2;
    res1[167] = -1;
    res1[168] = 0;
    res1[169] = -3;
    res1[170] = -2;
    res1[171] = 1;
    res1[172] = 0;
    res1[173] = -1;
    res1[174] = 3;
    res1[175] = -2;
    res1[176] = 0;
    res1[177] = -1;
    res1[178] = -3;
    res1[179] = 2;
    res1[180] = 0;
    res1[181] = 1;
    res1[182] = 2;
    res1[183] = -3;
    res1[184] = 0;
    res1[185] = 1;
    res1[186] = -2;
    res1[187] = 3;
    res1[188] = 0;
    res1[189] = -1;
    res1[190] = 2;
    res1[191] = 3;
    return res;
}
static double complex path9(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path11(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -mom_def_Im_tr_paths[9] + mom_def_Im_tr_paths[9] - mom_def_Im_tr_paths[11] + mom_def_Im_tr_paths[11];
}

static double complex path0(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path2(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex c0;
static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c0 * mom_def_Re_tr_paths[0] + (2.) * mom_def_Re_tr_paths[2];
}

static double complex path3(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path8(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex c1;
static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c1 * mom_def_Re_tr_paths[3] + (2.) * mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[9] + mom_def_Re_tr_paths[9] + mom_def_Re_tr_paths[11] + mom_def_Re_tr_paths[11];
}

static double complex path13(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path18(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[13] + (2.) * mom_def_Re_tr_paths[18];
}

static double complex path30(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path43(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[30] + (2.) * mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * mom_def_Im_tr_paths[0] + (2.) * mom_def_Im_tr_paths[0];
}

static double complex path5(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c0 * mom_def_Im_tr_paths[3] + (2.) * mom_def_Im_tr_paths[5];
}

static double complex path17(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path21(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path22(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex c2;
static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -c2 * mom_def_Im_tr_paths[17] - c2 * mom_def_Im_tr_paths[18] + mom_def_Im_tr_paths[21] + mom_def_Im_tr_paths[22];
}

static double complex path29(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path33(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path34(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[29] + mom_def_Im_tr_paths[30] - c0 * mom_def_Im_tr_paths[33] - c0 * mom_def_Im_tr_paths[34];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[0] + (2.) * mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c0 * mom_def_Re_tr_paths[3] + (2.) * mom_def_Re_tr_paths[5];
}

static double complex path12(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * c1 * mom_def_Re_tr_paths[11] + (2.) * c1 * mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +c2 * mom_def_Re_tr_paths[17] + c2 * mom_def_Re_tr_paths[18] + mom_def_Re_tr_paths[21] + mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[29] + mom_def_Re_tr_paths[30] + c0 * mom_def_Re_tr_paths[33] + c0 * mom_def_Re_tr_paths[34];
}

static double complex path14(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path15(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path19(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -c0 * mom_def_Im_tr_paths[14] + mom_def_Im_tr_paths[15] - c0 * mom_def_Im_tr_paths[18] + mom_def_Im_tr_paths[19];
}

static double complex path26(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path27(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path31(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[26] - c2 * mom_def_Im_tr_paths[27] + mom_def_Im_tr_paths[30] - c2 * mom_def_Im_tr_paths[31];
}

static double complex path1(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c0 * mom_def_Re_tr_paths[0] + (2.) * mom_def_Re_tr_paths[1];
}

static double complex path4(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c0 * mom_def_Re_tr_paths[3] + (2.) * mom_def_Re_tr_paths[4];
}

static double complex path10(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex c3;
static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[10] + (2.) * c3 * mom_def_Re_tr_paths[11];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +c0 * mom_def_Re_tr_paths[14] + mom_def_Re_tr_paths[15] + c0 * mom_def_Re_tr_paths[18] + mom_def_Re_tr_paths[19];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[26] + c2 * mom_def_Re_tr_paths[27] + mom_def_Re_tr_paths[30] + c2 * mom_def_Re_tr_paths[31];
}

static double complex path97(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path98(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path99(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path100(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path101(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path102(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path103(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path104(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path105(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path106(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path107(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path108(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path109(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path110(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path111(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path112(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path113(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path114(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path115(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path116(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path117(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path118(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path119(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path120(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path121(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path122(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path123(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path124(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path125(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path126(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path127(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path128(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path129(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path130(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path131(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path132(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path133(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path134(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path135(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path136(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path137(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path138(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path139(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path140(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path141(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path142(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path143(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path144(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[97] + mom_def_Im_tr_paths[98] + mom_def_Im_tr_paths[99] + mom_def_Im_tr_paths[100] + mom_def_Im_tr_paths[101] + mom_def_Im_tr_paths[102] + mom_def_Im_tr_paths[103] + mom_def_Im_tr_paths[104] + mom_def_Im_tr_paths[105] + mom_def_Im_tr_paths[106] + mom_def_Im_tr_paths[107] + mom_def_Im_tr_paths[108] + mom_def_Im_tr_paths[109] + mom_def_Im_tr_paths[110] + mom_def_Im_tr_paths[111] + mom_def_Im_tr_paths[112] + mom_def_Im_tr_paths[113] + mom_def_Im_tr_paths[114] + mom_def_Im_tr_paths[115] + mom_def_Im_tr_paths[116] + mom_def_Im_tr_paths[117] + mom_def_Im_tr_paths[118] + mom_def_Im_tr_paths[119] + mom_def_Im_tr_paths[120] + mom_def_Im_tr_paths[121] + mom_def_Im_tr_paths[122] + mom_def_Im_tr_paths[123] + mom_def_Im_tr_paths[124] + mom_def_Im_tr_paths[125] + mom_def_Im_tr_paths[126] + mom_def_Im_tr_paths[127] + mom_def_Im_tr_paths[128] + mom_def_Im_tr_paths[129] + mom_def_Im_tr_paths[130] + mom_def_Im_tr_paths[131] + mom_def_Im_tr_paths[132] + mom_def_Im_tr_paths[133] + mom_def_Im_tr_paths[134] + mom_def_Im_tr_paths[135] + mom_def_Im_tr_paths[136] + mom_def_Im_tr_paths[137] + mom_def_Im_tr_paths[138] + mom_def_Im_tr_paths[139] + mom_def_Im_tr_paths[140] + mom_def_Im_tr_paths[141] + mom_def_Im_tr_paths[142] + mom_def_Im_tr_paths[143] + mom_def_Im_tr_paths[144];
}

static double complex path145(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path146(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path147(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path148(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path149(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path150(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path151(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path152(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path153(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path154(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path155(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path156(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path157(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path158(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path159(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path160(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path161(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path162(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path163(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path164(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path165(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path166(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path167(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path168(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path169(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path170(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path171(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path172(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path173(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path174(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path175(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path176(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path177(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path178(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path179(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path180(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path181(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path182(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path183(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path184(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path185(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path186(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path187(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path188(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path189(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path190(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path191(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path192(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[145] + mom_def_Im_tr_paths[146] + mom_def_Im_tr_paths[147] + mom_def_Im_tr_paths[148] + mom_def_Im_tr_paths[149] + mom_def_Im_tr_paths[150] + mom_def_Im_tr_paths[151] + mom_def_Im_tr_paths[152] + mom_def_Im_tr_paths[153] + mom_def_Im_tr_paths[154] + mom_def_Im_tr_paths[155] + mom_def_Im_tr_paths[156] + mom_def_Im_tr_paths[157] + mom_def_Im_tr_paths[158] + mom_def_Im_tr_paths[159] + mom_def_Im_tr_paths[160] + mom_def_Im_tr_paths[161] + mom_def_Im_tr_paths[162] + mom_def_Im_tr_paths[163] + mom_def_Im_tr_paths[164] + mom_def_Im_tr_paths[165] + mom_def_Im_tr_paths[166] + mom_def_Im_tr_paths[167] + mom_def_Im_tr_paths[168] + mom_def_Im_tr_paths[169] + mom_def_Im_tr_paths[170] + mom_def_Im_tr_paths[171] + mom_def_Im_tr_paths[172] + mom_def_Im_tr_paths[173] + mom_def_Im_tr_paths[174] + mom_def_Im_tr_paths[175] + mom_def_Im_tr_paths[176] + mom_def_Im_tr_paths[177] + mom_def_Im_tr_paths[178] + mom_def_Im_tr_paths[179] + mom_def_Im_tr_paths[180] + mom_def_Im_tr_paths[181] + mom_def_Im_tr_paths[182] + mom_def_Im_tr_paths[183] + mom_def_Im_tr_paths[184] + mom_def_Im_tr_paths[185] + mom_def_Im_tr_paths[186] + mom_def_Im_tr_paths[187] + mom_def_Im_tr_paths[188] + mom_def_Im_tr_paths[189] + mom_def_Im_tr_paths[190] + mom_def_Im_tr_paths[191] + mom_def_Im_tr_paths[192];
}

static double complex path205(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path206(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path207(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path208(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path209(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path210(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path211(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path212(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path213(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path214(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path215(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path216(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path217(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path218(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path219(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path220(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path221(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path222(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path223(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path224(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path225(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path226(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path227(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path228(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Im_tr_paths[205] + (2.) * mom_def_Im_tr_paths[206] + (2.) * mom_def_Im_tr_paths[207] + (2.) * mom_def_Im_tr_paths[208] + (2.) * mom_def_Im_tr_paths[209] + (2.) * mom_def_Im_tr_paths[210] + (2.) * mom_def_Im_tr_paths[211] + (2.) * mom_def_Im_tr_paths[212] + (2.) * mom_def_Im_tr_paths[213] + (2.) * mom_def_Im_tr_paths[214] + (2.) * mom_def_Im_tr_paths[215] + (2.) * mom_def_Im_tr_paths[216] + (2.) * mom_def_Im_tr_paths[217] + (2.) * mom_def_Im_tr_paths[218] + (2.) * mom_def_Im_tr_paths[219] + (2.) * mom_def_Im_tr_paths[220] + (2.) * mom_def_Im_tr_paths[221] + (2.) * mom_def_Im_tr_paths[222] + (2.) * mom_def_Im_tr_paths[223] + (2.) * mom_def_Im_tr_paths[224] + (2.) * mom_def_Im_tr_paths[225] + (2.) * mom_def_Im_tr_paths[226] + (2.) * mom_def_Im_tr_paths[227] + (2.) * mom_def_Im_tr_paths[228];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(16.) * mom_def_Re_tr_paths[0] + (16.) * mom_def_Re_tr_paths[1] + (16.) * mom_def_Re_tr_paths[2];
}

static double complex path6(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path7(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(8.) * mom_def_Re_tr_paths[3] + (8.) * mom_def_Re_tr_paths[4] + (8.) * mom_def_Re_tr_paths[5] + (8.) * mom_def_Re_tr_paths[6] + (8.) * mom_def_Re_tr_paths[7] + (8.) * mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(12.) * mom_def_Re_tr_paths[9] + (12.) * mom_def_Re_tr_paths[10] + (12.) * mom_def_Re_tr_paths[11] + (12.) * mom_def_Re_tr_paths[12];
}

static double complex path16(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path20(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path23(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path24(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(4.) * mom_def_Re_tr_paths[13] + (4.) * mom_def_Re_tr_paths[14] + (4.) * mom_def_Re_tr_paths[15] + (4.) * mom_def_Re_tr_paths[16] + (4.) * mom_def_Re_tr_paths[17] + (4.) * mom_def_Re_tr_paths[18] + (4.) * mom_def_Re_tr_paths[19] + (4.) * mom_def_Re_tr_paths[20] + (4.) * mom_def_Re_tr_paths[21] + (4.) * mom_def_Re_tr_paths[22] + (4.) * mom_def_Re_tr_paths[23] + (4.) * mom_def_Re_tr_paths[24];
}

static double complex path25(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path28(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path32(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path35(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path36(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path37(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path38(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path39(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path40(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path41(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path42(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path44(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path45(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path46(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path47(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path48(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[25] + (2.) * mom_def_Re_tr_paths[26] + (2.) * mom_def_Re_tr_paths[27] + (2.) * mom_def_Re_tr_paths[28] + (2.) * mom_def_Re_tr_paths[29] + (2.) * mom_def_Re_tr_paths[30] + (2.) * mom_def_Re_tr_paths[31] + (2.) * mom_def_Re_tr_paths[32] + (2.) * mom_def_Re_tr_paths[33] + (2.) * mom_def_Re_tr_paths[34] + (2.) * mom_def_Re_tr_paths[35] + (2.) * mom_def_Re_tr_paths[36] + (2.) * mom_def_Re_tr_paths[37] + (2.) * mom_def_Re_tr_paths[38] + (2.) * mom_def_Re_tr_paths[39] + (2.) * mom_def_Re_tr_paths[40] + (2.) * mom_def_Re_tr_paths[41] + (2.) * mom_def_Re_tr_paths[42] + (2.) * mom_def_Re_tr_paths[43] + (2.) * mom_def_Re_tr_paths[44] + (2.) * mom_def_Re_tr_paths[45] + (2.) * mom_def_Re_tr_paths[46] + (2.) * mom_def_Re_tr_paths[47] + (2.) * mom_def_Re_tr_paths[48];
}

static double complex path49(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path50(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path51(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path52(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path53(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path54(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path55(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path56(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path57(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path58(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path59(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path60(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path61(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path62(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path63(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path64(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path65(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path66(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path67(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path68(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path69(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path70(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path71(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path72(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_6(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[49] + (2.) * mom_def_Re_tr_paths[50] + (2.) * mom_def_Re_tr_paths[51] + (2.) * mom_def_Re_tr_paths[52] + (2.) * mom_def_Re_tr_paths[53] + (2.) * mom_def_Re_tr_paths[54] + (2.) * mom_def_Re_tr_paths[55] + (2.) * mom_def_Re_tr_paths[56] + (2.) * mom_def_Re_tr_paths[57] + (2.) * mom_def_Re_tr_paths[58] + (2.) * mom_def_Re_tr_paths[59] + (2.) * mom_def_Re_tr_paths[60] + (2.) * mom_def_Re_tr_paths[61] + (2.) * mom_def_Re_tr_paths[62] + (2.) * mom_def_Re_tr_paths[63] + (2.) * mom_def_Re_tr_paths[64] + (2.) * mom_def_Re_tr_paths[65] + (2.) * mom_def_Re_tr_paths[66] + (2.) * mom_def_Re_tr_paths[67] + (2.) * mom_def_Re_tr_paths[68] + (2.) * mom_def_Re_tr_paths[69] + (2.) * mom_def_Re_tr_paths[70] + (2.) * mom_def_Re_tr_paths[71] + (2.) * mom_def_Re_tr_paths[72];
}

static double complex path73(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path74(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path75(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path76(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path77(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path78(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path79(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path80(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path81(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path82(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path83(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path84(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path85(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path86(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path87(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path88(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path89(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path90(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path91(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path92(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path93(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path94(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path95(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path96(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_7(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[73] + (2.) * mom_def_Re_tr_paths[74] + (2.) * mom_def_Re_tr_paths[75] + (2.) * mom_def_Re_tr_paths[76] + (2.) * mom_def_Re_tr_paths[77] + (2.) * mom_def_Re_tr_paths[78] + (2.) * mom_def_Re_tr_paths[79] + (2.) * mom_def_Re_tr_paths[80] + (2.) * mom_def_Re_tr_paths[81] + (2.) * mom_def_Re_tr_paths[82] + (2.) * mom_def_Re_tr_paths[83] + (2.) * mom_def_Re_tr_paths[84] + (2.) * mom_def_Re_tr_paths[85] + (2.) * mom_def_Re_tr_paths[86] + (2.) * mom_def_Re_tr_paths[87] + (2.) * mom_def_Re_tr_paths[88] + (2.) * mom_def_Re_tr_paths[89] + (2.) * mom_def_Re_tr_paths[90] + (2.) * mom_def_Re_tr_paths[91] + (2.) * mom_def_Re_tr_paths[92] + (2.) * mom_def_Re_tr_paths[93] + (2.) * mom_def_Re_tr_paths[94] + (2.) * mom_def_Re_tr_paths[95] + (2.) * mom_def_Re_tr_paths[96];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_8(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[97] + mom_def_Re_tr_paths[98] + mom_def_Re_tr_paths[99] + mom_def_Re_tr_paths[100] + mom_def_Re_tr_paths[101] + mom_def_Re_tr_paths[102] + mom_def_Re_tr_paths[103] + mom_def_Re_tr_paths[104] + mom_def_Re_tr_paths[105] + mom_def_Re_tr_paths[106] + mom_def_Re_tr_paths[107] + mom_def_Re_tr_paths[108] + mom_def_Re_tr_paths[109] + mom_def_Re_tr_paths[110] + mom_def_Re_tr_paths[111] + mom_def_Re_tr_paths[112] + mom_def_Re_tr_paths[113] + mom_def_Re_tr_paths[114] + mom_def_Re_tr_paths[115] + mom_def_Re_tr_paths[116] + mom_def_Re_tr_paths[117] + mom_def_Re_tr_paths[118] + mom_def_Re_tr_paths[119] + mom_def_Re_tr_paths[120] + mom_def_Re_tr_paths[121] + mom_def_Re_tr_paths[122] + mom_def_Re_tr_paths[123] + mom_def_Re_tr_paths[124] + mom_def_Re_tr_paths[125] + mom_def_Re_tr_paths[126] + mom_def_Re_tr_paths[127] + mom_def_Re_tr_paths[128] + mom_def_Re_tr_paths[129] + mom_def_Re_tr_paths[130] + mom_def_Re_tr_paths[131] + mom_def_Re_tr_paths[132] + mom_def_Re_tr_paths[133] + mom_def_Re_tr_paths[134] + mom_def_Re_tr_paths[135] + mom_def_Re_tr_paths[136] + mom_def_Re_tr_paths[137] + mom_def_Re_tr_paths[138] + mom_def_Re_tr_paths[139] + mom_def_Re_tr_paths[140] + mom_def_Re_tr_paths[141] + mom_def_Re_tr_paths[142] + mom_def_Re_tr_paths[143] + mom_def_Re_tr_paths[144];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_9(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[145] + mom_def_Re_tr_paths[146] + mom_def_Re_tr_paths[147] + mom_def_Re_tr_paths[148] + mom_def_Re_tr_paths[149] + mom_def_Re_tr_paths[150] + mom_def_Re_tr_paths[151] + mom_def_Re_tr_paths[152] + mom_def_Re_tr_paths[153] + mom_def_Re_tr_paths[154] + mom_def_Re_tr_paths[155] + mom_def_Re_tr_paths[156] + mom_def_Re_tr_paths[157] + mom_def_Re_tr_paths[158] + mom_def_Re_tr_paths[159] + mom_def_Re_tr_paths[160] + mom_def_Re_tr_paths[161] + mom_def_Re_tr_paths[162] + mom_def_Re_tr_paths[163] + mom_def_Re_tr_paths[164] + mom_def_Re_tr_paths[165] + mom_def_Re_tr_paths[166] + mom_def_Re_tr_paths[167] + mom_def_Re_tr_paths[168] + mom_def_Re_tr_paths[169] + mom_def_Re_tr_paths[170] + mom_def_Re_tr_paths[171] + mom_def_Re_tr_paths[172] + mom_def_Re_tr_paths[173] + mom_def_Re_tr_paths[174] + mom_def_Re_tr_paths[175] + mom_def_Re_tr_paths[176] + mom_def_Re_tr_paths[177] + mom_def_Re_tr_paths[178] + mom_def_Re_tr_paths[179] + mom_def_Re_tr_paths[180] + mom_def_Re_tr_paths[181] + mom_def_Re_tr_paths[182] + mom_def_Re_tr_paths[183] + mom_def_Re_tr_paths[184] + mom_def_Re_tr_paths[185] + mom_def_Re_tr_paths[186] + mom_def_Re_tr_paths[187] + mom_def_Re_tr_paths[188] + mom_def_Re_tr_paths[189] + mom_def_Re_tr_paths[190] + mom_def_Re_tr_paths[191] + mom_def_Re_tr_paths[192];
}

static double complex path193(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path194(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path195(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path196(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 1);
    w2 = pu_gauge_wrk(site, 1);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path197(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path198(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path199(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path200(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path201(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path202(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path203(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path204(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_10(double complex *op_out)
{
    *op_out = +(4.) * mom_def_Re_tr_paths[193] + (4.) * mom_def_Re_tr_paths[194] + (4.) * mom_def_Re_tr_paths[195] + (4.) * mom_def_Re_tr_paths[196] + (4.) * mom_def_Re_tr_paths[197] + (4.) * mom_def_Re_tr_paths[198] + (4.) * mom_def_Re_tr_paths[199] + (4.) * mom_def_Re_tr_paths[200] + (4.) * mom_def_Re_tr_paths[201] + (4.) * mom_def_Re_tr_paths[202] + (4.) * mom_def_Re_tr_paths[203] + (4.) * mom_def_Re_tr_paths[204];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_11(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[205] + (2.) * mom_def_Re_tr_paths[206] + (2.) * mom_def_Re_tr_paths[207] + (2.) * mom_def_Re_tr_paths[208] + (2.) * mom_def_Re_tr_paths[209] + (2.) * mom_def_Re_tr_paths[210] + (2.) * mom_def_Re_tr_paths[211] + (2.) * mom_def_Re_tr_paths[212] + (2.) * mom_def_Re_tr_paths[213] + (2.) * mom_def_Re_tr_paths[214] + (2.) * mom_def_Re_tr_paths[215] + (2.) * mom_def_Re_tr_paths[216] + (2.) * mom_def_Re_tr_paths[217] + (2.) * mom_def_Re_tr_paths[218] + (2.) * mom_def_Re_tr_paths[219] + (2.) * mom_def_Re_tr_paths[220] + (2.) * mom_def_Re_tr_paths[221] + (2.) * mom_def_Re_tr_paths[222] + (2.) * mom_def_Re_tr_paths[223] + (2.) * mom_def_Re_tr_paths[224] + (2.) * mom_def_Re_tr_paths[225] + (2.) * mom_def_Re_tr_paths[226] + (2.) * mom_def_Re_tr_paths[227] + (2.) * mom_def_Re_tr_paths[228];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(8.) * mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * mom_def_Re_tr_paths[3] + (4.) * mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[9] + (2.) * mom_def_Re_tr_paths[10] + (2.) * c2 * mom_def_Re_tr_paths[11] + (2.) * mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[14] + (2.) * mom_def_Re_tr_paths[18] + (2.) * mom_def_Re_tr_paths[22] + (2.) * mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[26] + (2.) * mom_def_Re_tr_paths[30] + (2.) * mom_def_Re_tr_paths[34] + (2.) * mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -c2 * mom_def_Im_tr_paths[14] + mom_def_Im_tr_paths[15] - c2 * mom_def_Im_tr_paths[18] + mom_def_Im_tr_paths[19];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[26] - c0 * mom_def_Im_tr_paths[27] + mom_def_Im_tr_paths[30] - c0 * mom_def_Im_tr_paths[31];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[0] + (2.) * mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[3] + (2.) * mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[10] + (2.) * c1 * mom_def_Re_tr_paths[11];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +c2 * mom_def_Re_tr_paths[14] + mom_def_Re_tr_paths[15] + c2 * mom_def_Re_tr_paths[18] + mom_def_Re_tr_paths[19];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[26] + c0 * mom_def_Re_tr_paths[27] + mom_def_Re_tr_paths[30] + c0 * mom_def_Re_tr_paths[31];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[25] + mom_def_Im_tr_paths[26] + mom_def_Im_tr_paths[29] + mom_def_Im_tr_paths[30] + mom_def_Im_tr_paths[43] + mom_def_Im_tr_paths[44] + mom_def_Im_tr_paths[47] + mom_def_Im_tr_paths[48];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * c2 * mom_def_Re_tr_paths[0] + (4.) * c2 * mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * c2 * mom_def_Re_tr_paths[3] + (4.) * c2 * mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[9] + (2.) * mom_def_Re_tr_paths[10] + (2.) * c0 * mom_def_Re_tr_paths[11] + (2.) * c0 * mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[13] + (2.) * c2 * mom_def_Re_tr_paths[14] + (2.) * c2 * mom_def_Re_tr_paths[17] + (2.) * c2 * mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[25] + mom_def_Re_tr_paths[26] + mom_def_Re_tr_paths[29] + mom_def_Re_tr_paths[30] + mom_def_Re_tr_paths[43] + mom_def_Re_tr_paths[44] + mom_def_Re_tr_paths[47] + mom_def_Re_tr_paths[48];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[10] - mom_def_Im_tr_paths[10] + mom_def_Im_tr_paths[11] - mom_def_Im_tr_paths[11];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = -c2 * mom_def_Im_tr_paths[14] - c2 * mom_def_Im_tr_paths[16] - c2 * mom_def_Im_tr_paths[18] - c2 * mom_def_Im_tr_paths[20];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[26] + mom_def_Im_tr_paths[28] + mom_def_Im_tr_paths[30] + mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[0] + (2.) * c2 * mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[3] + (2.) * c2 * mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[10] + mom_def_Re_tr_paths[10] + mom_def_Re_tr_paths[11] + mom_def_Re_tr_paths[11];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +c2 * mom_def_Re_tr_paths[14] + c2 * mom_def_Re_tr_paths[16] + c2 * mom_def_Re_tr_paths[18] + c2 * mom_def_Re_tr_paths[20];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[26] + mom_def_Re_tr_paths[28] + mom_def_Re_tr_paths[30] + mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * mom_def_Im_tr_paths[0] + (2.) * mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Im_tr_paths[3] + (2.) * mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -c0 * mom_def_Im_tr_paths[17] - c0 * mom_def_Im_tr_paths[18] + mom_def_Im_tr_paths[21] + mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[29] + mom_def_Im_tr_paths[30] - c2 * mom_def_Im_tr_paths[33] - c2 * mom_def_Im_tr_paths[34];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[0] + (2.) * mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[3] + (2.) * mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * c3 * mom_def_Re_tr_paths[11] + (2.) * c3 * mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +c0 * mom_def_Re_tr_paths[17] + c0 * mom_def_Re_tr_paths[18] + mom_def_Re_tr_paths[21] + mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[29] + mom_def_Re_tr_paths[30] + c2 * mom_def_Re_tr_paths[33] + c2 * mom_def_Re_tr_paths[34];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -mom_def_Im_tr_paths[9] + mom_def_Im_tr_paths[9] - mom_def_Im_tr_paths[11] + mom_def_Im_tr_paths[11];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[0] + (2.) * mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c3 * mom_def_Re_tr_paths[3] + (2.) * mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[9] + mom_def_Re_tr_paths[9] + mom_def_Re_tr_paths[11] + mom_def_Re_tr_paths[11];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[13] + (2.) * mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[30] + (2.) * mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * c2 * mom_def_Re_tr_paths[0] + (4.) * c2 * mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * c3 * mom_def_Re_tr_paths[3] + (4.) * c3 * mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[9] + (2.) * c2 * mom_def_Re_tr_paths[10] + (2.) * c2 * mom_def_Re_tr_paths[11] + (2.) * c2 * mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[17] + (2.) * mom_def_Re_tr_paths[18] + (2.) * mom_def_Re_tr_paths[19] + (2.) * mom_def_Re_tr_paths[20];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[29] + (2.) * mom_def_Re_tr_paths[30] + (2.) * mom_def_Re_tr_paths[31] + (2.) * mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c2 * mom_def_Re_tr_paths[0] + (2.) * c2 * mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c3 * mom_def_Re_tr_paths[3] + (2.) * c3 * mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * c3 * mom_def_Re_tr_paths[9] + (2.) * c3 * mom_def_Re_tr_paths[11];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +(4.) * mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * mom_def_Re_tr_paths[30] + (2.) * mom_def_Re_tr_paths[48];
}

static double complex c4;
static void OP_oneTr_p_1_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * c4 * mom_def_Im_tr_paths[3] + (2.) * c4 * mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = -mom_def_Im_tr_paths[11] + mom_def_Im_tr_paths[11] - mom_def_Im_tr_paths[12] + mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -c2 * mom_def_Im_tr_paths[17] - c2 * mom_def_Im_tr_paths[18] - c2 * mom_def_Im_tr_paths[23] - c2 * mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +mom_def_Im_tr_paths[29] + mom_def_Im_tr_paths[30] + mom_def_Im_tr_paths[37] + mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * c3 * mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * c4 * mom_def_Re_tr_paths[3] + (2.) * c4 * mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[11] + mom_def_Re_tr_paths[11] + mom_def_Re_tr_paths[12] + mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +c2 * mom_def_Re_tr_paths[17] + c2 * mom_def_Re_tr_paths[18] + c2 * mom_def_Re_tr_paths[23] + c2 * mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_5(double complex *op_out)
{
    *op_out = +mom_def_Re_tr_paths[29] + mom_def_Re_tr_paths[30] + mom_def_Re_tr_paths[37] + mom_def_Re_tr_paths[38];
}

static int last_t = -10;
void request_space_paths_evaluation() { last_t = -10; }

static void eval_time_momentum_glueball_paths(int t, int px, int py, int pz)
{
    int n_x, n_y, n_z, idx, in;
    double complex ce = 0.;
    if (path_storage == NULL)
    {
        c0 = cexp(I * PI * (2. / GLB_X));
        c1 = cexp(I * PI * (4. / GLB_X));
        c2 = cexp(I * PI * (-2. / GLB_X));
        c3 = cexp(I * PI * (-4. / GLB_X));
        c4 = cexp(I * PI * (-6. / GLB_X));
        path_storage = malloc(npaths * X * Y * Z * sizeof(double complex));
        mom_def_Re_tr_paths = malloc(npaths * sizeof(double complex));
        mom_def_Im_tr_paths = malloc(npaths * sizeof(double complex));
    }

    for (in = 0; in < npaths; in++)
    {
        mom_def_Re_tr_paths[in] = 0.;
        mom_def_Im_tr_paths[in] = 0.;
    }
    if (t != last_t)
    {
        last_t = t;
        for (n_x = 0; n_x < X; n_x++)
            for (n_y = 0; n_y < Y; n_y++)
                for (n_z = 0; n_z < Z; n_z++)
                {
                    in = ipt(t, n_x, n_y, n_z);
                    ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
                    idx = npaths * (n_x + X * (n_y + Y * n_z));
                    path_storage[0 + idx] = path0(in);
                    mom_def_Re_tr_paths[0] += ce * creal(path_storage[0 + idx]);
                    mom_def_Im_tr_paths[0] += I * ce * cimag(path_storage[0 + idx]);
                    path_storage[1 + idx] = path1(in);
                    mom_def_Re_tr_paths[1] += ce * creal(path_storage[1 + idx]);
                    mom_def_Im_tr_paths[1] += I * ce * cimag(path_storage[1 + idx]);
                    path_storage[2 + idx] = path2(in);
                    mom_def_Re_tr_paths[2] += ce * creal(path_storage[2 + idx]);
                    mom_def_Im_tr_paths[2] += I * ce * cimag(path_storage[2 + idx]);
                    path_storage[3 + idx] = path3(in);
                    mom_def_Re_tr_paths[3] += ce * creal(path_storage[3 + idx]);
                    mom_def_Im_tr_paths[3] += I * ce * cimag(path_storage[3 + idx]);
                    path_storage[4 + idx] = path4(in);
                    mom_def_Re_tr_paths[4] += ce * creal(path_storage[4 + idx]);
                    mom_def_Im_tr_paths[4] += I * ce * cimag(path_storage[4 + idx]);
                    path_storage[5 + idx] = path5(in);
                    mom_def_Re_tr_paths[5] += ce * creal(path_storage[5 + idx]);
                    mom_def_Im_tr_paths[5] += I * ce * cimag(path_storage[5 + idx]);
                    path_storage[6 + idx] = path6(in);
                    mom_def_Re_tr_paths[6] += ce * creal(path_storage[6 + idx]);
                    mom_def_Im_tr_paths[6] += I * ce * cimag(path_storage[6 + idx]);
                    path_storage[7 + idx] = path7(in);
                    mom_def_Re_tr_paths[7] += ce * creal(path_storage[7 + idx]);
                    mom_def_Im_tr_paths[7] += I * ce * cimag(path_storage[7 + idx]);
                    path_storage[8 + idx] = path8(in);
                    mom_def_Re_tr_paths[8] += ce * creal(path_storage[8 + idx]);
                    mom_def_Im_tr_paths[8] += I * ce * cimag(path_storage[8 + idx]);
                    path_storage[9 + idx] = path9(in);
                    mom_def_Re_tr_paths[9] += ce * creal(path_storage[9 + idx]);
                    mom_def_Im_tr_paths[9] += I * ce * cimag(path_storage[9 + idx]);
                    path_storage[10 + idx] = path10(in);
                    mom_def_Re_tr_paths[10] += ce * creal(path_storage[10 + idx]);
                    mom_def_Im_tr_paths[10] += I * ce * cimag(path_storage[10 + idx]);
                    path_storage[11 + idx] = path11(in);
                    mom_def_Re_tr_paths[11] += ce * creal(path_storage[11 + idx]);
                    mom_def_Im_tr_paths[11] += I * ce * cimag(path_storage[11 + idx]);
                    path_storage[12 + idx] = path12(in);
                    mom_def_Re_tr_paths[12] += ce * creal(path_storage[12 + idx]);
                    mom_def_Im_tr_paths[12] += I * ce * cimag(path_storage[12 + idx]);
                    path_storage[13 + idx] = path13(in);
                    mom_def_Re_tr_paths[13] += ce * creal(path_storage[13 + idx]);
                    mom_def_Im_tr_paths[13] += I * ce * cimag(path_storage[13 + idx]);
                    path_storage[14 + idx] = path14(in);
                    mom_def_Re_tr_paths[14] += ce * creal(path_storage[14 + idx]);
                    mom_def_Im_tr_paths[14] += I * ce * cimag(path_storage[14 + idx]);
                    path_storage[15 + idx] = path15(in);
                    mom_def_Re_tr_paths[15] += ce * creal(path_storage[15 + idx]);
                    mom_def_Im_tr_paths[15] += I * ce * cimag(path_storage[15 + idx]);
                    path_storage[16 + idx] = path16(in);
                    mom_def_Re_tr_paths[16] += ce * creal(path_storage[16 + idx]);
                    mom_def_Im_tr_paths[16] += I * ce * cimag(path_storage[16 + idx]);
                    path_storage[17 + idx] = path17(in);
                    mom_def_Re_tr_paths[17] += ce * creal(path_storage[17 + idx]);
                    mom_def_Im_tr_paths[17] += I * ce * cimag(path_storage[17 + idx]);
                    path_storage[18 + idx] = path18(in);
                    mom_def_Re_tr_paths[18] += ce * creal(path_storage[18 + idx]);
                    mom_def_Im_tr_paths[18] += I * ce * cimag(path_storage[18 + idx]);
                    path_storage[19 + idx] = path19(in);
                    mom_def_Re_tr_paths[19] += ce * creal(path_storage[19 + idx]);
                    mom_def_Im_tr_paths[19] += I * ce * cimag(path_storage[19 + idx]);
                    path_storage[20 + idx] = path20(in);
                    mom_def_Re_tr_paths[20] += ce * creal(path_storage[20 + idx]);
                    mom_def_Im_tr_paths[20] += I * ce * cimag(path_storage[20 + idx]);
                    path_storage[21 + idx] = path21(in);
                    mom_def_Re_tr_paths[21] += ce * creal(path_storage[21 + idx]);
                    mom_def_Im_tr_paths[21] += I * ce * cimag(path_storage[21 + idx]);
                    path_storage[22 + idx] = path22(in);
                    mom_def_Re_tr_paths[22] += ce * creal(path_storage[22 + idx]);
                    mom_def_Im_tr_paths[22] += I * ce * cimag(path_storage[22 + idx]);
                    path_storage[23 + idx] = path23(in);
                    mom_def_Re_tr_paths[23] += ce * creal(path_storage[23 + idx]);
                    mom_def_Im_tr_paths[23] += I * ce * cimag(path_storage[23 + idx]);
                    path_storage[24 + idx] = path24(in);
                    mom_def_Re_tr_paths[24] += ce * creal(path_storage[24 + idx]);
                    mom_def_Im_tr_paths[24] += I * ce * cimag(path_storage[24 + idx]);
                    path_storage[25 + idx] = path25(in);
                    mom_def_Re_tr_paths[25] += ce * creal(path_storage[25 + idx]);
                    mom_def_Im_tr_paths[25] += I * ce * cimag(path_storage[25 + idx]);
                    path_storage[26 + idx] = path26(in);
                    mom_def_Re_tr_paths[26] += ce * creal(path_storage[26 + idx]);
                    mom_def_Im_tr_paths[26] += I * ce * cimag(path_storage[26 + idx]);
                    path_storage[27 + idx] = path27(in);
                    mom_def_Re_tr_paths[27] += ce * creal(path_storage[27 + idx]);
                    mom_def_Im_tr_paths[27] += I * ce * cimag(path_storage[27 + idx]);
                    path_storage[28 + idx] = path28(in);
                    mom_def_Re_tr_paths[28] += ce * creal(path_storage[28 + idx]);
                    mom_def_Im_tr_paths[28] += I * ce * cimag(path_storage[28 + idx]);
                    path_storage[29 + idx] = path29(in);
                    mom_def_Re_tr_paths[29] += ce * creal(path_storage[29 + idx]);
                    mom_def_Im_tr_paths[29] += I * ce * cimag(path_storage[29 + idx]);
                    path_storage[30 + idx] = path30(in);
                    mom_def_Re_tr_paths[30] += ce * creal(path_storage[30 + idx]);
                    mom_def_Im_tr_paths[30] += I * ce * cimag(path_storage[30 + idx]);
                    path_storage[31 + idx] = path31(in);
                    mom_def_Re_tr_paths[31] += ce * creal(path_storage[31 + idx]);
                    mom_def_Im_tr_paths[31] += I * ce * cimag(path_storage[31 + idx]);
                    path_storage[32 + idx] = path32(in);
                    mom_def_Re_tr_paths[32] += ce * creal(path_storage[32 + idx]);
                    mom_def_Im_tr_paths[32] += I * ce * cimag(path_storage[32 + idx]);
                    path_storage[33 + idx] = path33(in);
                    mom_def_Re_tr_paths[33] += ce * creal(path_storage[33 + idx]);
                    mom_def_Im_tr_paths[33] += I * ce * cimag(path_storage[33 + idx]);
                    path_storage[34 + idx] = path34(in);
                    mom_def_Re_tr_paths[34] += ce * creal(path_storage[34 + idx]);
                    mom_def_Im_tr_paths[34] += I * ce * cimag(path_storage[34 + idx]);
                    path_storage[35 + idx] = path35(in);
                    mom_def_Re_tr_paths[35] += ce * creal(path_storage[35 + idx]);
                    mom_def_Im_tr_paths[35] += I * ce * cimag(path_storage[35 + idx]);
                    path_storage[36 + idx] = path36(in);
                    mom_def_Re_tr_paths[36] += ce * creal(path_storage[36 + idx]);
                    mom_def_Im_tr_paths[36] += I * ce * cimag(path_storage[36 + idx]);
                    path_storage[37 + idx] = path37(in);
                    mom_def_Re_tr_paths[37] += ce * creal(path_storage[37 + idx]);
                    mom_def_Im_tr_paths[37] += I * ce * cimag(path_storage[37 + idx]);
                    path_storage[38 + idx] = path38(in);
                    mom_def_Re_tr_paths[38] += ce * creal(path_storage[38 + idx]);
                    mom_def_Im_tr_paths[38] += I * ce * cimag(path_storage[38 + idx]);
                    path_storage[39 + idx] = path39(in);
                    mom_def_Re_tr_paths[39] += ce * creal(path_storage[39 + idx]);
                    mom_def_Im_tr_paths[39] += I * ce * cimag(path_storage[39 + idx]);
                    path_storage[40 + idx] = path40(in);
                    mom_def_Re_tr_paths[40] += ce * creal(path_storage[40 + idx]);
                    mom_def_Im_tr_paths[40] += I * ce * cimag(path_storage[40 + idx]);
                    path_storage[41 + idx] = path41(in);
                    mom_def_Re_tr_paths[41] += ce * creal(path_storage[41 + idx]);
                    mom_def_Im_tr_paths[41] += I * ce * cimag(path_storage[41 + idx]);
                    path_storage[42 + idx] = path42(in);
                    mom_def_Re_tr_paths[42] += ce * creal(path_storage[42 + idx]);
                    mom_def_Im_tr_paths[42] += I * ce * cimag(path_storage[42 + idx]);
                    path_storage[43 + idx] = path43(in);
                    mom_def_Re_tr_paths[43] += ce * creal(path_storage[43 + idx]);
                    mom_def_Im_tr_paths[43] += I * ce * cimag(path_storage[43 + idx]);
                    path_storage[44 + idx] = path44(in);
                    mom_def_Re_tr_paths[44] += ce * creal(path_storage[44 + idx]);
                    mom_def_Im_tr_paths[44] += I * ce * cimag(path_storage[44 + idx]);
                    path_storage[45 + idx] = path45(in);
                    mom_def_Re_tr_paths[45] += ce * creal(path_storage[45 + idx]);
                    mom_def_Im_tr_paths[45] += I * ce * cimag(path_storage[45 + idx]);
                    path_storage[46 + idx] = path46(in);
                    mom_def_Re_tr_paths[46] += ce * creal(path_storage[46 + idx]);
                    mom_def_Im_tr_paths[46] += I * ce * cimag(path_storage[46 + idx]);
                    path_storage[47 + idx] = path47(in);
                    mom_def_Re_tr_paths[47] += ce * creal(path_storage[47 + idx]);
                    mom_def_Im_tr_paths[47] += I * ce * cimag(path_storage[47 + idx]);
                    path_storage[48 + idx] = path48(in);
                    mom_def_Re_tr_paths[48] += ce * creal(path_storage[48 + idx]);
                    mom_def_Im_tr_paths[48] += I * ce * cimag(path_storage[48 + idx]);
                    path_storage[49 + idx] = path49(in);
                    mom_def_Re_tr_paths[49] += ce * creal(path_storage[49 + idx]);
                    mom_def_Im_tr_paths[49] += I * ce * cimag(path_storage[49 + idx]);
                    path_storage[50 + idx] = path50(in);
                    mom_def_Re_tr_paths[50] += ce * creal(path_storage[50 + idx]);
                    mom_def_Im_tr_paths[50] += I * ce * cimag(path_storage[50 + idx]);
                    path_storage[51 + idx] = path51(in);
                    mom_def_Re_tr_paths[51] += ce * creal(path_storage[51 + idx]);
                    mom_def_Im_tr_paths[51] += I * ce * cimag(path_storage[51 + idx]);
                    path_storage[52 + idx] = path52(in);
                    mom_def_Re_tr_paths[52] += ce * creal(path_storage[52 + idx]);
                    mom_def_Im_tr_paths[52] += I * ce * cimag(path_storage[52 + idx]);
                    path_storage[53 + idx] = path53(in);
                    mom_def_Re_tr_paths[53] += ce * creal(path_storage[53 + idx]);
                    mom_def_Im_tr_paths[53] += I * ce * cimag(path_storage[53 + idx]);
                    path_storage[54 + idx] = path54(in);
                    mom_def_Re_tr_paths[54] += ce * creal(path_storage[54 + idx]);
                    mom_def_Im_tr_paths[54] += I * ce * cimag(path_storage[54 + idx]);
                    path_storage[55 + idx] = path55(in);
                    mom_def_Re_tr_paths[55] += ce * creal(path_storage[55 + idx]);
                    mom_def_Im_tr_paths[55] += I * ce * cimag(path_storage[55 + idx]);
                    path_storage[56 + idx] = path56(in);
                    mom_def_Re_tr_paths[56] += ce * creal(path_storage[56 + idx]);
                    mom_def_Im_tr_paths[56] += I * ce * cimag(path_storage[56 + idx]);
                    path_storage[57 + idx] = path57(in);
                    mom_def_Re_tr_paths[57] += ce * creal(path_storage[57 + idx]);
                    mom_def_Im_tr_paths[57] += I * ce * cimag(path_storage[57 + idx]);
                    path_storage[58 + idx] = path58(in);
                    mom_def_Re_tr_paths[58] += ce * creal(path_storage[58 + idx]);
                    mom_def_Im_tr_paths[58] += I * ce * cimag(path_storage[58 + idx]);
                    path_storage[59 + idx] = path59(in);
                    mom_def_Re_tr_paths[59] += ce * creal(path_storage[59 + idx]);
                    mom_def_Im_tr_paths[59] += I * ce * cimag(path_storage[59 + idx]);
                    path_storage[60 + idx] = path60(in);
                    mom_def_Re_tr_paths[60] += ce * creal(path_storage[60 + idx]);
                    mom_def_Im_tr_paths[60] += I * ce * cimag(path_storage[60 + idx]);
                    path_storage[61 + idx] = path61(in);
                    mom_def_Re_tr_paths[61] += ce * creal(path_storage[61 + idx]);
                    mom_def_Im_tr_paths[61] += I * ce * cimag(path_storage[61 + idx]);
                    path_storage[62 + idx] = path62(in);
                    mom_def_Re_tr_paths[62] += ce * creal(path_storage[62 + idx]);
                    mom_def_Im_tr_paths[62] += I * ce * cimag(path_storage[62 + idx]);
                    path_storage[63 + idx] = path63(in);
                    mom_def_Re_tr_paths[63] += ce * creal(path_storage[63 + idx]);
                    mom_def_Im_tr_paths[63] += I * ce * cimag(path_storage[63 + idx]);
                    path_storage[64 + idx] = path64(in);
                    mom_def_Re_tr_paths[64] += ce * creal(path_storage[64 + idx]);
                    mom_def_Im_tr_paths[64] += I * ce * cimag(path_storage[64 + idx]);
                    path_storage[65 + idx] = path65(in);
                    mom_def_Re_tr_paths[65] += ce * creal(path_storage[65 + idx]);
                    mom_def_Im_tr_paths[65] += I * ce * cimag(path_storage[65 + idx]);
                    path_storage[66 + idx] = path66(in);
                    mom_def_Re_tr_paths[66] += ce * creal(path_storage[66 + idx]);
                    mom_def_Im_tr_paths[66] += I * ce * cimag(path_storage[66 + idx]);
                    path_storage[67 + idx] = path67(in);
                    mom_def_Re_tr_paths[67] += ce * creal(path_storage[67 + idx]);
                    mom_def_Im_tr_paths[67] += I * ce * cimag(path_storage[67 + idx]);
                    path_storage[68 + idx] = path68(in);
                    mom_def_Re_tr_paths[68] += ce * creal(path_storage[68 + idx]);
                    mom_def_Im_tr_paths[68] += I * ce * cimag(path_storage[68 + idx]);
                    path_storage[69 + idx] = path69(in);
                    mom_def_Re_tr_paths[69] += ce * creal(path_storage[69 + idx]);
                    mom_def_Im_tr_paths[69] += I * ce * cimag(path_storage[69 + idx]);
                    path_storage[70 + idx] = path70(in);
                    mom_def_Re_tr_paths[70] += ce * creal(path_storage[70 + idx]);
                    mom_def_Im_tr_paths[70] += I * ce * cimag(path_storage[70 + idx]);
                    path_storage[71 + idx] = path71(in);
                    mom_def_Re_tr_paths[71] += ce * creal(path_storage[71 + idx]);
                    mom_def_Im_tr_paths[71] += I * ce * cimag(path_storage[71 + idx]);
                    path_storage[72 + idx] = path72(in);
                    mom_def_Re_tr_paths[72] += ce * creal(path_storage[72 + idx]);
                    mom_def_Im_tr_paths[72] += I * ce * cimag(path_storage[72 + idx]);
                    path_storage[73 + idx] = path73(in);
                    mom_def_Re_tr_paths[73] += ce * creal(path_storage[73 + idx]);
                    mom_def_Im_tr_paths[73] += I * ce * cimag(path_storage[73 + idx]);
                    path_storage[74 + idx] = path74(in);
                    mom_def_Re_tr_paths[74] += ce * creal(path_storage[74 + idx]);
                    mom_def_Im_tr_paths[74] += I * ce * cimag(path_storage[74 + idx]);
                    path_storage[75 + idx] = path75(in);
                    mom_def_Re_tr_paths[75] += ce * creal(path_storage[75 + idx]);
                    mom_def_Im_tr_paths[75] += I * ce * cimag(path_storage[75 + idx]);
                    path_storage[76 + idx] = path76(in);
                    mom_def_Re_tr_paths[76] += ce * creal(path_storage[76 + idx]);
                    mom_def_Im_tr_paths[76] += I * ce * cimag(path_storage[76 + idx]);
                    path_storage[77 + idx] = path77(in);
                    mom_def_Re_tr_paths[77] += ce * creal(path_storage[77 + idx]);
                    mom_def_Im_tr_paths[77] += I * ce * cimag(path_storage[77 + idx]);
                    path_storage[78 + idx] = path78(in);
                    mom_def_Re_tr_paths[78] += ce * creal(path_storage[78 + idx]);
                    mom_def_Im_tr_paths[78] += I * ce * cimag(path_storage[78 + idx]);
                    path_storage[79 + idx] = path79(in);
                    mom_def_Re_tr_paths[79] += ce * creal(path_storage[79 + idx]);
                    mom_def_Im_tr_paths[79] += I * ce * cimag(path_storage[79 + idx]);
                    path_storage[80 + idx] = path80(in);
                    mom_def_Re_tr_paths[80] += ce * creal(path_storage[80 + idx]);
                    mom_def_Im_tr_paths[80] += I * ce * cimag(path_storage[80 + idx]);
                    path_storage[81 + idx] = path81(in);
                    mom_def_Re_tr_paths[81] += ce * creal(path_storage[81 + idx]);
                    mom_def_Im_tr_paths[81] += I * ce * cimag(path_storage[81 + idx]);
                    path_storage[82 + idx] = path82(in);
                    mom_def_Re_tr_paths[82] += ce * creal(path_storage[82 + idx]);
                    mom_def_Im_tr_paths[82] += I * ce * cimag(path_storage[82 + idx]);
                    path_storage[83 + idx] = path83(in);
                    mom_def_Re_tr_paths[83] += ce * creal(path_storage[83 + idx]);
                    mom_def_Im_tr_paths[83] += I * ce * cimag(path_storage[83 + idx]);
                    path_storage[84 + idx] = path84(in);
                    mom_def_Re_tr_paths[84] += ce * creal(path_storage[84 + idx]);
                    mom_def_Im_tr_paths[84] += I * ce * cimag(path_storage[84 + idx]);
                    path_storage[85 + idx] = path85(in);
                    mom_def_Re_tr_paths[85] += ce * creal(path_storage[85 + idx]);
                    mom_def_Im_tr_paths[85] += I * ce * cimag(path_storage[85 + idx]);
                    path_storage[86 + idx] = path86(in);
                    mom_def_Re_tr_paths[86] += ce * creal(path_storage[86 + idx]);
                    mom_def_Im_tr_paths[86] += I * ce * cimag(path_storage[86 + idx]);
                    path_storage[87 + idx] = path87(in);
                    mom_def_Re_tr_paths[87] += ce * creal(path_storage[87 + idx]);
                    mom_def_Im_tr_paths[87] += I * ce * cimag(path_storage[87 + idx]);
                    path_storage[88 + idx] = path88(in);
                    mom_def_Re_tr_paths[88] += ce * creal(path_storage[88 + idx]);
                    mom_def_Im_tr_paths[88] += I * ce * cimag(path_storage[88 + idx]);
                    path_storage[89 + idx] = path89(in);
                    mom_def_Re_tr_paths[89] += ce * creal(path_storage[89 + idx]);
                    mom_def_Im_tr_paths[89] += I * ce * cimag(path_storage[89 + idx]);
                    path_storage[90 + idx] = path90(in);
                    mom_def_Re_tr_paths[90] += ce * creal(path_storage[90 + idx]);
                    mom_def_Im_tr_paths[90] += I * ce * cimag(path_storage[90 + idx]);
                    path_storage[91 + idx] = path91(in);
                    mom_def_Re_tr_paths[91] += ce * creal(path_storage[91 + idx]);
                    mom_def_Im_tr_paths[91] += I * ce * cimag(path_storage[91 + idx]);
                    path_storage[92 + idx] = path92(in);
                    mom_def_Re_tr_paths[92] += ce * creal(path_storage[92 + idx]);
                    mom_def_Im_tr_paths[92] += I * ce * cimag(path_storage[92 + idx]);
                    path_storage[93 + idx] = path93(in);
                    mom_def_Re_tr_paths[93] += ce * creal(path_storage[93 + idx]);
                    mom_def_Im_tr_paths[93] += I * ce * cimag(path_storage[93 + idx]);
                    path_storage[94 + idx] = path94(in);
                    mom_def_Re_tr_paths[94] += ce * creal(path_storage[94 + idx]);
                    mom_def_Im_tr_paths[94] += I * ce * cimag(path_storage[94 + idx]);
                    path_storage[95 + idx] = path95(in);
                    mom_def_Re_tr_paths[95] += ce * creal(path_storage[95 + idx]);
                    mom_def_Im_tr_paths[95] += I * ce * cimag(path_storage[95 + idx]);
                    path_storage[96 + idx] = path96(in);
                    mom_def_Re_tr_paths[96] += ce * creal(path_storage[96 + idx]);
                    mom_def_Im_tr_paths[96] += I * ce * cimag(path_storage[96 + idx]);
                    path_storage[97 + idx] = path97(in);
                    mom_def_Re_tr_paths[97] += ce * creal(path_storage[97 + idx]);
                    mom_def_Im_tr_paths[97] += I * ce * cimag(path_storage[97 + idx]);
                    path_storage[98 + idx] = path98(in);
                    mom_def_Re_tr_paths[98] += ce * creal(path_storage[98 + idx]);
                    mom_def_Im_tr_paths[98] += I * ce * cimag(path_storage[98 + idx]);
                    path_storage[99 + idx] = path99(in);
                    mom_def_Re_tr_paths[99] += ce * creal(path_storage[99 + idx]);
                    mom_def_Im_tr_paths[99] += I * ce * cimag(path_storage[99 + idx]);
                    path_storage[100 + idx] = path100(in);
                    mom_def_Re_tr_paths[100] += ce * creal(path_storage[100 + idx]);
                    mom_def_Im_tr_paths[100] += I * ce * cimag(path_storage[100 + idx]);
                    path_storage[101 + idx] = path101(in);
                    mom_def_Re_tr_paths[101] += ce * creal(path_storage[101 + idx]);
                    mom_def_Im_tr_paths[101] += I * ce * cimag(path_storage[101 + idx]);
                    path_storage[102 + idx] = path102(in);
                    mom_def_Re_tr_paths[102] += ce * creal(path_storage[102 + idx]);
                    mom_def_Im_tr_paths[102] += I * ce * cimag(path_storage[102 + idx]);
                    path_storage[103 + idx] = path103(in);
                    mom_def_Re_tr_paths[103] += ce * creal(path_storage[103 + idx]);
                    mom_def_Im_tr_paths[103] += I * ce * cimag(path_storage[103 + idx]);
                    path_storage[104 + idx] = path104(in);
                    mom_def_Re_tr_paths[104] += ce * creal(path_storage[104 + idx]);
                    mom_def_Im_tr_paths[104] += I * ce * cimag(path_storage[104 + idx]);
                    path_storage[105 + idx] = path105(in);
                    mom_def_Re_tr_paths[105] += ce * creal(path_storage[105 + idx]);
                    mom_def_Im_tr_paths[105] += I * ce * cimag(path_storage[105 + idx]);
                    path_storage[106 + idx] = path106(in);
                    mom_def_Re_tr_paths[106] += ce * creal(path_storage[106 + idx]);
                    mom_def_Im_tr_paths[106] += I * ce * cimag(path_storage[106 + idx]);
                    path_storage[107 + idx] = path107(in);
                    mom_def_Re_tr_paths[107] += ce * creal(path_storage[107 + idx]);
                    mom_def_Im_tr_paths[107] += I * ce * cimag(path_storage[107 + idx]);
                    path_storage[108 + idx] = path108(in);
                    mom_def_Re_tr_paths[108] += ce * creal(path_storage[108 + idx]);
                    mom_def_Im_tr_paths[108] += I * ce * cimag(path_storage[108 + idx]);
                    path_storage[109 + idx] = path109(in);
                    mom_def_Re_tr_paths[109] += ce * creal(path_storage[109 + idx]);
                    mom_def_Im_tr_paths[109] += I * ce * cimag(path_storage[109 + idx]);
                    path_storage[110 + idx] = path110(in);
                    mom_def_Re_tr_paths[110] += ce * creal(path_storage[110 + idx]);
                    mom_def_Im_tr_paths[110] += I * ce * cimag(path_storage[110 + idx]);
                    path_storage[111 + idx] = path111(in);
                    mom_def_Re_tr_paths[111] += ce * creal(path_storage[111 + idx]);
                    mom_def_Im_tr_paths[111] += I * ce * cimag(path_storage[111 + idx]);
                    path_storage[112 + idx] = path112(in);
                    mom_def_Re_tr_paths[112] += ce * creal(path_storage[112 + idx]);
                    mom_def_Im_tr_paths[112] += I * ce * cimag(path_storage[112 + idx]);
                    path_storage[113 + idx] = path113(in);
                    mom_def_Re_tr_paths[113] += ce * creal(path_storage[113 + idx]);
                    mom_def_Im_tr_paths[113] += I * ce * cimag(path_storage[113 + idx]);
                    path_storage[114 + idx] = path114(in);
                    mom_def_Re_tr_paths[114] += ce * creal(path_storage[114 + idx]);
                    mom_def_Im_tr_paths[114] += I * ce * cimag(path_storage[114 + idx]);
                    path_storage[115 + idx] = path115(in);
                    mom_def_Re_tr_paths[115] += ce * creal(path_storage[115 + idx]);
                    mom_def_Im_tr_paths[115] += I * ce * cimag(path_storage[115 + idx]);
                    path_storage[116 + idx] = path116(in);
                    mom_def_Re_tr_paths[116] += ce * creal(path_storage[116 + idx]);
                    mom_def_Im_tr_paths[116] += I * ce * cimag(path_storage[116 + idx]);
                    path_storage[117 + idx] = path117(in);
                    mom_def_Re_tr_paths[117] += ce * creal(path_storage[117 + idx]);
                    mom_def_Im_tr_paths[117] += I * ce * cimag(path_storage[117 + idx]);
                    path_storage[118 + idx] = path118(in);
                    mom_def_Re_tr_paths[118] += ce * creal(path_storage[118 + idx]);
                    mom_def_Im_tr_paths[118] += I * ce * cimag(path_storage[118 + idx]);
                    path_storage[119 + idx] = path119(in);
                    mom_def_Re_tr_paths[119] += ce * creal(path_storage[119 + idx]);
                    mom_def_Im_tr_paths[119] += I * ce * cimag(path_storage[119 + idx]);
                    path_storage[120 + idx] = path120(in);
                    mom_def_Re_tr_paths[120] += ce * creal(path_storage[120 + idx]);
                    mom_def_Im_tr_paths[120] += I * ce * cimag(path_storage[120 + idx]);
                    path_storage[121 + idx] = path121(in);
                    mom_def_Re_tr_paths[121] += ce * creal(path_storage[121 + idx]);
                    mom_def_Im_tr_paths[121] += I * ce * cimag(path_storage[121 + idx]);
                    path_storage[122 + idx] = path122(in);
                    mom_def_Re_tr_paths[122] += ce * creal(path_storage[122 + idx]);
                    mom_def_Im_tr_paths[122] += I * ce * cimag(path_storage[122 + idx]);
                    path_storage[123 + idx] = path123(in);
                    mom_def_Re_tr_paths[123] += ce * creal(path_storage[123 + idx]);
                    mom_def_Im_tr_paths[123] += I * ce * cimag(path_storage[123 + idx]);
                    path_storage[124 + idx] = path124(in);
                    mom_def_Re_tr_paths[124] += ce * creal(path_storage[124 + idx]);
                    mom_def_Im_tr_paths[124] += I * ce * cimag(path_storage[124 + idx]);
                    path_storage[125 + idx] = path125(in);
                    mom_def_Re_tr_paths[125] += ce * creal(path_storage[125 + idx]);
                    mom_def_Im_tr_paths[125] += I * ce * cimag(path_storage[125 + idx]);
                    path_storage[126 + idx] = path126(in);
                    mom_def_Re_tr_paths[126] += ce * creal(path_storage[126 + idx]);
                    mom_def_Im_tr_paths[126] += I * ce * cimag(path_storage[126 + idx]);
                    path_storage[127 + idx] = path127(in);
                    mom_def_Re_tr_paths[127] += ce * creal(path_storage[127 + idx]);
                    mom_def_Im_tr_paths[127] += I * ce * cimag(path_storage[127 + idx]);
                    path_storage[128 + idx] = path128(in);
                    mom_def_Re_tr_paths[128] += ce * creal(path_storage[128 + idx]);
                    mom_def_Im_tr_paths[128] += I * ce * cimag(path_storage[128 + idx]);
                    path_storage[129 + idx] = path129(in);
                    mom_def_Re_tr_paths[129] += ce * creal(path_storage[129 + idx]);
                    mom_def_Im_tr_paths[129] += I * ce * cimag(path_storage[129 + idx]);
                    path_storage[130 + idx] = path130(in);
                    mom_def_Re_tr_paths[130] += ce * creal(path_storage[130 + idx]);
                    mom_def_Im_tr_paths[130] += I * ce * cimag(path_storage[130 + idx]);
                    path_storage[131 + idx] = path131(in);
                    mom_def_Re_tr_paths[131] += ce * creal(path_storage[131 + idx]);
                    mom_def_Im_tr_paths[131] += I * ce * cimag(path_storage[131 + idx]);
                    path_storage[132 + idx] = path132(in);
                    mom_def_Re_tr_paths[132] += ce * creal(path_storage[132 + idx]);
                    mom_def_Im_tr_paths[132] += I * ce * cimag(path_storage[132 + idx]);
                    path_storage[133 + idx] = path133(in);
                    mom_def_Re_tr_paths[133] += ce * creal(path_storage[133 + idx]);
                    mom_def_Im_tr_paths[133] += I * ce * cimag(path_storage[133 + idx]);
                    path_storage[134 + idx] = path134(in);
                    mom_def_Re_tr_paths[134] += ce * creal(path_storage[134 + idx]);
                    mom_def_Im_tr_paths[134] += I * ce * cimag(path_storage[134 + idx]);
                    path_storage[135 + idx] = path135(in);
                    mom_def_Re_tr_paths[135] += ce * creal(path_storage[135 + idx]);
                    mom_def_Im_tr_paths[135] += I * ce * cimag(path_storage[135 + idx]);
                    path_storage[136 + idx] = path136(in);
                    mom_def_Re_tr_paths[136] += ce * creal(path_storage[136 + idx]);
                    mom_def_Im_tr_paths[136] += I * ce * cimag(path_storage[136 + idx]);
                    path_storage[137 + idx] = path137(in);
                    mom_def_Re_tr_paths[137] += ce * creal(path_storage[137 + idx]);
                    mom_def_Im_tr_paths[137] += I * ce * cimag(path_storage[137 + idx]);
                    path_storage[138 + idx] = path138(in);
                    mom_def_Re_tr_paths[138] += ce * creal(path_storage[138 + idx]);
                    mom_def_Im_tr_paths[138] += I * ce * cimag(path_storage[138 + idx]);
                    path_storage[139 + idx] = path139(in);
                    mom_def_Re_tr_paths[139] += ce * creal(path_storage[139 + idx]);
                    mom_def_Im_tr_paths[139] += I * ce * cimag(path_storage[139 + idx]);
                    path_storage[140 + idx] = path140(in);
                    mom_def_Re_tr_paths[140] += ce * creal(path_storage[140 + idx]);
                    mom_def_Im_tr_paths[140] += I * ce * cimag(path_storage[140 + idx]);
                    path_storage[141 + idx] = path141(in);
                    mom_def_Re_tr_paths[141] += ce * creal(path_storage[141 + idx]);
                    mom_def_Im_tr_paths[141] += I * ce * cimag(path_storage[141 + idx]);
                    path_storage[142 + idx] = path142(in);
                    mom_def_Re_tr_paths[142] += ce * creal(path_storage[142 + idx]);
                    mom_def_Im_tr_paths[142] += I * ce * cimag(path_storage[142 + idx]);
                    path_storage[143 + idx] = path143(in);
                    mom_def_Re_tr_paths[143] += ce * creal(path_storage[143 + idx]);
                    mom_def_Im_tr_paths[143] += I * ce * cimag(path_storage[143 + idx]);
                    path_storage[144 + idx] = path144(in);
                    mom_def_Re_tr_paths[144] += ce * creal(path_storage[144 + idx]);
                    mom_def_Im_tr_paths[144] += I * ce * cimag(path_storage[144 + idx]);
                    path_storage[145 + idx] = path145(in);
                    mom_def_Re_tr_paths[145] += ce * creal(path_storage[145 + idx]);
                    mom_def_Im_tr_paths[145] += I * ce * cimag(path_storage[145 + idx]);
                    path_storage[146 + idx] = path146(in);
                    mom_def_Re_tr_paths[146] += ce * creal(path_storage[146 + idx]);
                    mom_def_Im_tr_paths[146] += I * ce * cimag(path_storage[146 + idx]);
                    path_storage[147 + idx] = path147(in);
                    mom_def_Re_tr_paths[147] += ce * creal(path_storage[147 + idx]);
                    mom_def_Im_tr_paths[147] += I * ce * cimag(path_storage[147 + idx]);
                    path_storage[148 + idx] = path148(in);
                    mom_def_Re_tr_paths[148] += ce * creal(path_storage[148 + idx]);
                    mom_def_Im_tr_paths[148] += I * ce * cimag(path_storage[148 + idx]);
                    path_storage[149 + idx] = path149(in);
                    mom_def_Re_tr_paths[149] += ce * creal(path_storage[149 + idx]);
                    mom_def_Im_tr_paths[149] += I * ce * cimag(path_storage[149 + idx]);
                    path_storage[150 + idx] = path150(in);
                    mom_def_Re_tr_paths[150] += ce * creal(path_storage[150 + idx]);
                    mom_def_Im_tr_paths[150] += I * ce * cimag(path_storage[150 + idx]);
                    path_storage[151 + idx] = path151(in);
                    mom_def_Re_tr_paths[151] += ce * creal(path_storage[151 + idx]);
                    mom_def_Im_tr_paths[151] += I * ce * cimag(path_storage[151 + idx]);
                    path_storage[152 + idx] = path152(in);
                    mom_def_Re_tr_paths[152] += ce * creal(path_storage[152 + idx]);
                    mom_def_Im_tr_paths[152] += I * ce * cimag(path_storage[152 + idx]);
                    path_storage[153 + idx] = path153(in);
                    mom_def_Re_tr_paths[153] += ce * creal(path_storage[153 + idx]);
                    mom_def_Im_tr_paths[153] += I * ce * cimag(path_storage[153 + idx]);
                    path_storage[154 + idx] = path154(in);
                    mom_def_Re_tr_paths[154] += ce * creal(path_storage[154 + idx]);
                    mom_def_Im_tr_paths[154] += I * ce * cimag(path_storage[154 + idx]);
                    path_storage[155 + idx] = path155(in);
                    mom_def_Re_tr_paths[155] += ce * creal(path_storage[155 + idx]);
                    mom_def_Im_tr_paths[155] += I * ce * cimag(path_storage[155 + idx]);
                    path_storage[156 + idx] = path156(in);
                    mom_def_Re_tr_paths[156] += ce * creal(path_storage[156 + idx]);
                    mom_def_Im_tr_paths[156] += I * ce * cimag(path_storage[156 + idx]);
                    path_storage[157 + idx] = path157(in);
                    mom_def_Re_tr_paths[157] += ce * creal(path_storage[157 + idx]);
                    mom_def_Im_tr_paths[157] += I * ce * cimag(path_storage[157 + idx]);
                    path_storage[158 + idx] = path158(in);
                    mom_def_Re_tr_paths[158] += ce * creal(path_storage[158 + idx]);
                    mom_def_Im_tr_paths[158] += I * ce * cimag(path_storage[158 + idx]);
                    path_storage[159 + idx] = path159(in);
                    mom_def_Re_tr_paths[159] += ce * creal(path_storage[159 + idx]);
                    mom_def_Im_tr_paths[159] += I * ce * cimag(path_storage[159 + idx]);
                    path_storage[160 + idx] = path160(in);
                    mom_def_Re_tr_paths[160] += ce * creal(path_storage[160 + idx]);
                    mom_def_Im_tr_paths[160] += I * ce * cimag(path_storage[160 + idx]);
                    path_storage[161 + idx] = path161(in);
                    mom_def_Re_tr_paths[161] += ce * creal(path_storage[161 + idx]);
                    mom_def_Im_tr_paths[161] += I * ce * cimag(path_storage[161 + idx]);
                    path_storage[162 + idx] = path162(in);
                    mom_def_Re_tr_paths[162] += ce * creal(path_storage[162 + idx]);
                    mom_def_Im_tr_paths[162] += I * ce * cimag(path_storage[162 + idx]);
                    path_storage[163 + idx] = path163(in);
                    mom_def_Re_tr_paths[163] += ce * creal(path_storage[163 + idx]);
                    mom_def_Im_tr_paths[163] += I * ce * cimag(path_storage[163 + idx]);
                    path_storage[164 + idx] = path164(in);
                    mom_def_Re_tr_paths[164] += ce * creal(path_storage[164 + idx]);
                    mom_def_Im_tr_paths[164] += I * ce * cimag(path_storage[164 + idx]);
                    path_storage[165 + idx] = path165(in);
                    mom_def_Re_tr_paths[165] += ce * creal(path_storage[165 + idx]);
                    mom_def_Im_tr_paths[165] += I * ce * cimag(path_storage[165 + idx]);
                    path_storage[166 + idx] = path166(in);
                    mom_def_Re_tr_paths[166] += ce * creal(path_storage[166 + idx]);
                    mom_def_Im_tr_paths[166] += I * ce * cimag(path_storage[166 + idx]);
                    path_storage[167 + idx] = path167(in);
                    mom_def_Re_tr_paths[167] += ce * creal(path_storage[167 + idx]);
                    mom_def_Im_tr_paths[167] += I * ce * cimag(path_storage[167 + idx]);
                    path_storage[168 + idx] = path168(in);
                    mom_def_Re_tr_paths[168] += ce * creal(path_storage[168 + idx]);
                    mom_def_Im_tr_paths[168] += I * ce * cimag(path_storage[168 + idx]);
                    path_storage[169 + idx] = path169(in);
                    mom_def_Re_tr_paths[169] += ce * creal(path_storage[169 + idx]);
                    mom_def_Im_tr_paths[169] += I * ce * cimag(path_storage[169 + idx]);
                    path_storage[170 + idx] = path170(in);
                    mom_def_Re_tr_paths[170] += ce * creal(path_storage[170 + idx]);
                    mom_def_Im_tr_paths[170] += I * ce * cimag(path_storage[170 + idx]);
                    path_storage[171 + idx] = path171(in);
                    mom_def_Re_tr_paths[171] += ce * creal(path_storage[171 + idx]);
                    mom_def_Im_tr_paths[171] += I * ce * cimag(path_storage[171 + idx]);
                    path_storage[172 + idx] = path172(in);
                    mom_def_Re_tr_paths[172] += ce * creal(path_storage[172 + idx]);
                    mom_def_Im_tr_paths[172] += I * ce * cimag(path_storage[172 + idx]);
                    path_storage[173 + idx] = path173(in);
                    mom_def_Re_tr_paths[173] += ce * creal(path_storage[173 + idx]);
                    mom_def_Im_tr_paths[173] += I * ce * cimag(path_storage[173 + idx]);
                    path_storage[174 + idx] = path174(in);
                    mom_def_Re_tr_paths[174] += ce * creal(path_storage[174 + idx]);
                    mom_def_Im_tr_paths[174] += I * ce * cimag(path_storage[174 + idx]);
                    path_storage[175 + idx] = path175(in);
                    mom_def_Re_tr_paths[175] += ce * creal(path_storage[175 + idx]);
                    mom_def_Im_tr_paths[175] += I * ce * cimag(path_storage[175 + idx]);
                    path_storage[176 + idx] = path176(in);
                    mom_def_Re_tr_paths[176] += ce * creal(path_storage[176 + idx]);
                    mom_def_Im_tr_paths[176] += I * ce * cimag(path_storage[176 + idx]);
                    path_storage[177 + idx] = path177(in);
                    mom_def_Re_tr_paths[177] += ce * creal(path_storage[177 + idx]);
                    mom_def_Im_tr_paths[177] += I * ce * cimag(path_storage[177 + idx]);
                    path_storage[178 + idx] = path178(in);
                    mom_def_Re_tr_paths[178] += ce * creal(path_storage[178 + idx]);
                    mom_def_Im_tr_paths[178] += I * ce * cimag(path_storage[178 + idx]);
                    path_storage[179 + idx] = path179(in);
                    mom_def_Re_tr_paths[179] += ce * creal(path_storage[179 + idx]);
                    mom_def_Im_tr_paths[179] += I * ce * cimag(path_storage[179 + idx]);
                    path_storage[180 + idx] = path180(in);
                    mom_def_Re_tr_paths[180] += ce * creal(path_storage[180 + idx]);
                    mom_def_Im_tr_paths[180] += I * ce * cimag(path_storage[180 + idx]);
                    path_storage[181 + idx] = path181(in);
                    mom_def_Re_tr_paths[181] += ce * creal(path_storage[181 + idx]);
                    mom_def_Im_tr_paths[181] += I * ce * cimag(path_storage[181 + idx]);
                    path_storage[182 + idx] = path182(in);
                    mom_def_Re_tr_paths[182] += ce * creal(path_storage[182 + idx]);
                    mom_def_Im_tr_paths[182] += I * ce * cimag(path_storage[182 + idx]);
                    path_storage[183 + idx] = path183(in);
                    mom_def_Re_tr_paths[183] += ce * creal(path_storage[183 + idx]);
                    mom_def_Im_tr_paths[183] += I * ce * cimag(path_storage[183 + idx]);
                    path_storage[184 + idx] = path184(in);
                    mom_def_Re_tr_paths[184] += ce * creal(path_storage[184 + idx]);
                    mom_def_Im_tr_paths[184] += I * ce * cimag(path_storage[184 + idx]);
                    path_storage[185 + idx] = path185(in);
                    mom_def_Re_tr_paths[185] += ce * creal(path_storage[185 + idx]);
                    mom_def_Im_tr_paths[185] += I * ce * cimag(path_storage[185 + idx]);
                    path_storage[186 + idx] = path186(in);
                    mom_def_Re_tr_paths[186] += ce * creal(path_storage[186 + idx]);
                    mom_def_Im_tr_paths[186] += I * ce * cimag(path_storage[186 + idx]);
                    path_storage[187 + idx] = path187(in);
                    mom_def_Re_tr_paths[187] += ce * creal(path_storage[187 + idx]);
                    mom_def_Im_tr_paths[187] += I * ce * cimag(path_storage[187 + idx]);
                    path_storage[188 + idx] = path188(in);
                    mom_def_Re_tr_paths[188] += ce * creal(path_storage[188 + idx]);
                    mom_def_Im_tr_paths[188] += I * ce * cimag(path_storage[188 + idx]);
                    path_storage[189 + idx] = path189(in);
                    mom_def_Re_tr_paths[189] += ce * creal(path_storage[189 + idx]);
                    mom_def_Im_tr_paths[189] += I * ce * cimag(path_storage[189 + idx]);
                    path_storage[190 + idx] = path190(in);
                    mom_def_Re_tr_paths[190] += ce * creal(path_storage[190 + idx]);
                    mom_def_Im_tr_paths[190] += I * ce * cimag(path_storage[190 + idx]);
                    path_storage[191 + idx] = path191(in);
                    mom_def_Re_tr_paths[191] += ce * creal(path_storage[191 + idx]);
                    mom_def_Im_tr_paths[191] += I * ce * cimag(path_storage[191 + idx]);
                    path_storage[192 + idx] = path192(in);
                    mom_def_Re_tr_paths[192] += ce * creal(path_storage[192 + idx]);
                    mom_def_Im_tr_paths[192] += I * ce * cimag(path_storage[192 + idx]);
                    path_storage[193 + idx] = path193(in);
                    mom_def_Re_tr_paths[193] += ce * creal(path_storage[193 + idx]);
                    mom_def_Im_tr_paths[193] += I * ce * cimag(path_storage[193 + idx]);
                    path_storage[194 + idx] = path194(in);
                    mom_def_Re_tr_paths[194] += ce * creal(path_storage[194 + idx]);
                    mom_def_Im_tr_paths[194] += I * ce * cimag(path_storage[194 + idx]);
                    path_storage[195 + idx] = path195(in);
                    mom_def_Re_tr_paths[195] += ce * creal(path_storage[195 + idx]);
                    mom_def_Im_tr_paths[195] += I * ce * cimag(path_storage[195 + idx]);
                    path_storage[196 + idx] = path196(in);
                    mom_def_Re_tr_paths[196] += ce * creal(path_storage[196 + idx]);
                    mom_def_Im_tr_paths[196] += I * ce * cimag(path_storage[196 + idx]);
                    path_storage[197 + idx] = path197(in);
                    mom_def_Re_tr_paths[197] += ce * creal(path_storage[197 + idx]);
                    mom_def_Im_tr_paths[197] += I * ce * cimag(path_storage[197 + idx]);
                    path_storage[198 + idx] = path198(in);
                    mom_def_Re_tr_paths[198] += ce * creal(path_storage[198 + idx]);
                    mom_def_Im_tr_paths[198] += I * ce * cimag(path_storage[198 + idx]);
                    path_storage[199 + idx] = path199(in);
                    mom_def_Re_tr_paths[199] += ce * creal(path_storage[199 + idx]);
                    mom_def_Im_tr_paths[199] += I * ce * cimag(path_storage[199 + idx]);
                    path_storage[200 + idx] = path200(in);
                    mom_def_Re_tr_paths[200] += ce * creal(path_storage[200 + idx]);
                    mom_def_Im_tr_paths[200] += I * ce * cimag(path_storage[200 + idx]);
                    path_storage[201 + idx] = path201(in);
                    mom_def_Re_tr_paths[201] += ce * creal(path_storage[201 + idx]);
                    mom_def_Im_tr_paths[201] += I * ce * cimag(path_storage[201 + idx]);
                    path_storage[202 + idx] = path202(in);
                    mom_def_Re_tr_paths[202] += ce * creal(path_storage[202 + idx]);
                    mom_def_Im_tr_paths[202] += I * ce * cimag(path_storage[202 + idx]);
                    path_storage[203 + idx] = path203(in);
                    mom_def_Re_tr_paths[203] += ce * creal(path_storage[203 + idx]);
                    mom_def_Im_tr_paths[203] += I * ce * cimag(path_storage[203 + idx]);
                    path_storage[204 + idx] = path204(in);
                    mom_def_Re_tr_paths[204] += ce * creal(path_storage[204 + idx]);
                    mom_def_Im_tr_paths[204] += I * ce * cimag(path_storage[204 + idx]);
                    path_storage[205 + idx] = path205(in);
                    mom_def_Re_tr_paths[205] += ce * creal(path_storage[205 + idx]);
                    mom_def_Im_tr_paths[205] += I * ce * cimag(path_storage[205 + idx]);
                    path_storage[206 + idx] = path206(in);
                    mom_def_Re_tr_paths[206] += ce * creal(path_storage[206 + idx]);
                    mom_def_Im_tr_paths[206] += I * ce * cimag(path_storage[206 + idx]);
                    path_storage[207 + idx] = path207(in);
                    mom_def_Re_tr_paths[207] += ce * creal(path_storage[207 + idx]);
                    mom_def_Im_tr_paths[207] += I * ce * cimag(path_storage[207 + idx]);
                    path_storage[208 + idx] = path208(in);
                    mom_def_Re_tr_paths[208] += ce * creal(path_storage[208 + idx]);
                    mom_def_Im_tr_paths[208] += I * ce * cimag(path_storage[208 + idx]);
                    path_storage[209 + idx] = path209(in);
                    mom_def_Re_tr_paths[209] += ce * creal(path_storage[209 + idx]);
                    mom_def_Im_tr_paths[209] += I * ce * cimag(path_storage[209 + idx]);
                    path_storage[210 + idx] = path210(in);
                    mom_def_Re_tr_paths[210] += ce * creal(path_storage[210 + idx]);
                    mom_def_Im_tr_paths[210] += I * ce * cimag(path_storage[210 + idx]);
                    path_storage[211 + idx] = path211(in);
                    mom_def_Re_tr_paths[211] += ce * creal(path_storage[211 + idx]);
                    mom_def_Im_tr_paths[211] += I * ce * cimag(path_storage[211 + idx]);
                    path_storage[212 + idx] = path212(in);
                    mom_def_Re_tr_paths[212] += ce * creal(path_storage[212 + idx]);
                    mom_def_Im_tr_paths[212] += I * ce * cimag(path_storage[212 + idx]);
                    path_storage[213 + idx] = path213(in);
                    mom_def_Re_tr_paths[213] += ce * creal(path_storage[213 + idx]);
                    mom_def_Im_tr_paths[213] += I * ce * cimag(path_storage[213 + idx]);
                    path_storage[214 + idx] = path214(in);
                    mom_def_Re_tr_paths[214] += ce * creal(path_storage[214 + idx]);
                    mom_def_Im_tr_paths[214] += I * ce * cimag(path_storage[214 + idx]);
                    path_storage[215 + idx] = path215(in);
                    mom_def_Re_tr_paths[215] += ce * creal(path_storage[215 + idx]);
                    mom_def_Im_tr_paths[215] += I * ce * cimag(path_storage[215 + idx]);
                    path_storage[216 + idx] = path216(in);
                    mom_def_Re_tr_paths[216] += ce * creal(path_storage[216 + idx]);
                    mom_def_Im_tr_paths[216] += I * ce * cimag(path_storage[216 + idx]);
                    path_storage[217 + idx] = path217(in);
                    mom_def_Re_tr_paths[217] += ce * creal(path_storage[217 + idx]);
                    mom_def_Im_tr_paths[217] += I * ce * cimag(path_storage[217 + idx]);
                    path_storage[218 + idx] = path218(in);
                    mom_def_Re_tr_paths[218] += ce * creal(path_storage[218 + idx]);
                    mom_def_Im_tr_paths[218] += I * ce * cimag(path_storage[218 + idx]);
                    path_storage[219 + idx] = path219(in);
                    mom_def_Re_tr_paths[219] += ce * creal(path_storage[219 + idx]);
                    mom_def_Im_tr_paths[219] += I * ce * cimag(path_storage[219 + idx]);
                    path_storage[220 + idx] = path220(in);
                    mom_def_Re_tr_paths[220] += ce * creal(path_storage[220 + idx]);
                    mom_def_Im_tr_paths[220] += I * ce * cimag(path_storage[220 + idx]);
                    path_storage[221 + idx] = path221(in);
                    mom_def_Re_tr_paths[221] += ce * creal(path_storage[221 + idx]);
                    mom_def_Im_tr_paths[221] += I * ce * cimag(path_storage[221 + idx]);
                    path_storage[222 + idx] = path222(in);
                    mom_def_Re_tr_paths[222] += ce * creal(path_storage[222 + idx]);
                    mom_def_Im_tr_paths[222] += I * ce * cimag(path_storage[222 + idx]);
                    path_storage[223 + idx] = path223(in);
                    mom_def_Re_tr_paths[223] += ce * creal(path_storage[223 + idx]);
                    mom_def_Im_tr_paths[223] += I * ce * cimag(path_storage[223 + idx]);
                    path_storage[224 + idx] = path224(in);
                    mom_def_Re_tr_paths[224] += ce * creal(path_storage[224 + idx]);
                    mom_def_Im_tr_paths[224] += I * ce * cimag(path_storage[224 + idx]);
                    path_storage[225 + idx] = path225(in);
                    mom_def_Re_tr_paths[225] += ce * creal(path_storage[225 + idx]);
                    mom_def_Im_tr_paths[225] += I * ce * cimag(path_storage[225 + idx]);
                    path_storage[226 + idx] = path226(in);
                    mom_def_Re_tr_paths[226] += ce * creal(path_storage[226 + idx]);
                    mom_def_Im_tr_paths[226] += I * ce * cimag(path_storage[226 + idx]);
                    path_storage[227 + idx] = path227(in);
                    mom_def_Re_tr_paths[227] += ce * creal(path_storage[227 + idx]);
                    mom_def_Im_tr_paths[227] += I * ce * cimag(path_storage[227 + idx]);
                    path_storage[228 + idx] = path228(in);
                    mom_def_Re_tr_paths[228] += ce * creal(path_storage[228 + idx]);
                    mom_def_Im_tr_paths[228] += I * ce * cimag(path_storage[228 + idx]);
                }
    }
    else
    {
        for (n_x = 0; n_x < X; n_x++)
            for (n_y = 0; n_y < Y; n_y++)
                for (n_z = 0; n_z < Z; n_z++)
                {
                    ce = cexp(I * 2.0 * PI * (double)(n_x * px + n_y * py + n_z * pz) / GLB_X);
                    idx = npaths * (n_x + X * (n_y + Y * n_z));
                    for (int i = 0; i < 229; i++)
                    {
                        mom_def_Re_tr_paths[i] += ce * creal(path_storage[i + idx]);
                        mom_def_Im_tr_paths[i] += I * ce * cimag(path_storage[i + idx]);
                    }
                }
    }
};
void eval_all_glueball_ops(int t, double complex *numerical_op_out)
{
    static double complex *numerical_op = NULL;
    if (numerical_op == NULL)
    {
        numerical_op = malloc(total_n_glue_op * sizeof(double complex));
    }
    request_space_paths_evaluation();
    eval_time_momentum_glueball_paths(t, -1, 0, 1);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_5(numerical_op + 4);
    eval_time_momentum_glueball_paths(t, -1, 1, 0);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_1(numerical_op + 5);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_2(numerical_op + 6);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_3(numerical_op + 7);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_4(numerical_op + 8);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_5(numerical_op + 9);
    eval_time_momentum_glueball_paths(t, 0, -1, 1);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_1(numerical_op + 10);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_2(numerical_op + 11);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_3(numerical_op + 12);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_4(numerical_op + 13);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_5(numerical_op + 14);
    eval_time_momentum_glueball_paths(t, 0, 0, 0);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_1(numerical_op + 15);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_2(numerical_op + 16);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_3(numerical_op + 17);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_4(numerical_op + 18);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_5(numerical_op + 19);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_6(numerical_op + 20);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_7(numerical_op + 21);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_8(numerical_op + 22);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_9(numerical_op + 23);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_10(numerical_op + 24);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_11(numerical_op + 25);
    eval_time_momentum_glueball_paths(t, 0, 0, 1);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_1(numerical_op + 36);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_2(numerical_op + 37);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_3(numerical_op + 38);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_4(numerical_op + 39);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_5(numerical_op + 40);
    eval_time_momentum_glueball_paths(t, 0, 1, -1);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_1(numerical_op + 41);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_2(numerical_op + 42);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_3(numerical_op + 43);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_4(numerical_op + 44);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_5(numerical_op + 45);
    eval_time_momentum_glueball_paths(t, 0, 1, 0);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_1(numerical_op + 46);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_2(numerical_op + 47);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_3(numerical_op + 48);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_4(numerical_op + 49);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_5(numerical_op + 50);
    eval_time_momentum_glueball_paths(t, 0, 1, 1);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_1(numerical_op + 51);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_2(numerical_op + 52);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_3(numerical_op + 53);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_4(numerical_op + 54);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_5(numerical_op + 55);
    eval_time_momentum_glueball_paths(t, 1, -1, 0);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_1(numerical_op + 56);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_2(numerical_op + 57);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_3(numerical_op + 58);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_4(numerical_op + 59);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_5(numerical_op + 60);
    eval_time_momentum_glueball_paths(t, 1, 0, -1);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_1(numerical_op + 61);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_2(numerical_op + 62);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_3(numerical_op + 63);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_4(numerical_op + 64);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_5(numerical_op + 65);
    eval_time_momentum_glueball_paths(t, 1, 0, 0);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_1(numerical_op + 66);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_2(numerical_op + 67);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_3(numerical_op + 68);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_4(numerical_op + 69);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_5(numerical_op + 70);
    eval_time_momentum_glueball_paths(t, 1, 0, 1);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_1(numerical_op + 71);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_2(numerical_op + 72);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_3(numerical_op + 73);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_4(numerical_op + 74);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_5(numerical_op + 75);
    eval_time_momentum_glueball_paths(t, 1, 1, 0);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_1(numerical_op + 76);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_2(numerical_op + 77);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_3(numerical_op + 78);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_4(numerical_op + 79);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_5(numerical_op + 80);
    *(numerical_op + 26) = +(0.816496580927726033) * (conj(*(numerical_op + 36))) * ((*(numerical_op + 36))) + (0.816496580927726033) * (conj(*(numerical_op + 46))) * ((*(numerical_op + 46))) + (0.816496580927726033) * (conj(*(numerical_op + 66))) * ((*(numerical_op + 66)));
    *(numerical_op + 27) = +(0.816496580927726033) * (conj(*(numerical_op + 37))) * ((*(numerical_op + 37))) + (0.816496580927726033) * (conj(*(numerical_op + 47))) * ((*(numerical_op + 47))) + (0.816496580927726033) * (conj(*(numerical_op + 67))) * ((*(numerical_op + 67)));
    *(numerical_op + 28) = +(0.816496580927726033) * (conj(*(numerical_op + 38))) * ((*(numerical_op + 38))) + (0.816496580927726033) * (conj(*(numerical_op + 48))) * ((*(numerical_op + 48))) + (0.816496580927726033) * (conj(*(numerical_op + 68))) * ((*(numerical_op + 68)));
    *(numerical_op + 29) = +(0.816496580927726033) * (conj(*(numerical_op + 39))) * ((*(numerical_op + 39))) + (0.816496580927726033) * (conj(*(numerical_op + 49))) * ((*(numerical_op + 49))) + (0.816496580927726033) * (conj(*(numerical_op + 69))) * ((*(numerical_op + 69)));
    *(numerical_op + 30) = +(0.816496580927726033) * (conj(*(numerical_op + 40))) * ((*(numerical_op + 40))) + (0.816496580927726033) * (conj(*(numerical_op + 50))) * ((*(numerical_op + 50))) + (0.816496580927726033) * (conj(*(numerical_op + 70))) * ((*(numerical_op + 70)));
    *(numerical_op + 31) = +(0.577350269189625765) * ((*(numerical_op + 10))) * ((*(numerical_op + 41))) + (0.577350269189625765) * (conj(*(numerical_op + 51))) * ((*(numerical_op + 51))) + (0.577350269189625765) * ((*(numerical_op + 5))) * ((*(numerical_op + 56))) + (0.577350269189625765) * ((*(numerical_op + 0))) * ((*(numerical_op + 61))) + (0.577350269189625765) * (conj(*(numerical_op + 71))) * ((*(numerical_op + 71))) + (0.577350269189625765) * (conj(*(numerical_op + 76))) * ((*(numerical_op + 76)));
    *(numerical_op + 32) = +(0.577350269189625765) * ((*(numerical_op + 11))) * ((*(numerical_op + 42))) + (0.577350269189625765) * (conj(*(numerical_op + 52))) * ((*(numerical_op + 52))) + (0.577350269189625765) * ((*(numerical_op + 6))) * ((*(numerical_op + 57))) + (0.577350269189625765) * ((*(numerical_op + 1))) * ((*(numerical_op + 62))) + (0.577350269189625765) * (conj(*(numerical_op + 72))) * ((*(numerical_op + 72))) + (0.577350269189625765) * (conj(*(numerical_op + 77))) * ((*(numerical_op + 77)));
    *(numerical_op + 33) = +(0.577350269189625765) * ((*(numerical_op + 12))) * ((*(numerical_op + 43))) + (0.577350269189625765) * (conj(*(numerical_op + 53))) * ((*(numerical_op + 53))) + (0.577350269189625765) * ((*(numerical_op + 7))) * ((*(numerical_op + 58))) + (0.577350269189625765) * ((*(numerical_op + 2))) * ((*(numerical_op + 63))) + (0.577350269189625765) * (conj(*(numerical_op + 73))) * ((*(numerical_op + 73))) + (0.577350269189625765) * (conj(*(numerical_op + 78))) * ((*(numerical_op + 78)));
    *(numerical_op + 34) = +(0.577350269189625765) * ((*(numerical_op + 13))) * ((*(numerical_op + 44))) + (0.577350269189625765) * (conj(*(numerical_op + 54))) * ((*(numerical_op + 54))) + (0.577350269189625765) * ((*(numerical_op + 8))) * ((*(numerical_op + 59))) + (0.577350269189625765) * ((*(numerical_op + 3))) * ((*(numerical_op + 64))) + (0.577350269189625765) * (conj(*(numerical_op + 74))) * ((*(numerical_op + 74))) + (0.577350269189625765) * (conj(*(numerical_op + 79))) * ((*(numerical_op + 79)));
    *(numerical_op + 35) = +(0.577350269189625765) * ((*(numerical_op + 14))) * ((*(numerical_op + 45))) + (0.577350269189625765) * (conj(*(numerical_op + 55))) * ((*(numerical_op + 55))) + (0.577350269189625765) * ((*(numerical_op + 9))) * ((*(numerical_op + 60))) + (0.577350269189625765) * ((*(numerical_op + 4))) * ((*(numerical_op + 65))) + (0.577350269189625765) * (conj(*(numerical_op + 75))) * ((*(numerical_op + 75))) + (0.577350269189625765) * (conj(*(numerical_op + 80))) * ((*(numerical_op + 80)));
    for (int i = 0; i < total_n_glue_op; i++)
        *(numerical_op_out + i) += *(numerical_op + i);
}
void eval_all_glueball_oneTr_ops_p_m1_0_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_oneTr_ops_p_m1_0_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_m1_1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_oneTr_ops_p_m1_1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_0_m1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_oneTr_ops_p_0_m1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_0_0_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_1_C_m1_n_3(numerical_op + 2);
}

void eval_all_glueball_oneTr_ops_p_0_0_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_8(numerical_op + 7);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_9(numerical_op + 8);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_10(numerical_op + 9);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_11(numerical_op + 10);
}

void eval_all_glueball_oneTr_ops_p_0_0_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_0_1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_oneTr_ops_p_0_1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_0_1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_oneTr_ops_p_0_1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_0_1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_3(numerical_op + 2);
}

void eval_all_glueball_oneTr_ops_p_0_1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_1_m1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_oneTr_ops_p_1_m1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_1_0_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_oneTr_ops_p_1_0_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_1_0_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_1_0_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_5(numerical_op + 4);
}

void eval_all_glueball_oneTr_ops_p_1_1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_0_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_0_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_0_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_oneTr_ops_p_1_1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_5(numerical_op + 4);
}

void evaluate_correlators(cor_list *lcor, int nblocking, double complex *gb_storage, double *cor_storage)
{
    int totalsize, i1, i2, t1, t2, id, n1, n2, b1, b2, i;
    double norm;
    double *cor_pointer, tmp;
    static double complex *gb1_bf;
    static int n_total_active_slices = 0;
    static int *listactive = NULL;
    if (listactive == NULL)
    {
        listactive = malloc(sizeof(int) * GLB_T);
        for (i = 0; i < GLB_T; i++)
            listactive[i] = -1;

        for (i = 0; i < lcor->n_entries; i++)
        {
            listactive[lcor->list[i].t1] = 1;
            listactive[lcor->list[i].t2] = 1;
        }

        for (i = 0; i < GLB_T; i++)
            if (listactive[i] == 1)
            {
                listactive[i] = n_total_active_slices;
                n_total_active_slices++;
            }
    }

#ifdef WITH_MPI
    static int *t_to_proc = NULL;
    static int *listsent = NULL;

    if (t_to_proc == NULL)
    {
        listsent = malloc(sizeof(int) * GLB_T);

        gb1_bf = malloc(sizeof(double complex) * total_n_glue_op * nblocking * n_total_active_slices);

        t_to_proc = malloc(sizeof(int) * GLB_T);
        for (i = 0; i < GLB_T; i++)
        {
            int x[4] = {i, 0, 0, 0}, p[4];
            glb_to_proc(x, p);
            t_to_proc[i] = p[0];
        }
    }

    for (i = 0; i < GLB_T; i++)
        listsent[i] = -1;

    MPI_Status r1, r2;
#else
    gb1_bf = gb_storage;
#endif

    static double complex *gb1;
    static double complex *gb2;

    for (int icor = 0; icor < lcor->n_entries; icor++)
    {
        t1 = glbT_to_active_slices[lcor->list[icor].t1];
        t2 = glbT_to_active_slices[lcor->list[icor].t2];
        id = lcor->list[icor].id;
        gb1 = gb1_bf + total_n_glue_op * nblocking * listactive[lcor->list[icor].t1];
        gb2 = gb1_bf + total_n_glue_op * nblocking * listactive[lcor->list[icor].t2];

#ifdef WITH_MPI

        if (listsent[lcor->list[icor].t1] == -1)
        {
            if (t1 != -1)
            {
                if (PID == 0)
                {
                    memcpy(gb1, gb_storage + t1 * total_n_glue_op * nblocking, sizeof(double complex) * total_n_glue_op * nblocking);
                }
                else
                {
                    MPI_Send(gb_storage + t1 * total_n_glue_op * nblocking, total_n_glue_op * nblocking * 2, MPI_DOUBLE, 0, lcor->list[icor].t1, cart_comm);
                }
            }

            if (PID == 0 && t1 == -1)
            {
                MPI_Recv(gb1, total_n_glue_op * nblocking * 2, MPI_DOUBLE, t_to_proc[lcor->list[icor].t1], lcor->list[icor].t1, cart_comm, &r1);
            }
            listsent[lcor->list[icor].t1] = 0;
        }

        if (lcor->list[icor].t1 != lcor->list[icor].t2)
        {
            if (listsent[lcor->list[icor].t2] == -1)
            {
                if (t2 != -1)
                {
                    if (PID == 0)
                    {
                        memcpy(gb2, gb_storage + t2 * total_n_glue_op * nblocking, sizeof(double complex) * total_n_glue_op * nblocking);
                    }
                    else
                    {
                        MPI_Send(gb_storage + t2 * total_n_glue_op * nblocking, total_n_glue_op * nblocking * 2, MPI_DOUBLE, 0, GLB_T + lcor->list[icor].t2, cart_comm);
                    }
                }

                if (PID == 0 && t2 == -1)
                {
                    MPI_Recv(gb2, total_n_glue_op * nblocking * 2, MPI_DOUBLE, t_to_proc[lcor->list[icor].t2], GLB_T + lcor->list[icor].t2, cart_comm, &r2);
                }
                listsent[lcor->list[icor].t2] = 0;
            }
        }
#endif

        totalsize = 0;
        norm = 0.5 / lcor->list[icor].n_pairs;

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 0 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 0 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 0 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 0 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 5 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 5 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 5 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 5 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 10 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 10 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 10 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 10 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 441;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 21; b1++)
            {
                i1 = 15 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 15 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 21 * (n1 + nblocking * (b2 + 21 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 21; b2++)
                {
                    i2 = 15 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 21 * (n1 + nblocking * (b2 + 21 * n2))] += tmp;
                    cor_pointer[b2 + 21 * (n2 + nblocking * (b1 + 21 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 21; b2++)
                    {
                        i2 = 15 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 21 * (n1 + nblocking * (b2 + 21 * n2))] += tmp;
                        cor_pointer[b2 + 21 * (n2 + nblocking * (b1 + 21 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 441 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 36 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 36 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 36 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 36 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 41 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 41 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 41 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 41 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 46 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 46 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 46 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 46 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 51 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 51 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 51 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 51 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 56 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 56 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 56 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 56 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 61 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 61 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 61 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 61 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 66 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 66 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 66 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 66 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 71 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 71 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 71 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 71 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));

        cor_pointer = cor_storage + totalsize + id * nblocking * nblocking * 25;

        for (n1 = 0; n1 < nblocking; n1++)
            for (b1 = 0; b1 < 5; b1++)
            {
                i1 = 76 + b1 + total_n_glue_op * n1;

                n2 = n1;
                b2 = b1;
                i2 = 76 + b2 + total_n_glue_op * n2;
                tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;

                for (b2 = b1 + 1; b2 < 5; b2++)
                {
                    i2 = 76 + b2 + total_n_glue_op * n2;
                    tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                    cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                    cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                }

                for (n2 = n1 + 1; n2 < nblocking; n2++)
                    for (b2 = 0; b2 < 5; b2++)
                    {
                        i2 = 76 + b2 + total_n_glue_op * n2;
                        tmp = norm * creal(conj(gb1[i1]) * gb2[i2] + gb1[i2] * conj(gb2[i1]));
                        cor_pointer[b1 + 5 * (n1 + nblocking * (b2 + 5 * n2))] += tmp;
                        cor_pointer[b2 + 5 * (n2 + nblocking * (b1 + 5 * n1))] += tmp;
                    }
            }
        totalsize += (nblocking * nblocking * 25 * (lcor->n_corrs));
    }

    totalsize = 0;

    lprintf("Measure ML", 0, "\n1pt function P=(-1,0,1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 0 1 2 3 4 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 0; i < 5; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(-1,1,0) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 5 6 7 8 9 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 5; i < 10; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,-1,1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 10 11 12 13 14 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 10; i < 15; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,0) Irrep=A1plusOhP Irrep ev=1/1 Charge=+ nop=%d\n", 21 * nblocking);
    lprintf("Measure ML", 0, "Op id= 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 15; i < 36; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,1) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 36 37 38 39 40 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 36; i < 41; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,1,-1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 41 42 43 44 45 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 41; i < 46; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,1,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 46 47 48 49 50 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 46; i < 51; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,1,1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 51 52 53 54 55 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 51; i < 56; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,-1,0) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 56 57 58 59 60 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 56; i < 61; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,0,-1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 61 62 63 64 65 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 61; i < 66; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,0,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 66 67 68 69 70 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 66; i < 71; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,0,1) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 71 72 73 74 75 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 71; i < 76; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,1,0) Irrep=A1Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 5 * nblocking);
    lprintf("Measure ML", 0, "Op id= 76 77 78 79 80 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 76; i < 81; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }
}
void report_op_group_setup()
{
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(-1,0,1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |0|1|2|3|4|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(-1,1,0) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |5|6|7|8|9|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,-1,1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |10|11|12|13|14|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,0,0) Irrep=A1plusOhP Charge=+");
    lprintf("INIT Measure ML", 0, " |15|16|17|18|19|20|21|22|23|24|25|2tr|2tr|2tr|2tr|2tr|2tr|2tr|2tr|2tr|2tr|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,0,1) Irrep=A1Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |36|37|38|39|40|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,1,-1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |41|42|43|44|45|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,1,0) Irrep=A1Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |46|47|48|49|50|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,1,1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |51|52|53|54|55|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,-1,0) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |56|57|58|59|60|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,0,-1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |61|62|63|64|65|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,0,0) Irrep=A1Dic4 Charge=+");
    lprintf("INIT Measure ML", 0, " |66|67|68|69|70|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,0,1) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |71|72|73|74|75|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,1,0) Irrep=A1Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |76|77|78|79|80|");
}
