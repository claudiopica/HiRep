#include <stdlib.h>
#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "utils.h"
#include "glueballs.h"

#include <string.h>
#define npaths 74
static double PI = 3.141592653589793238462643383279502884197;
static double complex *mom_def_Cp_tr_paths = NULL;
static double complex *mom_def_Cm_tr_paths = NULL;
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
static double complex path14(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path17(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path16(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path10(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path11(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path15(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path12(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path13(int in)
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex c0;
static void OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +c0 * mom_def_Cm_tr_paths[14] + c0 * mom_def_Cm_tr_paths[17] + c0 * mom_def_Cm_tr_paths[16] + mom_def_Cm_tr_paths[10] + c0 * mom_def_Cm_tr_paths[11] + c0 * mom_def_Cm_tr_paths[15] + mom_def_Cm_tr_paths[12] + c0 * mom_def_Cm_tr_paths[13];
}

static double complex path48(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path49(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path46(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

static double complex path23(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path22(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static double complex path47(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path36(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static double complex path37(int in)
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

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static void OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +c0 * mom_def_Cm_tr_paths[48] + c0 * mom_def_Cm_tr_paths[49] + c0 * mom_def_Cm_tr_paths[46] - mom_def_Cm_tr_paths[23] - c0 * mom_def_Cm_tr_paths[22] + mom_def_Cm_tr_paths[47] - c0 * mom_def_Cm_tr_paths[36] - mom_def_Cm_tr_paths[37];
}

static double complex path19(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

static double complex path20(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path21(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path18(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path7(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static double complex path9(int in)
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

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static void OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +c0 * mom_def_Cm_tr_paths[19] + c0 * mom_def_Cm_tr_paths[20] + c0 * mom_def_Cm_tr_paths[21] + c0 * mom_def_Cm_tr_paths[23] + mom_def_Cm_tr_paths[18] + c0 * mom_def_Cm_tr_paths[22] - c0 * mom_def_Cm_tr_paths[7] - mom_def_Cm_tr_paths[9];
}

static double complex path44(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path41(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path42(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path38(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path40(int in)
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +c0 * mom_def_Cm_tr_paths[44] + mom_def_Cm_tr_paths[43] + c0 * mom_def_Cm_tr_paths[39] + c0 * mom_def_Cm_tr_paths[41] + c0 * mom_def_Cm_tr_paths[45] + c0 * mom_def_Cm_tr_paths[42] + mom_def_Cm_tr_paths[38] + c0 * mom_def_Cm_tr_paths[40];
}

static double complex path34(int in)
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

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path32(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

static double complex path3(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

static double complex path5(int in)
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

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static void OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -mom_def_Cm_tr_paths[34] - c0 * mom_def_Cm_tr_paths[32] - mom_def_Cm_tr_paths[48] - c0 * mom_def_Cm_tr_paths[19] - mom_def_Cm_tr_paths[20] - c0 * mom_def_Cm_tr_paths[3] - mom_def_Cm_tr_paths[5] - c0 * mom_def_Cm_tr_paths[46];
}

static double complex path66(int in)
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path72(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path71(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path68(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path70(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path73(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +c0 * mom_def_Cm_tr_paths[66] + mom_def_Cm_tr_paths[72] + c0 * mom_def_Cm_tr_paths[67] + c0 * mom_def_Cm_tr_paths[71] + mom_def_Cm_tr_paths[68] + c0 * mom_def_Cm_tr_paths[69] + c0 * mom_def_Cm_tr_paths[70] + c0 * mom_def_Cm_tr_paths[73];
}

static double complex path58(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path60(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path64(int in)
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    _suNg_trace(p, res);
    return p;
}

static double complex path62(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Cm_tr_paths[58] + mom_def_Cm_tr_paths[61] + mom_def_Cm_tr_paths[63] + mom_def_Cm_tr_paths[60] + mom_def_Cm_tr_paths[64] + mom_def_Cm_tr_paths[65] + mom_def_Cm_tr_paths[59] + mom_def_Cm_tr_paths[62];
}

static double complex path35(int in)
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
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path6(int in)
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
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static double complex path33(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path4(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_0_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = -mom_def_Cm_tr_paths[35] - mom_def_Cm_tr_paths[49] - mom_def_Cm_tr_paths[21] - mom_def_Cm_tr_paths[6] - mom_def_Cm_tr_paths[18] - mom_def_Cm_tr_paths[33] - mom_def_Cm_tr_paths[47] - mom_def_Cm_tr_paths[4];
}

static double complex path26(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path29(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    _suNg_trace(p, res);
    return p;
}

static double complex path28(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path25(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path30(int in)
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path31(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Cm_tr_paths[26] + mom_def_Cm_tr_paths[27] + mom_def_Cm_tr_paths[29] + mom_def_Cm_tr_paths[24] + mom_def_Cm_tr_paths[28] + mom_def_Cm_tr_paths[25] + mom_def_Cm_tr_paths[30] + mom_def_Cm_tr_paths[31];
}

static double complex path8(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path2(int in)
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
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

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

static void OP_oneTr_p_0_1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +mom_def_Cm_tr_paths[34] + mom_def_Cm_tr_paths[35] + mom_def_Cm_tr_paths[32] - mom_def_Cm_tr_paths[8] + mom_def_Cm_tr_paths[33] + mom_def_Cm_tr_paths[36] + mom_def_Cm_tr_paths[37] - mom_def_Cm_tr_paths[2];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Cm_tr_paths[3] + mom_def_Cm_tr_paths[5] + mom_def_Cm_tr_paths[6] + mom_def_Cm_tr_paths[8] + mom_def_Cm_tr_paths[7] + mom_def_Cm_tr_paths[4] + mom_def_Cm_tr_paths[9] + mom_def_Cm_tr_paths[2];
}

static double complex path51(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path54(int in)
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path56(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path52(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path50(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_1_0_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +mom_def_Cm_tr_paths[51] + mom_def_Cm_tr_paths[55] + mom_def_Cm_tr_paths[53] + mom_def_Cm_tr_paths[54] + mom_def_Cm_tr_paths[56] + mom_def_Cm_tr_paths[52] + mom_def_Cm_tr_paths[57] + mom_def_Cm_tr_paths[50];
}

static double complex path0(int in)
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

static double complex path1(int in)
{
    suNg *w1, *w2;
    suNg res, res1;
    int site = in;
    double complex p;

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
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

static double complex c1;
static void OP_oneTr_p_1_1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = +mom_def_Cp_tr_paths[0] + mom_def_Cp_tr_paths[0] - c1 * mom_def_Cp_tr_paths[1] - c1 * mom_def_Cp_tr_paths[1];
}

static int last_t = -10;
void request_space_paths_evaluation() { last_t = -10; }
static void eval_time_momentum_glueball_paths(int t, int px, int py, int pz)
{
    int n_x, n_y, n_z, idx, in;
    double complex ce = 0.;
    if (path_storage == NULL)
    {
        c0 = cexp(I * PI * (-2. / GLB_X));
        c1 = cexp(I * PI * (2. / GLB_X));
        path_storage = malloc(npaths * X * Y * Z * sizeof(double complex));
        mom_def_Cp_tr_paths = malloc(npaths * sizeof(double complex));
        mom_def_Cm_tr_paths = malloc(npaths * sizeof(double complex));
        for (in = 0; in < npaths * X * Y * Z; in++)
            path_storage[in] = 0.;
    }
    for (in = 0; in < npaths; in++)
    {
        mom_def_Cp_tr_paths[in] = 0.;
        mom_def_Cm_tr_paths[in] = 0.;
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
                    mom_def_Cp_tr_paths[0] += ce * creal(path_storage[0 + idx]);
                    mom_def_Cm_tr_paths[0] += I * ce * cimag(path_storage[0 + idx]);
                    path_storage[1 + idx] = path1(in);
                    mom_def_Cp_tr_paths[1] += ce * creal(path_storage[1 + idx]);
                    mom_def_Cm_tr_paths[1] += I * ce * cimag(path_storage[1 + idx]);
                    path_storage[2 + idx] = path2(in);
                    mom_def_Cp_tr_paths[2] += ce * creal(path_storage[2 + idx]);
                    mom_def_Cm_tr_paths[2] += I * ce * cimag(path_storage[2 + idx]);
                    path_storage[3 + idx] = path3(in);
                    mom_def_Cp_tr_paths[3] += ce * creal(path_storage[3 + idx]);
                    mom_def_Cm_tr_paths[3] += I * ce * cimag(path_storage[3 + idx]);
                    path_storage[4 + idx] = path4(in);
                    mom_def_Cp_tr_paths[4] += ce * creal(path_storage[4 + idx]);
                    mom_def_Cm_tr_paths[4] += I * ce * cimag(path_storage[4 + idx]);
                    path_storage[5 + idx] = path5(in);
                    mom_def_Cp_tr_paths[5] += ce * creal(path_storage[5 + idx]);
                    mom_def_Cm_tr_paths[5] += I * ce * cimag(path_storage[5 + idx]);
                    path_storage[6 + idx] = path6(in);
                    mom_def_Cp_tr_paths[6] += ce * creal(path_storage[6 + idx]);
                    mom_def_Cm_tr_paths[6] += I * ce * cimag(path_storage[6 + idx]);
                    path_storage[7 + idx] = path7(in);
                    mom_def_Cp_tr_paths[7] += ce * creal(path_storage[7 + idx]);
                    mom_def_Cm_tr_paths[7] += I * ce * cimag(path_storage[7 + idx]);
                    path_storage[8 + idx] = path8(in);
                    mom_def_Cp_tr_paths[8] += ce * creal(path_storage[8 + idx]);
                    mom_def_Cm_tr_paths[8] += I * ce * cimag(path_storage[8 + idx]);
                    path_storage[9 + idx] = path9(in);
                    mom_def_Cp_tr_paths[9] += ce * creal(path_storage[9 + idx]);
                    mom_def_Cm_tr_paths[9] += I * ce * cimag(path_storage[9 + idx]);
                    path_storage[10 + idx] = path10(in);
                    mom_def_Cp_tr_paths[10] += ce * creal(path_storage[10 + idx]);
                    mom_def_Cm_tr_paths[10] += I * ce * cimag(path_storage[10 + idx]);
                    path_storage[11 + idx] = path11(in);
                    mom_def_Cp_tr_paths[11] += ce * creal(path_storage[11 + idx]);
                    mom_def_Cm_tr_paths[11] += I * ce * cimag(path_storage[11 + idx]);
                    path_storage[12 + idx] = path12(in);
                    mom_def_Cp_tr_paths[12] += ce * creal(path_storage[12 + idx]);
                    mom_def_Cm_tr_paths[12] += I * ce * cimag(path_storage[12 + idx]);
                    path_storage[13 + idx] = path13(in);
                    mom_def_Cp_tr_paths[13] += ce * creal(path_storage[13 + idx]);
                    mom_def_Cm_tr_paths[13] += I * ce * cimag(path_storage[13 + idx]);
                    path_storage[14 + idx] = path14(in);
                    mom_def_Cp_tr_paths[14] += ce * creal(path_storage[14 + idx]);
                    mom_def_Cm_tr_paths[14] += I * ce * cimag(path_storage[14 + idx]);
                    path_storage[15 + idx] = path15(in);
                    mom_def_Cp_tr_paths[15] += ce * creal(path_storage[15 + idx]);
                    mom_def_Cm_tr_paths[15] += I * ce * cimag(path_storage[15 + idx]);
                    path_storage[16 + idx] = path16(in);
                    mom_def_Cp_tr_paths[16] += ce * creal(path_storage[16 + idx]);
                    mom_def_Cm_tr_paths[16] += I * ce * cimag(path_storage[16 + idx]);
                    path_storage[17 + idx] = path17(in);
                    mom_def_Cp_tr_paths[17] += ce * creal(path_storage[17 + idx]);
                    mom_def_Cm_tr_paths[17] += I * ce * cimag(path_storage[17 + idx]);
                    path_storage[18 + idx] = path18(in);
                    mom_def_Cp_tr_paths[18] += ce * creal(path_storage[18 + idx]);
                    mom_def_Cm_tr_paths[18] += I * ce * cimag(path_storage[18 + idx]);
                    path_storage[19 + idx] = path19(in);
                    mom_def_Cp_tr_paths[19] += ce * creal(path_storage[19 + idx]);
                    mom_def_Cm_tr_paths[19] += I * ce * cimag(path_storage[19 + idx]);
                    path_storage[20 + idx] = path20(in);
                    mom_def_Cp_tr_paths[20] += ce * creal(path_storage[20 + idx]);
                    mom_def_Cm_tr_paths[20] += I * ce * cimag(path_storage[20 + idx]);
                    path_storage[21 + idx] = path21(in);
                    mom_def_Cp_tr_paths[21] += ce * creal(path_storage[21 + idx]);
                    mom_def_Cm_tr_paths[21] += I * ce * cimag(path_storage[21 + idx]);
                    path_storage[22 + idx] = path22(in);
                    mom_def_Cp_tr_paths[22] += ce * creal(path_storage[22 + idx]);
                    mom_def_Cm_tr_paths[22] += I * ce * cimag(path_storage[22 + idx]);
                    path_storage[23 + idx] = path23(in);
                    mom_def_Cp_tr_paths[23] += ce * creal(path_storage[23 + idx]);
                    mom_def_Cm_tr_paths[23] += I * ce * cimag(path_storage[23 + idx]);
                    path_storage[24 + idx] = path24(in);
                    mom_def_Cp_tr_paths[24] += ce * creal(path_storage[24 + idx]);
                    mom_def_Cm_tr_paths[24] += I * ce * cimag(path_storage[24 + idx]);
                    path_storage[25 + idx] = path25(in);
                    mom_def_Cp_tr_paths[25] += ce * creal(path_storage[25 + idx]);
                    mom_def_Cm_tr_paths[25] += I * ce * cimag(path_storage[25 + idx]);
                    path_storage[26 + idx] = path26(in);
                    mom_def_Cp_tr_paths[26] += ce * creal(path_storage[26 + idx]);
                    mom_def_Cm_tr_paths[26] += I * ce * cimag(path_storage[26 + idx]);
                    path_storage[27 + idx] = path27(in);
                    mom_def_Cp_tr_paths[27] += ce * creal(path_storage[27 + idx]);
                    mom_def_Cm_tr_paths[27] += I * ce * cimag(path_storage[27 + idx]);
                    path_storage[28 + idx] = path28(in);
                    mom_def_Cp_tr_paths[28] += ce * creal(path_storage[28 + idx]);
                    mom_def_Cm_tr_paths[28] += I * ce * cimag(path_storage[28 + idx]);
                    path_storage[29 + idx] = path29(in);
                    mom_def_Cp_tr_paths[29] += ce * creal(path_storage[29 + idx]);
                    mom_def_Cm_tr_paths[29] += I * ce * cimag(path_storage[29 + idx]);
                    path_storage[30 + idx] = path30(in);
                    mom_def_Cp_tr_paths[30] += ce * creal(path_storage[30 + idx]);
                    mom_def_Cm_tr_paths[30] += I * ce * cimag(path_storage[30 + idx]);
                    path_storage[31 + idx] = path31(in);
                    mom_def_Cp_tr_paths[31] += ce * creal(path_storage[31 + idx]);
                    mom_def_Cm_tr_paths[31] += I * ce * cimag(path_storage[31 + idx]);
                    path_storage[32 + idx] = path32(in);
                    mom_def_Cp_tr_paths[32] += ce * creal(path_storage[32 + idx]);
                    mom_def_Cm_tr_paths[32] += I * ce * cimag(path_storage[32 + idx]);
                    path_storage[33 + idx] = path33(in);
                    mom_def_Cp_tr_paths[33] += ce * creal(path_storage[33 + idx]);
                    mom_def_Cm_tr_paths[33] += I * ce * cimag(path_storage[33 + idx]);
                    path_storage[34 + idx] = path34(in);
                    mom_def_Cp_tr_paths[34] += ce * creal(path_storage[34 + idx]);
                    mom_def_Cm_tr_paths[34] += I * ce * cimag(path_storage[34 + idx]);
                    path_storage[35 + idx] = path35(in);
                    mom_def_Cp_tr_paths[35] += ce * creal(path_storage[35 + idx]);
                    mom_def_Cm_tr_paths[35] += I * ce * cimag(path_storage[35 + idx]);
                    path_storage[36 + idx] = path36(in);
                    mom_def_Cp_tr_paths[36] += ce * creal(path_storage[36 + idx]);
                    mom_def_Cm_tr_paths[36] += I * ce * cimag(path_storage[36 + idx]);
                    path_storage[37 + idx] = path37(in);
                    mom_def_Cp_tr_paths[37] += ce * creal(path_storage[37 + idx]);
                    mom_def_Cm_tr_paths[37] += I * ce * cimag(path_storage[37 + idx]);
                    path_storage[38 + idx] = path38(in);
                    mom_def_Cp_tr_paths[38] += ce * creal(path_storage[38 + idx]);
                    mom_def_Cm_tr_paths[38] += I * ce * cimag(path_storage[38 + idx]);
                    path_storage[39 + idx] = path39(in);
                    mom_def_Cp_tr_paths[39] += ce * creal(path_storage[39 + idx]);
                    mom_def_Cm_tr_paths[39] += I * ce * cimag(path_storage[39 + idx]);
                    path_storage[40 + idx] = path40(in);
                    mom_def_Cp_tr_paths[40] += ce * creal(path_storage[40 + idx]);
                    mom_def_Cm_tr_paths[40] += I * ce * cimag(path_storage[40 + idx]);
                    path_storage[41 + idx] = path41(in);
                    mom_def_Cp_tr_paths[41] += ce * creal(path_storage[41 + idx]);
                    mom_def_Cm_tr_paths[41] += I * ce * cimag(path_storage[41 + idx]);
                    path_storage[42 + idx] = path42(in);
                    mom_def_Cp_tr_paths[42] += ce * creal(path_storage[42 + idx]);
                    mom_def_Cm_tr_paths[42] += I * ce * cimag(path_storage[42 + idx]);
                    path_storage[43 + idx] = path43(in);
                    mom_def_Cp_tr_paths[43] += ce * creal(path_storage[43 + idx]);
                    mom_def_Cm_tr_paths[43] += I * ce * cimag(path_storage[43 + idx]);
                    path_storage[44 + idx] = path44(in);
                    mom_def_Cp_tr_paths[44] += ce * creal(path_storage[44 + idx]);
                    mom_def_Cm_tr_paths[44] += I * ce * cimag(path_storage[44 + idx]);
                    path_storage[45 + idx] = path45(in);
                    mom_def_Cp_tr_paths[45] += ce * creal(path_storage[45 + idx]);
                    mom_def_Cm_tr_paths[45] += I * ce * cimag(path_storage[45 + idx]);
                    path_storage[46 + idx] = path46(in);
                    mom_def_Cp_tr_paths[46] += ce * creal(path_storage[46 + idx]);
                    mom_def_Cm_tr_paths[46] += I * ce * cimag(path_storage[46 + idx]);
                    path_storage[47 + idx] = path47(in);
                    mom_def_Cp_tr_paths[47] += ce * creal(path_storage[47 + idx]);
                    mom_def_Cm_tr_paths[47] += I * ce * cimag(path_storage[47 + idx]);
                    path_storage[48 + idx] = path48(in);
                    mom_def_Cp_tr_paths[48] += ce * creal(path_storage[48 + idx]);
                    mom_def_Cm_tr_paths[48] += I * ce * cimag(path_storage[48 + idx]);
                    path_storage[49 + idx] = path49(in);
                    mom_def_Cp_tr_paths[49] += ce * creal(path_storage[49 + idx]);
                    mom_def_Cm_tr_paths[49] += I * ce * cimag(path_storage[49 + idx]);
                    path_storage[50 + idx] = path50(in);
                    mom_def_Cp_tr_paths[50] += ce * creal(path_storage[50 + idx]);
                    mom_def_Cm_tr_paths[50] += I * ce * cimag(path_storage[50 + idx]);
                    path_storage[51 + idx] = path51(in);
                    mom_def_Cp_tr_paths[51] += ce * creal(path_storage[51 + idx]);
                    mom_def_Cm_tr_paths[51] += I * ce * cimag(path_storage[51 + idx]);
                    path_storage[52 + idx] = path52(in);
                    mom_def_Cp_tr_paths[52] += ce * creal(path_storage[52 + idx]);
                    mom_def_Cm_tr_paths[52] += I * ce * cimag(path_storage[52 + idx]);
                    path_storage[53 + idx] = path53(in);
                    mom_def_Cp_tr_paths[53] += ce * creal(path_storage[53 + idx]);
                    mom_def_Cm_tr_paths[53] += I * ce * cimag(path_storage[53 + idx]);
                    path_storage[54 + idx] = path54(in);
                    mom_def_Cp_tr_paths[54] += ce * creal(path_storage[54 + idx]);
                    mom_def_Cm_tr_paths[54] += I * ce * cimag(path_storage[54 + idx]);
                    path_storage[55 + idx] = path55(in);
                    mom_def_Cp_tr_paths[55] += ce * creal(path_storage[55 + idx]);
                    mom_def_Cm_tr_paths[55] += I * ce * cimag(path_storage[55 + idx]);
                    path_storage[56 + idx] = path56(in);
                    mom_def_Cp_tr_paths[56] += ce * creal(path_storage[56 + idx]);
                    mom_def_Cm_tr_paths[56] += I * ce * cimag(path_storage[56 + idx]);
                    path_storage[57 + idx] = path57(in);
                    mom_def_Cp_tr_paths[57] += ce * creal(path_storage[57 + idx]);
                    mom_def_Cm_tr_paths[57] += I * ce * cimag(path_storage[57 + idx]);
                    path_storage[58 + idx] = path58(in);
                    mom_def_Cp_tr_paths[58] += ce * creal(path_storage[58 + idx]);
                    mom_def_Cm_tr_paths[58] += I * ce * cimag(path_storage[58 + idx]);
                    path_storage[59 + idx] = path59(in);
                    mom_def_Cp_tr_paths[59] += ce * creal(path_storage[59 + idx]);
                    mom_def_Cm_tr_paths[59] += I * ce * cimag(path_storage[59 + idx]);
                    path_storage[60 + idx] = path60(in);
                    mom_def_Cp_tr_paths[60] += ce * creal(path_storage[60 + idx]);
                    mom_def_Cm_tr_paths[60] += I * ce * cimag(path_storage[60 + idx]);
                    path_storage[61 + idx] = path61(in);
                    mom_def_Cp_tr_paths[61] += ce * creal(path_storage[61 + idx]);
                    mom_def_Cm_tr_paths[61] += I * ce * cimag(path_storage[61 + idx]);
                    path_storage[62 + idx] = path62(in);
                    mom_def_Cp_tr_paths[62] += ce * creal(path_storage[62 + idx]);
                    mom_def_Cm_tr_paths[62] += I * ce * cimag(path_storage[62 + idx]);
                    path_storage[63 + idx] = path63(in);
                    mom_def_Cp_tr_paths[63] += ce * creal(path_storage[63 + idx]);
                    mom_def_Cm_tr_paths[63] += I * ce * cimag(path_storage[63 + idx]);
                    path_storage[64 + idx] = path64(in);
                    mom_def_Cp_tr_paths[64] += ce * creal(path_storage[64 + idx]);
                    mom_def_Cm_tr_paths[64] += I * ce * cimag(path_storage[64 + idx]);
                    path_storage[65 + idx] = path65(in);
                    mom_def_Cp_tr_paths[65] += ce * creal(path_storage[65 + idx]);
                    mom_def_Cm_tr_paths[65] += I * ce * cimag(path_storage[65 + idx]);
                    path_storage[66 + idx] = path66(in);
                    mom_def_Cp_tr_paths[66] += ce * creal(path_storage[66 + idx]);
                    mom_def_Cm_tr_paths[66] += I * ce * cimag(path_storage[66 + idx]);
                    path_storage[67 + idx] = path67(in);
                    mom_def_Cp_tr_paths[67] += ce * creal(path_storage[67 + idx]);
                    mom_def_Cm_tr_paths[67] += I * ce * cimag(path_storage[67 + idx]);
                    path_storage[68 + idx] = path68(in);
                    mom_def_Cp_tr_paths[68] += ce * creal(path_storage[68 + idx]);
                    mom_def_Cm_tr_paths[68] += I * ce * cimag(path_storage[68 + idx]);
                    path_storage[69 + idx] = path69(in);
                    mom_def_Cp_tr_paths[69] += ce * creal(path_storage[69 + idx]);
                    mom_def_Cm_tr_paths[69] += I * ce * cimag(path_storage[69 + idx]);
                    path_storage[70 + idx] = path70(in);
                    mom_def_Cp_tr_paths[70] += ce * creal(path_storage[70 + idx]);
                    mom_def_Cm_tr_paths[70] += I * ce * cimag(path_storage[70 + idx]);
                    path_storage[71 + idx] = path71(in);
                    mom_def_Cp_tr_paths[71] += ce * creal(path_storage[71 + idx]);
                    mom_def_Cm_tr_paths[71] += I * ce * cimag(path_storage[71 + idx]);
                    path_storage[72 + idx] = path72(in);
                    mom_def_Cp_tr_paths[72] += ce * creal(path_storage[72 + idx]);
                    mom_def_Cm_tr_paths[72] += I * ce * cimag(path_storage[72 + idx]);
                    path_storage[73 + idx] = path73(in);
                    mom_def_Cp_tr_paths[73] += ce * creal(path_storage[73 + idx]);
                    mom_def_Cm_tr_paths[73] += I * ce * cimag(path_storage[73 + idx]);
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
                    for (int i = 0; i < 74; i++)
                    {
                        mom_def_Cp_tr_paths[i] += ce * creal(path_storage[i + idx]);
                        mom_def_Cm_tr_paths[i] += I * ce * cimag(path_storage[i + idx]);
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
    eval_time_momentum_glueball_paths(t, -1, 0, 0);
    OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_2(numerical_op + 1);
    eval_time_momentum_glueball_paths(t, 0, -1, 0);
    OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_1(numerical_op + 2);
    OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_2(numerical_op + 3);
    eval_time_momentum_glueball_paths(t, 0, 0, -1);
    OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_1(numerical_op + 4);
    OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_2(numerical_op + 5);
    eval_time_momentum_glueball_paths(t, 0, 0, 1);
    OP_oneTr_p_0_0_1_Ir_1_C_m1_n_1(numerical_op + 9);
    OP_oneTr_p_0_0_1_Ir_1_C_m1_n_2(numerical_op + 10);
    eval_time_momentum_glueball_paths(t, 0, 1, 0);
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(numerical_op + 11);
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_2(numerical_op + 12);
    eval_time_momentum_glueball_paths(t, 1, 0, 0);
    OP_oneTr_p_1_0_0_Ir_1_C_m1_n_1(numerical_op + 13);
    OP_oneTr_p_1_0_0_Ir_1_C_m1_n_2(numerical_op + 14);
    eval_time_momentum_glueball_paths(t, 1, 1, 0);
    OP_oneTr_p_1_1_0_Ir_4_C_1_n_1(numerical_op + 15);
    *(numerical_op + 6) = +(-0.5) * ((*(numerical_op + 13))) * ((*(numerical_op + 0))) + (-I * 0.5) * ((*(numerical_op + 2))) * ((*(numerical_op + 11))) + (+I * 0.5) * ((*(numerical_op + 12))) * ((*(numerical_op + 3))) + (0.5) * ((*(numerical_op + 1))) * ((*(numerical_op + 14)));
    *(numerical_op + 7) = +(-0.707106781186547524) * ((*(numerical_op + 4))) * ((*(numerical_op + 9))) + (0.707106781186547524) * ((*(numerical_op + 10))) * ((*(numerical_op + 5)));
    *(numerical_op + 8) = +(-0.5) * ((*(numerical_op + 1))) * ((*(numerical_op + 14))) + (-I * 0.5) * ((*(numerical_op + 2))) * ((*(numerical_op + 11))) + (+I * 0.5) * ((*(numerical_op + 12))) * ((*(numerical_op + 3))) + (0.5) * ((*(numerical_op + 13))) * ((*(numerical_op + 0)));
    for (int i = 0; i < total_n_glue_op; i++)
        *(numerical_op_out + i) += *(numerical_op + i);
}

void evaluate_1pt_functions(cor_list *lcor, int nblocking, double complex *gb_storage)
{
    int n1, n2, i;
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
    int t1, t2;

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

    static double complex *gb2;
    static double complex *gb1;
    MPI_Request req_1pt[GLB_T];

    for (int icor = 0; icor < lcor->n_entries; icor++)
    {

        t1 = glbT_to_active_slices[lcor->list[icor].t1];
        t2 = glbT_to_active_slices[lcor->list[icor].t2];
        gb1 = gb1_bf + total_n_glue_op * nblocking * listactive[lcor->list[icor].t1];
        gb2 = gb1_bf + total_n_glue_op * nblocking * listactive[lcor->list[icor].t2];

        if (listsent[lcor->list[icor].t1] == -1)
        {
            if (t1 != -1)
            {
                listsent[lcor->list[icor].t1] = 0;
                if (PID == 0)
                {
                    memcpy(gb1, gb_storage + t1 * total_n_glue_op * nblocking, sizeof(double complex) * total_n_glue_op * nblocking);
                }
                else
                {
                    MPI_Isend((double *)(gb_storage + t1 * total_n_glue_op * nblocking), total_n_glue_op * nblocking * 2, MPI_DOUBLE, 0, lcor->list[icor].t1, cart_comm, req_1pt + lcor->list[icor].t1);
                }
            }

            if (PID == 0 && t1 == -1)
            {
                MPI_Irecv((double *)(gb1), total_n_glue_op * nblocking * 2, MPI_DOUBLE, t_to_proc[lcor->list[icor].t1], lcor->list[icor].t1, cart_comm, req_1pt + lcor->list[icor].t1);
                listsent[lcor->list[icor].t1] = 1;
            }
        }

        if (listsent[lcor->list[icor].t2] == -1)
        {
            listsent[lcor->list[icor].t2] = 0;
            if (t2 != -1)
            {
                if (PID == 0)
                {
                    memcpy(gb2, gb_storage + t2 * total_n_glue_op * nblocking, sizeof(double complex) * total_n_glue_op * nblocking);
                }
                else
                {
                    MPI_Isend((double *)(gb_storage + t2 * total_n_glue_op * nblocking), total_n_glue_op * nblocking * 2, MPI_DOUBLE, 0, lcor->list[icor].t2, cart_comm, req_1pt + lcor->list[icor].t2);
                }
            }

            if (PID == 0 && t2 == -1)
            {
                MPI_Irecv((double *)(gb2), total_n_glue_op * nblocking * 2, MPI_DOUBLE, t_to_proc[lcor->list[icor].t2], lcor->list[icor].t2, cart_comm, req_1pt + lcor->list[icor].t2);
                listsent[lcor->list[icor].t2] = 1;
            }
        }
    }
    if (PID == 0)
        for (i = 0; i < GLB_T; i++)
            if (listsent[i] == 1)
            {
                MPI_Wait(req_1pt + i, MPI_STATUS_IGNORE);
            }

#else
    gb1_bf = gb_storage;
#endif

    lprintf("Measure ML", 0, "\n1pt function P=(-1,0,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n", 2 * nblocking);
    lprintf("Measure ML", 0, "Op id= 0 1 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 0; i < 2; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,-1,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n", 2 * nblocking);
    lprintf("Measure ML", 0, "Op id= 2 3 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 2; i < 4; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,-1) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n", 2 * nblocking);
    lprintf("Measure ML", 0, "Op id= 4 5 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 4; i < 6; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,0) Irrep=T1minusOhP Irrep ev=1/3 Charge=+ nop=%d\n", 1 * nblocking);
    lprintf("Measure ML", 0, "Op id= 6 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 6; i < 7; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,0) Irrep=T1minusOhP Irrep ev=2/3 Charge=+ nop=%d\n", 1 * nblocking);
    lprintf("Measure ML", 0, "Op id= 7 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 7; i < 8; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,0) Irrep=T1minusOhP Irrep ev=3/3 Charge=+ nop=%d\n", 1 * nblocking);
    lprintf("Measure ML", 0, "Op id= 8 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 8; i < 9; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,0,1) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n", 2 * nblocking);
    lprintf("Measure ML", 0, "Op id= 9 10 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 9; i < 11; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(0,1,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n", 2 * nblocking);
    lprintf("Measure ML", 0, "Op id= 11 12 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 11; i < 13; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,0,0) Irrep=A1Dic4 Irrep ev=1/1 Charge=- nop=%d\n", 2 * nblocking);
    lprintf("Measure ML", 0, "Op id= 13 14 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 13; i < 15; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }

    lprintf("Measure ML", 0, "\n1pt function P=(1,1,0) Irrep=B2Dic2 Irrep ev=1/1 Charge=+ nop=%d\n", 1 * nblocking);
    lprintf("Measure ML", 0, "Op id= 15 (repeated nblocking times)\n");
    for (n1 = 0; n1 < GLB_T; n1++)
        if (listactive[n1] > -1)
        {
            lprintf("Measure ML", 0, " t=%d", n1);
            for (n2 = 0; n2 < nblocking; n2++)
                for (i = 15; i < 16; i++)
                    lprintf("Measure ML", 0, " ( %.10e %.10e )", creal(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]),
                            cimag(gb1_bf[i + total_n_glue_op * (n2 + nblocking * listactive[n1])]));
            lprintf("Measure ML", 0, "\n");
        }
}
void report_op_group_setup()
{
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(-1,0,0) Irrep=A1Dic4 Charge=-");
    lprintf("INIT Measure ML", 0, " |0=xyyz-yz-x-y-z-z|1=xyyy-z-y-x-yz-y|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,-1,0) Irrep=A1Dic4 Charge=-");
    lprintf("INIT Measure ML", 0, " |2=xzxy-x-x-x-zx-y|3=xxyzz-x-z-x-y-z|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,0,-1) Irrep=A1Dic4 Charge=-");
    lprintf("INIT Measure ML", 0, " |4=xxx-z-xy-xz-x-y|5=xxzyy-x-y-x-z-y|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,0,0) Irrep=T1minusOhP Charge=+");
    lprintf("INIT Measure ML", 0, " |6=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x,7=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x,8=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,0,1) Irrep=A1Dic4 Charge=-");
    lprintf("INIT Measure ML", 0, " |9=xx-y-x-yz-xyy-z|10=xxxz-xy-x-z-x-y|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(0,1,0) Irrep=A1Dic4 Charge=-");
    lprintf("INIT Measure ML", 0, " |11=xx-z-x-zy-xzz-y|12=xxx-z-xy-xz-x-y|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,0,0) Irrep=A1Dic4 Charge=-");
    lprintf("INIT Measure ML", 0, " |13=yxyzy-x-y-y-y-z|14=xyzz-x-y-y-zy-z|");
    lprintf("INIT Measure ML", 0, "\n1pt Irrep multiplets Total P=(1,1,0) Irrep=B2Dic2 Charge=+");
    lprintf("INIT Measure ML", 0, " |15=x-yz-xy-z|");
}
