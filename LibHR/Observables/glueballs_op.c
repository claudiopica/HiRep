#include <stdlib.h>
#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "utils.h"
#include "glueballs.h"

#define npaths 66
#define PI 3.141592653589793238462643383279502884197
static double complex *Mom_def_Re_tr_paths = NULL;
static double complex *Mom_def_Im_tr_paths = NULL;
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path10(int in)
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    w2 = pu_gauge_wrk(site, 2);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path12(int in)
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    _suNg_trace(p, res);
    return p;
}

static double complex path14(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_m1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[9] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[13] + Mom_def_Im_tr_paths[14];
}

static double complex path0(int in)
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

static double complex path1(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

static double complex path2(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

static void OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[1] + (2.) * Mom_def_Re_tr_paths[2];
}

static double complex path3(int in)
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

static double complex path4(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[4] + (2.) * Mom_def_Re_tr_paths[5];
}

static double complex path6(int in)
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
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

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
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[7] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[9] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[13] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[1] + (2.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[4] + (2.) * Mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[7] + (2.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[9] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[13] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_m1_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[9] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[13] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[13];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[9] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[13] + (-2.) * Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[0] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[3] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[6] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[13];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (-4.) * Mom_def_Re_tr_paths[1] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (-4.) * Mom_def_Re_tr_paths[4] + (2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (-4.) * Mom_def_Re_tr_paths[7] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[9] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[13] + (-2.) * Mom_def_Re_tr_paths[14];
}

static double complex path15(int in)
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_m1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[3];
}

static void OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-4.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-4.) * Mom_def_Im_tr_paths[3];
}

static void OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-4.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_m1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static double complex path19(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

static double complex path20(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

static double complex path21(int in)
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_m1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[21] + Mom_def_Im_tr_paths[22];
}

static double complex c0;
static void OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * c0 * Mom_def_Re_tr_paths[1] + (2.) * c0 * Mom_def_Re_tr_paths[2];
}

static double complex path17(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[17] + (2.) * Mom_def_Re_tr_paths[18];
}

static double complex c1;
static void OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * c1 * Mom_def_Re_tr_paths[7] + (2.) * c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[21] + Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (2.) * c0 * Mom_def_Im_tr_paths[1] + (-2.) * c0 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[17] + (2.) * Mom_def_Im_tr_paths[18];
}

static void OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (2.) * c1 * Mom_def_Im_tr_paths[7] + (-2.) * c1 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[19] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[21] - Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_m1_m1_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[19] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[21] - Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[21];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[21] + (-2.) * Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[0] + (+I * 3.46410161513775459) * c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[3] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[6] + (+I * 3.46410161513775459) * c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[21];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (-4.) * c0 * Mom_def_Re_tr_paths[1] + (2.) * c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (-4.) * Mom_def_Re_tr_paths[17] + (2.) * Mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (-4.) * c1 * Mom_def_Re_tr_paths[7] + (2.) * c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[21] + (-2.) * Mom_def_Re_tr_paths[22];
}

static double complex path23(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] + Mom_def_Im_tr_paths[24];
}

static double complex path25(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
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

    _suNg_trace(p, res);
    return p;
}

static double complex path26(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[25] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[25] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[25] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] - Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[25] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] + Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[25] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[25] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[25] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_m1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] - Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_0_m1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[25] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[17] + (2.) * Mom_def_Im_tr_paths[4];
}

static double complex path27(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path29(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (4.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_m1_0_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[17] + (2.) * Mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_m1_0_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (4.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_m1_0_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[26] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[26] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 5.6568542494923802) * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(+I * 2.8284271247461901) * Mom_def_Im_tr_paths[17] + (-I * 2.8284271247461901) * Mom_def_Im_tr_paths[4];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-I * 5.6568542494923802) * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Im_tr_paths[22] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[28] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[29] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_5(double complex *op_out)
{
    *op_out = +(5.6568542494923802) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_6(double complex *op_out)
{
    *op_out = +(-2.8284271247461901) * Mom_def_Im_tr_paths[23] + (2.8284271247461901) * Mom_def_Im_tr_paths[3];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_7(double complex *op_out)
{
    *op_out = +(5.6568542494923802) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_8(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[27] + (-1.41421356237309505) * Mom_def_Im_tr_paths[26] + (1.41421356237309505) * Mom_def_Im_tr_paths[16] + (1.41421356237309505) * Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(+I * 2.8284271247461901) * Mom_def_Re_tr_paths[17] + (-I * 2.8284271247461901) * Mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Re_tr_paths[22] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[28] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[29] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.8284271247461901) * Mom_def_Re_tr_paths[23] + (2.8284271247461901) * Mom_def_Re_tr_paths[3];
}

static void OP_oneTr_p_m1_0_0_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[27] + (-1.41421356237309505) * Mom_def_Re_tr_paths[26] + (1.41421356237309505) * Mom_def_Re_tr_paths[16] + (1.41421356237309505) * Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[17] + (-2.) * Mom_def_Im_tr_paths[4];
}

static void OP_oneTr_p_m1_0_0_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (-4.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_m1_0_0_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[17] + (-2.) * Mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_m1_0_0_Ir_4_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (-4.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_m1_0_0_Ir_4_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[26] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_0_0_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[26] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[14];
}

static double complex path30(int in)
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
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[30];
}

static double complex path31(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[31] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_0_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[31] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (2.) * c0 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] - Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (2.) * c1 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[31] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] - Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[31] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[31] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_0_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_0_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[31] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (2.) * c0 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] - Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (2.) * c1 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[31] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_0_1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] - Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_m1_0_1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[31] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12];
}

static double complex c2;
static void OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + c2 * Mom_def_Im_tr_paths[1] - c2 * Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[2] - Mom_def_Im_tr_paths[2];
}

static double complex path32(int in)
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
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[30] + Mom_def_Im_tr_paths[33];
}

static double complex c3;
static void OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + c3 * Mom_def_Im_tr_paths[7] - c3 * Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[8] - Mom_def_Im_tr_paths[8];
}

static double complex path34(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static double complex path35(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static double complex path36(int in)
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

static void OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[36] + Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + c2 * Mom_def_Re_tr_paths[1] + c2 * Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[2] + Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[30] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + c3 * Mom_def_Re_tr_paths[7] + c3 * Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[8] + Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[36] + Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + c2 * Mom_def_Im_tr_paths[1] + c2 * Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[2] + Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[17] - Mom_def_Im_tr_paths[32] + Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[30] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + c3 * Mom_def_Im_tr_paths[7] + c3 * Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[8] + Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[36] - Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + c2 * Mom_def_Re_tr_paths[1] - c2 * Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[2] - Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] - Mom_def_Re_tr_paths[17] - Mom_def_Re_tr_paths[32] - Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[30] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + c3 * Mom_def_Re_tr_paths[7] - c3 * Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[8] - Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[36] - Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[2] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[32] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[5] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[8] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[36];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_5(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + (-2.) * c2 * Mom_def_Im_tr_paths[1] + (2.) * c2 * Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[2] - Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[30] + (-2.) * Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_7(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + (-2.) * c3 * Mom_def_Im_tr_paths[7] + (2.) * c3 * Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[8] - Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[34] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[36] + (-2.) * Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[2] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[32] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[5] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[8] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[36];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + (-2.) * c2 * Mom_def_Re_tr_paths[1] + (-2.) * c2 * Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[2] + Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[30] + (-2.) * Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + (-2.) * c3 * Mom_def_Re_tr_paths[7] + (-2.) * c3 * Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[8] + Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[34] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[36] + (-2.) * Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[6];
}

static double complex path37(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_m1_1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (-2.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (-2.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_m1_1_0_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_m1_1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] - Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[1] - c0 * Mom_def_Im_tr_paths[2] + c0 * Mom_def_Im_tr_paths[2];
}

static double complex path38(int in)
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
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[38] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] - Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[7] - c1 * Mom_def_Im_tr_paths[8] + c1 * Mom_def_Im_tr_paths[8];
}

static double complex path39(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

static void OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[1] + c0 * Mom_def_Re_tr_paths[2] + c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[38] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[7] + c1 * Mom_def_Re_tr_paths[8] + c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] - Mom_def_Im_tr_paths[1] - Mom_def_Im_tr_paths[1] - c0 * Mom_def_Im_tr_paths[2] - c0 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[4] - Mom_def_Im_tr_paths[32] + Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[38] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] - Mom_def_Im_tr_paths[7] - Mom_def_Im_tr_paths[7] - c1 * Mom_def_Im_tr_paths[8] - c1 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[39] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[40] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[1] - Mom_def_Re_tr_paths[1] + c0 * Mom_def_Re_tr_paths[2] - c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] - Mom_def_Re_tr_paths[4] - Mom_def_Re_tr_paths[32] - Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[38] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[7] - Mom_def_Re_tr_paths[7] + c1 * Mom_def_Re_tr_paths[8] - c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_2_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[39] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[40] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (-I * 1.73205080756887729) * c0 * Mom_def_Im_tr_paths[2] + (+I * 1.73205080756887729) * c0 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[32] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[18] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (-I * 1.73205080756887729) * c1 * Mom_def_Im_tr_paths[8] + (+I * 1.73205080756887729) * c1 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[40];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_5(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[1] + (-2.) * Mom_def_Im_tr_paths[1] - c0 * Mom_def_Im_tr_paths[2] + c0 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[18] + (-2.) * Mom_def_Im_tr_paths[38] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_7(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[7] + (-2.) * Mom_def_Im_tr_paths[7] - c1 * Mom_def_Im_tr_paths[8] + c1 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[40] + (-2.) * Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (+I * 1.73205080756887729) * c0 * Mom_def_Re_tr_paths[2] + (+I * 1.73205080756887729) * c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[32] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[18] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (+I * 1.73205080756887729) * c1 * Mom_def_Re_tr_paths[8] + (+I * 1.73205080756887729) * c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[40];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + (-2.) * Mom_def_Re_tr_paths[1] + (-2.) * Mom_def_Re_tr_paths[1] + c0 * Mom_def_Re_tr_paths[2] + c0 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[18] + (-2.) * Mom_def_Re_tr_paths[38] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + (-2.) * Mom_def_Re_tr_paths[7] + (-2.) * Mom_def_Re_tr_paths[7] + c1 * Mom_def_Re_tr_paths[8] + c1 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_m1_1_1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[40] + (-2.) * Mom_def_Re_tr_paths[14];
}

static double complex path41(int in)
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

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] - Mom_def_Im_tr_paths[33];
}

static double complex path42(int in)
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

static double complex path43(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static void OP_oneTr_p_0_m1_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] - Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] - Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] - Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] - Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] - Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_m1_m1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] - Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_m1_m1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] - Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[41] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[18] + (-2.) * Mom_def_Im_tr_paths[5];
}

static double complex path44(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

static void OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[44] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (4.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_m1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[41] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[18] + (2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_0_m1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (4.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_m1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[44] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[44] - Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[44] - Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Im_tr_paths[44] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[20] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[11] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[45];
}

static void OP_oneTr_p_0_m1_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[16] + (-1.41421356237309505) * Mom_def_Im_tr_paths[46] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (1.41421356237309505) * Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Re_tr_paths[44] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[20] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[11] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[45];
}

static void OP_oneTr_p_0_m1_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[16] + (-1.41421356237309505) * Mom_def_Re_tr_paths[46] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (1.41421356237309505) * Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[41] + (2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[18] + (2.) * Mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_0_m1_0_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[44] - Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (-4.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_m1_0_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[41] + (2.) * Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[18] + (-2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_0_m1_0_Ir_4_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (-4.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_m1_0_Ir_4_C_1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[44] - Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[44] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_m1_0_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[44] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] - Mom_def_Im_tr_paths[38];
}

static double complex path47(int in)
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    _suNg_trace(p, res);
    return p;
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * c0 * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * c1 * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_m1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (-2.) * c0 * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (-2.) * c1 * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] - Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] - Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * c0 * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_m1_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * c1 * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_m1_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (-2.) * c0 * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (-2.) * c1 * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_m1_1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] - Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_m1_1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[47];
}

static double complex path48(int in)
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

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

static double complex path50(int in)
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static void OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_0_0_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[41] + (2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_0_0_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[41] + (-2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Im_tr_paths[41] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[23] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[3] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Im_tr_paths[48] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[35] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[10] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[49];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[41] + (1.41421356237309505) * Mom_def_Im_tr_paths[23] + (1.41421356237309505) * Mom_def_Im_tr_paths[3] + (1.41421356237309505) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Im_tr_paths[26] + (-1.41421356237309505) * Mom_def_Im_tr_paths[50] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (-1.41421356237309505) * Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Re_tr_paths[41] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[23] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[3] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Re_tr_paths[48] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[35] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[10] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[49];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[41] + (1.41421356237309505) * Mom_def_Re_tr_paths[23] + (1.41421356237309505) * Mom_def_Re_tr_paths[3] + (-1.41421356237309505) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Re_tr_paths[26] + (-1.41421356237309505) * Mom_def_Re_tr_paths[50] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (-1.41421356237309505) * Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[41] + (2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[49] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[49] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[49] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_m1_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[41] + (-2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_m1_Ir_5_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[49] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    w2 = pu_gauge_wrk(site, 1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
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

    site = idn_wrk(site, 2);
    w2 = pu_gauge_wrk(site, 2);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg_dagger(res, res1, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

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

    _suNg_trace(p, res);
    return p;
}

static double complex path60(int in)
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg_dagger(res1, res, *w1);

    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    _suNg_trace(p, res);
    return p;
}

static double complex path61(int in)
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w2 = pu_gauge_wrk(site, 3);

    _suNg_dagger(res1, *w2);
    w2 = &res1;
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, *w2, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

static double complex path63(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
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

    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res1, res, *w1);

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
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
    w1 = pu_gauge_wrk(site, 2);
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

    w2 = pu_gauge_wrk(site, 3);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, *w2, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

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

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg_dagger(res, res1, *w1);

    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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

static double complex path65(int in)
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
    w1 = pu_gauge_wrk(site, 3);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 3);
    site = idn_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
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

    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 1);
    w1 = pu_gauge_wrk(site, 1);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 1);
    site = idn_wrk(site, 3);
    w1 = pu_gauge_wrk(site, 3);
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
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res, res1, *w1);

    site = iup_wrk(site, 2);
    w1 = pu_gauge_wrk(site, 2);
    _suNg_times_suNg(res1, res, *w1);

    site = iup_wrk(site, 2);
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

static void OP_oneTr_p_0_0_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] + Mom_def_Im_tr_paths[51] + Mom_def_Im_tr_paths[52] + Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[53] + Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[54] + Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[9] + Mom_def_Im_tr_paths[55] + Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[56] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] + Mom_def_Im_tr_paths[57] + Mom_def_Im_tr_paths[31] + Mom_def_Im_tr_paths[25] + Mom_def_Im_tr_paths[58] + Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[59] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[44] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[21] + Mom_def_Im_tr_paths[60] + Mom_def_Im_tr_paths[61] + Mom_def_Im_tr_paths[36] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47] + Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[62] + Mom_def_Im_tr_paths[63] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[13] + Mom_def_Im_tr_paths[64] + Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[65] + Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(16.) * Mom_def_Re_tr_paths[0] + (16.) * Mom_def_Re_tr_paths[1] + (16.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[41] + (4.) * Mom_def_Re_tr_paths[23] + (4.) * Mom_def_Re_tr_paths[3] + (4.) * Mom_def_Re_tr_paths[17] + (4.) * Mom_def_Re_tr_paths[4] + (4.) * Mom_def_Re_tr_paths[32] + (4.) * Mom_def_Re_tr_paths[18] + (4.) * Mom_def_Re_tr_paths[5] + (4.) * Mom_def_Re_tr_paths[38] + (4.) * Mom_def_Re_tr_paths[30] + (4.) * Mom_def_Re_tr_paths[33] + (4.) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(16.) * Mom_def_Re_tr_paths[6] + (16.) * Mom_def_Re_tr_paths[7] + (16.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_0_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] + Mom_def_Re_tr_paths[51] + Mom_def_Re_tr_paths[52] + Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[53] + Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[54] + Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[9] + Mom_def_Re_tr_paths[55] + Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[56] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] + Mom_def_Re_tr_paths[57] + Mom_def_Re_tr_paths[31] + Mom_def_Re_tr_paths[25] + Mom_def_Re_tr_paths[58] + Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[59] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[44] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[21] + Mom_def_Re_tr_paths[60] + Mom_def_Re_tr_paths[61] + Mom_def_Re_tr_paths[36] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47] + Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[62] + Mom_def_Re_tr_paths[63] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[13] + Mom_def_Re_tr_paths[64] + Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[65] + Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Im_tr_paths[41] + (4.) * Mom_def_Im_tr_paths[23] + (4.) * Mom_def_Im_tr_paths[3] + (-4.) * Mom_def_Im_tr_paths[17] + (-4.) * Mom_def_Im_tr_paths[4] + (-4.) * Mom_def_Im_tr_paths[32] + (4.) * Mom_def_Im_tr_paths[18] + (4.) * Mom_def_Im_tr_paths[5] + (4.) * Mom_def_Im_tr_paths[38] + (-4.) * Mom_def_Im_tr_paths[30] + (4.) * Mom_def_Im_tr_paths[33] + (-4.) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] + Mom_def_Im_tr_paths[51] + Mom_def_Im_tr_paths[52] + Mom_def_Im_tr_paths[39] - Mom_def_Im_tr_paths[53] - Mom_def_Im_tr_paths[37] - Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[54] + Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[9] + Mom_def_Im_tr_paths[55] - Mom_def_Im_tr_paths[15] - Mom_def_Im_tr_paths[56] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[49] - Mom_def_Im_tr_paths[57] - Mom_def_Im_tr_paths[31] - Mom_def_Im_tr_paths[25] - Mom_def_Im_tr_paths[58] + Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[59] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] - Mom_def_Im_tr_paths[44] - Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[21] + Mom_def_Im_tr_paths[60] + Mom_def_Im_tr_paths[61] + Mom_def_Im_tr_paths[36] - Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[47] - Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[62] + Mom_def_Im_tr_paths[63] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[13] + Mom_def_Im_tr_paths[64] - Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[65] - Mom_def_Im_tr_paths[14] - Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] + Mom_def_Re_tr_paths[51] + Mom_def_Re_tr_paths[52] + Mom_def_Re_tr_paths[39] - Mom_def_Re_tr_paths[53] - Mom_def_Re_tr_paths[37] - Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[54] + Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[9] + Mom_def_Re_tr_paths[55] - Mom_def_Re_tr_paths[15] - Mom_def_Re_tr_paths[56] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[49] - Mom_def_Re_tr_paths[57] - Mom_def_Re_tr_paths[31] - Mom_def_Re_tr_paths[25] - Mom_def_Re_tr_paths[58] + Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[59] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] - Mom_def_Re_tr_paths[44] - Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[21] + Mom_def_Re_tr_paths[60] + Mom_def_Re_tr_paths[61] + Mom_def_Re_tr_paths[36] - Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[47] - Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[62] + Mom_def_Re_tr_paths[63] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[13] + Mom_def_Re_tr_paths[64] - Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[65] - Mom_def_Re_tr_paths[14] - Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(6.) * Mom_def_Im_tr_paths[17] + (6.) * Mom_def_Im_tr_paths[4] + (6.) * Mom_def_Im_tr_paths[18] + (6.) * Mom_def_Im_tr_paths[5] + (-6.) * Mom_def_Im_tr_paths[38] + (-6.) * Mom_def_Im_tr_paths[30] + (-6.) * Mom_def_Im_tr_paths[33] + (-6.) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[34] + (-2.) * Mom_def_Im_tr_paths[51] + (-2.) * Mom_def_Im_tr_paths[52] + (-2.) * Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[53] + Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] + (-2.) * Mom_def_Im_tr_paths[54] + (-2.) * Mom_def_Im_tr_paths[19] + (-2.) * Mom_def_Im_tr_paths[9] + (-2.) * Mom_def_Im_tr_paths[55] + Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[56] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] + (-2.) * Mom_def_Im_tr_paths[57] + (-2.) * Mom_def_Im_tr_paths[31] + (-2.) * Mom_def_Im_tr_paths[25] + (-2.) * Mom_def_Im_tr_paths[58] + Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[59] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + (-2.) * Mom_def_Im_tr_paths[44] + (-2.) * Mom_def_Im_tr_paths[20] + (-2.) * Mom_def_Im_tr_paths[11] + (-2.) * Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[21] + Mom_def_Im_tr_paths[60] + Mom_def_Im_tr_paths[61] + Mom_def_Im_tr_paths[36] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47] + Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[62] + Mom_def_Im_tr_paths[63] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[13] + Mom_def_Im_tr_paths[64] + Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[65] + Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-6.92820323027550917) * Mom_def_Im_tr_paths[41] + (-6.92820323027550917) * Mom_def_Im_tr_paths[23] + (-6.92820323027550917) * Mom_def_Im_tr_paths[3] + (-3.46410161513775459) * Mom_def_Im_tr_paths[17] + (-3.46410161513775459) * Mom_def_Im_tr_paths[4] + (6.92820323027550917) * Mom_def_Im_tr_paths[32] + (3.46410161513775459) * Mom_def_Im_tr_paths[18] + (3.46410161513775459) * Mom_def_Im_tr_paths[5] + (3.46410161513775459) * Mom_def_Im_tr_paths[38] + (-3.46410161513775459) * Mom_def_Im_tr_paths[30] + (3.46410161513775459) * Mom_def_Im_tr_paths[33] + (-3.46410161513775459) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(1.73205080756887729) * Mom_def_Im_tr_paths[53] + (1.73205080756887729) * Mom_def_Im_tr_paths[37] + (1.73205080756887729) * Mom_def_Im_tr_paths[48] + (1.73205080756887729) * Mom_def_Im_tr_paths[35] + (1.73205080756887729) * Mom_def_Im_tr_paths[15] + (1.73205080756887729) * Mom_def_Im_tr_paths[56] + (1.73205080756887729) * Mom_def_Im_tr_paths[10] + (1.73205080756887729) * Mom_def_Im_tr_paths[49] + (-1.73205080756887729) * Mom_def_Im_tr_paths[27] + (-1.73205080756887729) * Mom_def_Im_tr_paths[59] + (-1.73205080756887729) * Mom_def_Im_tr_paths[26] + (-1.73205080756887729) * Mom_def_Im_tr_paths[50] + (-1.73205080756887729) * Mom_def_Im_tr_paths[16] + (-1.73205080756887729) * Mom_def_Im_tr_paths[46] + (-1.73205080756887729) * Mom_def_Im_tr_paths[12] + (-1.73205080756887729) * Mom_def_Im_tr_paths[42] + (1.73205080756887729) * Mom_def_Im_tr_paths[21] + (1.73205080756887729) * Mom_def_Im_tr_paths[60] + (1.73205080756887729) * Mom_def_Im_tr_paths[61] + (1.73205080756887729) * Mom_def_Im_tr_paths[36] + (-1.73205080756887729) * Mom_def_Im_tr_paths[22] + (-1.73205080756887729) * Mom_def_Im_tr_paths[47] + (-1.73205080756887729) * Mom_def_Im_tr_paths[28] + (-1.73205080756887729) * Mom_def_Im_tr_paths[62] + (1.73205080756887729) * Mom_def_Im_tr_paths[63] + (1.73205080756887729) * Mom_def_Im_tr_paths[40] + (1.73205080756887729) * Mom_def_Im_tr_paths[13] + (1.73205080756887729) * Mom_def_Im_tr_paths[64] + (-1.73205080756887729) * Mom_def_Im_tr_paths[29] + (-1.73205080756887729) * Mom_def_Im_tr_paths[65] + (-1.73205080756887729) * Mom_def_Im_tr_paths[14] + (-1.73205080756887729) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(16.) * Mom_def_Re_tr_paths[0] + (-8.) * Mom_def_Re_tr_paths[1] + (-8.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[41] + (4.) * Mom_def_Re_tr_paths[23] + (4.) * Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[17] + (-2.) * Mom_def_Re_tr_paths[4] + (4.) * Mom_def_Re_tr_paths[32] + (-2.) * Mom_def_Re_tr_paths[18] + (-2.) * Mom_def_Re_tr_paths[5] + (-2.) * Mom_def_Re_tr_paths[38] + (-2.) * Mom_def_Re_tr_paths[30] + (-2.) * Mom_def_Re_tr_paths[33] + (-2.) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(16.) * Mom_def_Re_tr_paths[6] + (-8.) * Mom_def_Re_tr_paths[7] + (-8.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[34] + (-2.) * Mom_def_Re_tr_paths[51] + (-2.) * Mom_def_Re_tr_paths[52] + (-2.) * Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[53] + Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] + (-2.) * Mom_def_Re_tr_paths[54] + (-2.) * Mom_def_Re_tr_paths[19] + (-2.) * Mom_def_Re_tr_paths[9] + (-2.) * Mom_def_Re_tr_paths[55] + Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[56] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] + (-2.) * Mom_def_Re_tr_paths[57] + (-2.) * Mom_def_Re_tr_paths[31] + (-2.) * Mom_def_Re_tr_paths[25] + (-2.) * Mom_def_Re_tr_paths[58] + Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[59] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + (-2.) * Mom_def_Re_tr_paths[44] + (-2.) * Mom_def_Re_tr_paths[20] + (-2.) * Mom_def_Re_tr_paths[11] + (-2.) * Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[21] + Mom_def_Re_tr_paths[60] + Mom_def_Re_tr_paths[61] + Mom_def_Re_tr_paths[36] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47] + Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[62] + Mom_def_Re_tr_paths[63] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[13] + Mom_def_Re_tr_paths[64] + Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[65] + Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +(-13.8564064605510183) * Mom_def_Re_tr_paths[1] + (13.8564064605510183) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +(-3.46410161513775459) * Mom_def_Re_tr_paths[17] + (-3.46410161513775459) * Mom_def_Re_tr_paths[4] + (3.46410161513775459) * Mom_def_Re_tr_paths[18] + (3.46410161513775459) * Mom_def_Re_tr_paths[5] + (-3.46410161513775459) * Mom_def_Re_tr_paths[38] + (3.46410161513775459) * Mom_def_Re_tr_paths[30] + (-3.46410161513775459) * Mom_def_Re_tr_paths[33] + (3.46410161513775459) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +(-13.8564064605510183) * Mom_def_Re_tr_paths[7] + (13.8564064605510183) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_0_0_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(1.73205080756887729) * Mom_def_Re_tr_paths[53] + (1.73205080756887729) * Mom_def_Re_tr_paths[37] + (1.73205080756887729) * Mom_def_Re_tr_paths[48] + (1.73205080756887729) * Mom_def_Re_tr_paths[35] + (1.73205080756887729) * Mom_def_Re_tr_paths[15] + (1.73205080756887729) * Mom_def_Re_tr_paths[56] + (1.73205080756887729) * Mom_def_Re_tr_paths[10] + (1.73205080756887729) * Mom_def_Re_tr_paths[49] + (-1.73205080756887729) * Mom_def_Re_tr_paths[27] + (-1.73205080756887729) * Mom_def_Re_tr_paths[59] + (-1.73205080756887729) * Mom_def_Re_tr_paths[26] + (-1.73205080756887729) * Mom_def_Re_tr_paths[50] + (-1.73205080756887729) * Mom_def_Re_tr_paths[16] + (-1.73205080756887729) * Mom_def_Re_tr_paths[46] + (-1.73205080756887729) * Mom_def_Re_tr_paths[12] + (-1.73205080756887729) * Mom_def_Re_tr_paths[42] + (1.73205080756887729) * Mom_def_Re_tr_paths[21] + (1.73205080756887729) * Mom_def_Re_tr_paths[60] + (1.73205080756887729) * Mom_def_Re_tr_paths[61] + (1.73205080756887729) * Mom_def_Re_tr_paths[36] + (-1.73205080756887729) * Mom_def_Re_tr_paths[22] + (-1.73205080756887729) * Mom_def_Re_tr_paths[47] + (-1.73205080756887729) * Mom_def_Re_tr_paths[28] + (-1.73205080756887729) * Mom_def_Re_tr_paths[62] + (1.73205080756887729) * Mom_def_Re_tr_paths[63] + (1.73205080756887729) * Mom_def_Re_tr_paths[40] + (1.73205080756887729) * Mom_def_Re_tr_paths[13] + (1.73205080756887729) * Mom_def_Re_tr_paths[64] + (-1.73205080756887729) * Mom_def_Re_tr_paths[29] + (-1.73205080756887729) * Mom_def_Re_tr_paths[65] + (-1.73205080756887729) * Mom_def_Re_tr_paths[14] + (-1.73205080756887729) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 16.) * Mom_def_Im_tr_paths[1] + (-16.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = +(+I * 4.) * Mom_def_Im_tr_paths[17] + (-I * 4.) * Mom_def_Im_tr_paths[4] + (4.) * Mom_def_Im_tr_paths[18] + (-4.) * Mom_def_Im_tr_paths[5] + (+I * 4.) * Mom_def_Im_tr_paths[38] + (4.) * Mom_def_Im_tr_paths[30] + (-I * 4.) * Mom_def_Im_tr_paths[33] + (-4.) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-I * 16.) * Mom_def_Im_tr_paths[7] + (-16.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = +(1. - I * 1.) * Mom_def_Im_tr_paths[34] + (-1. + I * 1.) * Mom_def_Im_tr_paths[51] + (-1. - I * 1.) * Mom_def_Im_tr_paths[52] + (1. + I * 1.) * Mom_def_Im_tr_paths[39] + (-1. - I * 1.) * Mom_def_Im_tr_paths[53] + (1. - I * 1.) * Mom_def_Im_tr_paths[37] + (1. + I * 1.) * Mom_def_Im_tr_paths[48] + (-1. + I * 1.) * Mom_def_Im_tr_paths[35] + (1. + I * 1.) * Mom_def_Im_tr_paths[54] + (-1. - I * 1.) * Mom_def_Im_tr_paths[19] + (-1. + I * 1.) * Mom_def_Im_tr_paths[9] + (1. - I * 1.) * Mom_def_Im_tr_paths[55] + (-1. + I * 1.) * Mom_def_Im_tr_paths[15] + (1. + I * 1.) * Mom_def_Im_tr_paths[56] + (1. - I * 1.) * Mom_def_Im_tr_paths[10] + (-1. - I * 1.) * Mom_def_Im_tr_paths[49] + (-1. + I * 1.) * Mom_def_Im_tr_paths[57] + (1. - I * 1.) * Mom_def_Im_tr_paths[31] + (-1. - I * 1.) * Mom_def_Im_tr_paths[25] + (1. + I * 1.) * Mom_def_Im_tr_paths[58] + (-1. - I * 1.) * Mom_def_Im_tr_paths[27] + (-1. + I * 1.) * Mom_def_Im_tr_paths[59] + (1. + I * 1.) * Mom_def_Im_tr_paths[26] + (1. - I * 1.) * Mom_def_Im_tr_paths[50] + (1. + I * 1.) * Mom_def_Im_tr_paths[44] + (-1. - I * 1.) * Mom_def_Im_tr_paths[20] + (1. - I * 1.) * Mom_def_Im_tr_paths[11] + (-1. + I * 1.) * Mom_def_Im_tr_paths[45] + (1. - I * 1.) * Mom_def_Im_tr_paths[16] + (1. + I * 1.) * Mom_def_Im_tr_paths[46] + (-1. + I * 1.) * Mom_def_Im_tr_paths[12] + (-1. - I * 1.) * Mom_def_Im_tr_paths[42] + (1. + I * 1.) * Mom_def_Im_tr_paths[21] + (-1. + I * 1.) * Mom_def_Im_tr_paths[60] + (1. - I * 1.) * Mom_def_Im_tr_paths[61] + (-1. - I * 1.) * Mom_def_Im_tr_paths[36] + (1. + I * 1.) * Mom_def_Im_tr_paths[22] + (1. - I * 1.) * Mom_def_Im_tr_paths[47] + (-1. + I * 1.) * Mom_def_Im_tr_paths[28] + (-1. - I * 1.) * Mom_def_Im_tr_paths[62] + (-1. - I * 1.) * Mom_def_Im_tr_paths[63] + (1. - I * 1.) * Mom_def_Im_tr_paths[40] + (-1. + I * 1.) * Mom_def_Im_tr_paths[13] + (1. + I * 1.) * Mom_def_Im_tr_paths[64] + (-1. - I * 1.) * Mom_def_Im_tr_paths[29] + (-1. + I * 1.) * Mom_def_Im_tr_paths[65] + (1. - I * 1.) * Mom_def_Im_tr_paths[14] + (1. + I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_5(double complex *op_out)
{
    *op_out = +(22.6274169979695208) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_6(double complex *op_out)
{
    *op_out = +(-5.6568542494923802) * Mom_def_Im_tr_paths[41] + (-5.6568542494923802) * Mom_def_Im_tr_paths[23] + (5.6568542494923802) * Mom_def_Im_tr_paths[3] + (-5.6568542494923802) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_7(double complex *op_out)
{
    *op_out = +(22.6274169979695208) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_8(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Im_tr_paths[34] + (1.41421356237309505) * Mom_def_Im_tr_paths[51] + (-1.41421356237309505) * Mom_def_Im_tr_paths[52] + (-1.41421356237309505) * Mom_def_Im_tr_paths[39] + (1.41421356237309505) * Mom_def_Im_tr_paths[53] + (-1.41421356237309505) * Mom_def_Im_tr_paths[37] + (1.41421356237309505) * Mom_def_Im_tr_paths[48] + (-1.41421356237309505) * Mom_def_Im_tr_paths[35] + (-1.41421356237309505) * Mom_def_Im_tr_paths[54] + (-1.41421356237309505) * Mom_def_Im_tr_paths[19] + (1.41421356237309505) * Mom_def_Im_tr_paths[9] + (1.41421356237309505) * Mom_def_Im_tr_paths[55] + (-1.41421356237309505) * Mom_def_Im_tr_paths[15] + (1.41421356237309505) * Mom_def_Im_tr_paths[56] + (-1.41421356237309505) * Mom_def_Im_tr_paths[10] + (1.41421356237309505) * Mom_def_Im_tr_paths[49] + (-1.41421356237309505) * Mom_def_Im_tr_paths[57] + (-1.41421356237309505) * Mom_def_Im_tr_paths[31] + (1.41421356237309505) * Mom_def_Im_tr_paths[25] + (1.41421356237309505) * Mom_def_Im_tr_paths[58] + (-1.41421356237309505) * Mom_def_Im_tr_paths[27] + (1.41421356237309505) * Mom_def_Im_tr_paths[59] + (-1.41421356237309505) * Mom_def_Im_tr_paths[26] + (1.41421356237309505) * Mom_def_Im_tr_paths[50] + (1.41421356237309505) * Mom_def_Im_tr_paths[44] + (1.41421356237309505) * Mom_def_Im_tr_paths[20] + (-1.41421356237309505) * Mom_def_Im_tr_paths[11] + (-1.41421356237309505) * Mom_def_Im_tr_paths[45] + (1.41421356237309505) * Mom_def_Im_tr_paths[16] + (-1.41421356237309505) * Mom_def_Im_tr_paths[46] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (-1.41421356237309505) * Mom_def_Im_tr_paths[42] + (-1.41421356237309505) * Mom_def_Im_tr_paths[21] + (1.41421356237309505) * Mom_def_Im_tr_paths[60] + (1.41421356237309505) * Mom_def_Im_tr_paths[61] + (-1.41421356237309505) * Mom_def_Im_tr_paths[36] + (1.41421356237309505) * Mom_def_Im_tr_paths[22] + (-1.41421356237309505) * Mom_def_Im_tr_paths[47] + (-1.41421356237309505) * Mom_def_Im_tr_paths[28] + (1.41421356237309505) * Mom_def_Im_tr_paths[62] + (-1.41421356237309505) * Mom_def_Im_tr_paths[63] + (1.41421356237309505) * Mom_def_Im_tr_paths[40] + (1.41421356237309505) * Mom_def_Im_tr_paths[13] + (-1.41421356237309505) * Mom_def_Im_tr_paths[64] + (1.41421356237309505) * Mom_def_Im_tr_paths[29] + (-1.41421356237309505) * Mom_def_Im_tr_paths[65] + (-1.41421356237309505) * Mom_def_Im_tr_paths[14] + (1.41421356237309505) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_9(double complex *op_out)
{
    *op_out = +(-I * 16.) * Mom_def_Im_tr_paths[1] + (16.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_10(double complex *op_out)
{
    *op_out = +(+I * 4.) * Mom_def_Im_tr_paths[17] + (-I * 4.) * Mom_def_Im_tr_paths[4] + (-4.) * Mom_def_Im_tr_paths[18] + (4.) * Mom_def_Im_tr_paths[5] + (+I * 4.) * Mom_def_Im_tr_paths[38] + (-4.) * Mom_def_Im_tr_paths[30] + (-I * 4.) * Mom_def_Im_tr_paths[33] + (4.) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_11(double complex *op_out)
{
    *op_out = +(-I * 16.) * Mom_def_Im_tr_paths[7] + (16.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_m1_n_12(double complex *op_out)
{
    *op_out = +(-1. - I * 1.) * Mom_def_Im_tr_paths[34] + (1. + I * 1.) * Mom_def_Im_tr_paths[51] + (1. - I * 1.) * Mom_def_Im_tr_paths[52] + (-1. + I * 1.) * Mom_def_Im_tr_paths[39] + (1. - I * 1.) * Mom_def_Im_tr_paths[53] + (-1. - I * 1.) * Mom_def_Im_tr_paths[37] + (-1. + I * 1.) * Mom_def_Im_tr_paths[48] + (1. + I * 1.) * Mom_def_Im_tr_paths[35] + (-1. + I * 1.) * Mom_def_Im_tr_paths[54] + (1. - I * 1.) * Mom_def_Im_tr_paths[19] + (1. + I * 1.) * Mom_def_Im_tr_paths[9] + (-1. - I * 1.) * Mom_def_Im_tr_paths[55] + (1. + I * 1.) * Mom_def_Im_tr_paths[15] + (-1. + I * 1.) * Mom_def_Im_tr_paths[56] + (-1. - I * 1.) * Mom_def_Im_tr_paths[10] + (1. - I * 1.) * Mom_def_Im_tr_paths[49] + (1. + I * 1.) * Mom_def_Im_tr_paths[57] + (-1. - I * 1.) * Mom_def_Im_tr_paths[31] + (1. - I * 1.) * Mom_def_Im_tr_paths[25] + (-1. + I * 1.) * Mom_def_Im_tr_paths[58] + (1. - I * 1.) * Mom_def_Im_tr_paths[27] + (1. + I * 1.) * Mom_def_Im_tr_paths[59] + (-1. + I * 1.) * Mom_def_Im_tr_paths[26] + (-1. - I * 1.) * Mom_def_Im_tr_paths[50] + (-1. + I * 1.) * Mom_def_Im_tr_paths[44] + (1. - I * 1.) * Mom_def_Im_tr_paths[20] + (-1. - I * 1.) * Mom_def_Im_tr_paths[11] + (1. + I * 1.) * Mom_def_Im_tr_paths[45] + (-1. - I * 1.) * Mom_def_Im_tr_paths[16] + (-1. + I * 1.) * Mom_def_Im_tr_paths[46] + (1. + I * 1.) * Mom_def_Im_tr_paths[12] + (1. - I * 1.) * Mom_def_Im_tr_paths[42] + (-1. + I * 1.) * Mom_def_Im_tr_paths[21] + (1. + I * 1.) * Mom_def_Im_tr_paths[60] + (-1. - I * 1.) * Mom_def_Im_tr_paths[61] + (1. - I * 1.) * Mom_def_Im_tr_paths[36] + (-1. + I * 1.) * Mom_def_Im_tr_paths[22] + (-1. - I * 1.) * Mom_def_Im_tr_paths[47] + (1. + I * 1.) * Mom_def_Im_tr_paths[28] + (1. - I * 1.) * Mom_def_Im_tr_paths[62] + (1. - I * 1.) * Mom_def_Im_tr_paths[63] + (-1. - I * 1.) * Mom_def_Im_tr_paths[40] + (1. + I * 1.) * Mom_def_Im_tr_paths[13] + (-1. + I * 1.) * Mom_def_Im_tr_paths[64] + (1. - I * 1.) * Mom_def_Im_tr_paths[29] + (1. + I * 1.) * Mom_def_Im_tr_paths[65] + (-1. - I * 1.) * Mom_def_Im_tr_paths[14] + (-1. + I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = +(1. - I * 1.) * Mom_def_Re_tr_paths[34] + (-1. + I * 1.) * Mom_def_Re_tr_paths[51] + (-1. - I * 1.) * Mom_def_Re_tr_paths[52] + (1. + I * 1.) * Mom_def_Re_tr_paths[39] + (-1. - I * 1.) * Mom_def_Re_tr_paths[53] + (1. - I * 1.) * Mom_def_Re_tr_paths[37] + (1. + I * 1.) * Mom_def_Re_tr_paths[48] + (-1. + I * 1.) * Mom_def_Re_tr_paths[35] + (1. + I * 1.) * Mom_def_Re_tr_paths[54] + (-1. - I * 1.) * Mom_def_Re_tr_paths[19] + (-1. + I * 1.) * Mom_def_Re_tr_paths[9] + (1. - I * 1.) * Mom_def_Re_tr_paths[55] + (-1. + I * 1.) * Mom_def_Re_tr_paths[15] + (1. + I * 1.) * Mom_def_Re_tr_paths[56] + (1. - I * 1.) * Mom_def_Re_tr_paths[10] + (-1. - I * 1.) * Mom_def_Re_tr_paths[49] + (-1. + I * 1.) * Mom_def_Re_tr_paths[57] + (1. - I * 1.) * Mom_def_Re_tr_paths[31] + (-1. - I * 1.) * Mom_def_Re_tr_paths[25] + (1. + I * 1.) * Mom_def_Re_tr_paths[58] + (-1. - I * 1.) * Mom_def_Re_tr_paths[27] + (-1. + I * 1.) * Mom_def_Re_tr_paths[59] + (1. + I * 1.) * Mom_def_Re_tr_paths[26] + (1. - I * 1.) * Mom_def_Re_tr_paths[50] + (1. + I * 1.) * Mom_def_Re_tr_paths[44] + (-1. - I * 1.) * Mom_def_Re_tr_paths[20] + (1. - I * 1.) * Mom_def_Re_tr_paths[11] + (-1. + I * 1.) * Mom_def_Re_tr_paths[45] + (1. - I * 1.) * Mom_def_Re_tr_paths[16] + (1. + I * 1.) * Mom_def_Re_tr_paths[46] + (-1. + I * 1.) * Mom_def_Re_tr_paths[12] + (-1. - I * 1.) * Mom_def_Re_tr_paths[42] + (1. + I * 1.) * Mom_def_Re_tr_paths[21] + (-1. + I * 1.) * Mom_def_Re_tr_paths[60] + (1. - I * 1.) * Mom_def_Re_tr_paths[61] + (-1. - I * 1.) * Mom_def_Re_tr_paths[36] + (1. + I * 1.) * Mom_def_Re_tr_paths[22] + (1. - I * 1.) * Mom_def_Re_tr_paths[47] + (-1. + I * 1.) * Mom_def_Re_tr_paths[28] + (-1. - I * 1.) * Mom_def_Re_tr_paths[62] + (-1. - I * 1.) * Mom_def_Re_tr_paths[63] + (1. - I * 1.) * Mom_def_Re_tr_paths[40] + (-1. + I * 1.) * Mom_def_Re_tr_paths[13] + (1. + I * 1.) * Mom_def_Re_tr_paths[64] + (-1. - I * 1.) * Mom_def_Re_tr_paths[29] + (-1. + I * 1.) * Mom_def_Re_tr_paths[65] + (1. - I * 1.) * Mom_def_Re_tr_paths[14] + (1. + I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Re_tr_paths[34] + (1.41421356237309505) * Mom_def_Re_tr_paths[51] + (-1.41421356237309505) * Mom_def_Re_tr_paths[52] + (-1.41421356237309505) * Mom_def_Re_tr_paths[39] + (1.41421356237309505) * Mom_def_Re_tr_paths[53] + (-1.41421356237309505) * Mom_def_Re_tr_paths[37] + (1.41421356237309505) * Mom_def_Re_tr_paths[48] + (-1.41421356237309505) * Mom_def_Re_tr_paths[35] + (-1.41421356237309505) * Mom_def_Re_tr_paths[54] + (-1.41421356237309505) * Mom_def_Re_tr_paths[19] + (1.41421356237309505) * Mom_def_Re_tr_paths[9] + (1.41421356237309505) * Mom_def_Re_tr_paths[55] + (-1.41421356237309505) * Mom_def_Re_tr_paths[15] + (1.41421356237309505) * Mom_def_Re_tr_paths[56] + (-1.41421356237309505) * Mom_def_Re_tr_paths[10] + (1.41421356237309505) * Mom_def_Re_tr_paths[49] + (-1.41421356237309505) * Mom_def_Re_tr_paths[57] + (-1.41421356237309505) * Mom_def_Re_tr_paths[31] + (1.41421356237309505) * Mom_def_Re_tr_paths[25] + (1.41421356237309505) * Mom_def_Re_tr_paths[58] + (-1.41421356237309505) * Mom_def_Re_tr_paths[27] + (1.41421356237309505) * Mom_def_Re_tr_paths[59] + (-1.41421356237309505) * Mom_def_Re_tr_paths[26] + (1.41421356237309505) * Mom_def_Re_tr_paths[50] + (1.41421356237309505) * Mom_def_Re_tr_paths[44] + (1.41421356237309505) * Mom_def_Re_tr_paths[20] + (-1.41421356237309505) * Mom_def_Re_tr_paths[11] + (-1.41421356237309505) * Mom_def_Re_tr_paths[45] + (1.41421356237309505) * Mom_def_Re_tr_paths[16] + (-1.41421356237309505) * Mom_def_Re_tr_paths[46] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (-1.41421356237309505) * Mom_def_Re_tr_paths[42] + (-1.41421356237309505) * Mom_def_Re_tr_paths[21] + (1.41421356237309505) * Mom_def_Re_tr_paths[60] + (1.41421356237309505) * Mom_def_Re_tr_paths[61] + (-1.41421356237309505) * Mom_def_Re_tr_paths[36] + (1.41421356237309505) * Mom_def_Re_tr_paths[22] + (-1.41421356237309505) * Mom_def_Re_tr_paths[47] + (-1.41421356237309505) * Mom_def_Re_tr_paths[28] + (1.41421356237309505) * Mom_def_Re_tr_paths[62] + (-1.41421356237309505) * Mom_def_Re_tr_paths[63] + (1.41421356237309505) * Mom_def_Re_tr_paths[40] + (1.41421356237309505) * Mom_def_Re_tr_paths[13] + (-1.41421356237309505) * Mom_def_Re_tr_paths[64] + (1.41421356237309505) * Mom_def_Re_tr_paths[29] + (-1.41421356237309505) * Mom_def_Re_tr_paths[65] + (-1.41421356237309505) * Mom_def_Re_tr_paths[14] + (1.41421356237309505) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_4_C_1_n_3(double complex *op_out)
{
    *op_out = +(-1. - I * 1.) * Mom_def_Re_tr_paths[34] + (1. + I * 1.) * Mom_def_Re_tr_paths[51] + (1. - I * 1.) * Mom_def_Re_tr_paths[52] + (-1. + I * 1.) * Mom_def_Re_tr_paths[39] + (1. - I * 1.) * Mom_def_Re_tr_paths[53] + (-1. - I * 1.) * Mom_def_Re_tr_paths[37] + (-1. + I * 1.) * Mom_def_Re_tr_paths[48] + (1. + I * 1.) * Mom_def_Re_tr_paths[35] + (-1. + I * 1.) * Mom_def_Re_tr_paths[54] + (1. - I * 1.) * Mom_def_Re_tr_paths[19] + (1. + I * 1.) * Mom_def_Re_tr_paths[9] + (-1. - I * 1.) * Mom_def_Re_tr_paths[55] + (1. + I * 1.) * Mom_def_Re_tr_paths[15] + (-1. + I * 1.) * Mom_def_Re_tr_paths[56] + (-1. - I * 1.) * Mom_def_Re_tr_paths[10] + (1. - I * 1.) * Mom_def_Re_tr_paths[49] + (1. + I * 1.) * Mom_def_Re_tr_paths[57] + (-1. - I * 1.) * Mom_def_Re_tr_paths[31] + (1. - I * 1.) * Mom_def_Re_tr_paths[25] + (-1. + I * 1.) * Mom_def_Re_tr_paths[58] + (1. - I * 1.) * Mom_def_Re_tr_paths[27] + (1. + I * 1.) * Mom_def_Re_tr_paths[59] + (-1. + I * 1.) * Mom_def_Re_tr_paths[26] + (-1. - I * 1.) * Mom_def_Re_tr_paths[50] + (-1. + I * 1.) * Mom_def_Re_tr_paths[44] + (1. - I * 1.) * Mom_def_Re_tr_paths[20] + (-1. - I * 1.) * Mom_def_Re_tr_paths[11] + (1. + I * 1.) * Mom_def_Re_tr_paths[45] + (-1. - I * 1.) * Mom_def_Re_tr_paths[16] + (-1. + I * 1.) * Mom_def_Re_tr_paths[46] + (1. + I * 1.) * Mom_def_Re_tr_paths[12] + (1. - I * 1.) * Mom_def_Re_tr_paths[42] + (-1. + I * 1.) * Mom_def_Re_tr_paths[21] + (1. + I * 1.) * Mom_def_Re_tr_paths[60] + (-1. - I * 1.) * Mom_def_Re_tr_paths[61] + (1. - I * 1.) * Mom_def_Re_tr_paths[36] + (-1. + I * 1.) * Mom_def_Re_tr_paths[22] + (-1. - I * 1.) * Mom_def_Re_tr_paths[47] + (1. + I * 1.) * Mom_def_Re_tr_paths[28] + (1. - I * 1.) * Mom_def_Re_tr_paths[62] + (1. - I * 1.) * Mom_def_Re_tr_paths[63] + (-1. - I * 1.) * Mom_def_Re_tr_paths[40] + (1. + I * 1.) * Mom_def_Re_tr_paths[13] + (-1. + I * 1.) * Mom_def_Re_tr_paths[64] + (1. - I * 1.) * Mom_def_Re_tr_paths[29] + (1. + I * 1.) * Mom_def_Re_tr_paths[65] + (-1. - I * 1.) * Mom_def_Re_tr_paths[14] + (-1. + I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-1. + I * 1.) * Mom_def_Im_tr_paths[34] + (1. - I * 1.) * Mom_def_Im_tr_paths[51] + (1. + I * 1.) * Mom_def_Im_tr_paths[52] + (-1. - I * 1.) * Mom_def_Im_tr_paths[39] + (-1. + I * 1.) * Mom_def_Im_tr_paths[53] + (1. + I * 1.) * Mom_def_Im_tr_paths[37] + (1. - I * 1.) * Mom_def_Im_tr_paths[48] + (-1. - I * 1.) * Mom_def_Im_tr_paths[35] + (-1. - I * 1.) * Mom_def_Im_tr_paths[54] + (1. + I * 1.) * Mom_def_Im_tr_paths[19] + (1. - I * 1.) * Mom_def_Im_tr_paths[9] + (-1. + I * 1.) * Mom_def_Im_tr_paths[55] + (-1. - I * 1.) * Mom_def_Im_tr_paths[15] + (1. - I * 1.) * Mom_def_Im_tr_paths[56] + (1. + I * 1.) * Mom_def_Im_tr_paths[10] + (-1. + I * 1.) * Mom_def_Im_tr_paths[49] + (1. - I * 1.) * Mom_def_Im_tr_paths[57] + (-1. + I * 1.) * Mom_def_Im_tr_paths[31] + (1. + I * 1.) * Mom_def_Im_tr_paths[25] + (-1. - I * 1.) * Mom_def_Im_tr_paths[58] + (1. - I * 1.) * Mom_def_Im_tr_paths[27] + (1. + I * 1.) * Mom_def_Im_tr_paths[59] + (-1. + I * 1.) * Mom_def_Im_tr_paths[26] + (-1. - I * 1.) * Mom_def_Im_tr_paths[50] + (-1. - I * 1.) * Mom_def_Im_tr_paths[44] + (1. + I * 1.) * Mom_def_Im_tr_paths[20] + (-1. + I * 1.) * Mom_def_Im_tr_paths[11] + (1. - I * 1.) * Mom_def_Im_tr_paths[45] + (-1. - I * 1.) * Mom_def_Im_tr_paths[16] + (-1. + I * 1.) * Mom_def_Im_tr_paths[46] + (1. + I * 1.) * Mom_def_Im_tr_paths[12] + (1. - I * 1.) * Mom_def_Im_tr_paths[42] + (1. + I * 1.) * Mom_def_Im_tr_paths[21] + (-1. + I * 1.) * Mom_def_Im_tr_paths[60] + (1. - I * 1.) * Mom_def_Im_tr_paths[61] + (-1. - I * 1.) * Mom_def_Im_tr_paths[36] + (1. + I * 1.) * Mom_def_Im_tr_paths[22] + (1. - I * 1.) * Mom_def_Im_tr_paths[47] + (-1. + I * 1.) * Mom_def_Im_tr_paths[28] + (-1. - I * 1.) * Mom_def_Im_tr_paths[62] + (-1. - I * 1.) * Mom_def_Im_tr_paths[63] + (1. - I * 1.) * Mom_def_Im_tr_paths[40] + (-1. + I * 1.) * Mom_def_Im_tr_paths[13] + (1. + I * 1.) * Mom_def_Im_tr_paths[64] + (-1. - I * 1.) * Mom_def_Im_tr_paths[29] + (-1. + I * 1.) * Mom_def_Im_tr_paths[65] + (1. - I * 1.) * Mom_def_Im_tr_paths[14] + (1. + I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[34] + (-1.41421356237309505) * Mom_def_Im_tr_paths[51] + (1.41421356237309505) * Mom_def_Im_tr_paths[52] + (1.41421356237309505) * Mom_def_Im_tr_paths[39] + (-1.41421356237309505) * Mom_def_Im_tr_paths[53] + (1.41421356237309505) * Mom_def_Im_tr_paths[37] + (-1.41421356237309505) * Mom_def_Im_tr_paths[48] + (1.41421356237309505) * Mom_def_Im_tr_paths[35] + (1.41421356237309505) * Mom_def_Im_tr_paths[54] + (1.41421356237309505) * Mom_def_Im_tr_paths[19] + (-1.41421356237309505) * Mom_def_Im_tr_paths[9] + (-1.41421356237309505) * Mom_def_Im_tr_paths[55] + (1.41421356237309505) * Mom_def_Im_tr_paths[15] + (-1.41421356237309505) * Mom_def_Im_tr_paths[56] + (1.41421356237309505) * Mom_def_Im_tr_paths[10] + (-1.41421356237309505) * Mom_def_Im_tr_paths[49] + (-1.41421356237309505) * Mom_def_Im_tr_paths[57] + (-1.41421356237309505) * Mom_def_Im_tr_paths[31] + (1.41421356237309505) * Mom_def_Im_tr_paths[25] + (1.41421356237309505) * Mom_def_Im_tr_paths[58] + (-1.41421356237309505) * Mom_def_Im_tr_paths[27] + (1.41421356237309505) * Mom_def_Im_tr_paths[59] + (-1.41421356237309505) * Mom_def_Im_tr_paths[26] + (1.41421356237309505) * Mom_def_Im_tr_paths[50] + (1.41421356237309505) * Mom_def_Im_tr_paths[44] + (1.41421356237309505) * Mom_def_Im_tr_paths[20] + (-1.41421356237309505) * Mom_def_Im_tr_paths[11] + (-1.41421356237309505) * Mom_def_Im_tr_paths[45] + (1.41421356237309505) * Mom_def_Im_tr_paths[16] + (-1.41421356237309505) * Mom_def_Im_tr_paths[46] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (-1.41421356237309505) * Mom_def_Im_tr_paths[42] + (1.41421356237309505) * Mom_def_Im_tr_paths[21] + (-1.41421356237309505) * Mom_def_Im_tr_paths[60] + (-1.41421356237309505) * Mom_def_Im_tr_paths[61] + (1.41421356237309505) * Mom_def_Im_tr_paths[36] + (1.41421356237309505) * Mom_def_Im_tr_paths[22] + (-1.41421356237309505) * Mom_def_Im_tr_paths[47] + (-1.41421356237309505) * Mom_def_Im_tr_paths[28] + (1.41421356237309505) * Mom_def_Im_tr_paths[62] + (1.41421356237309505) * Mom_def_Im_tr_paths[63] + (-1.41421356237309505) * Mom_def_Im_tr_paths[40] + (-1.41421356237309505) * Mom_def_Im_tr_paths[13] + (1.41421356237309505) * Mom_def_Im_tr_paths[64] + (1.41421356237309505) * Mom_def_Im_tr_paths[29] + (-1.41421356237309505) * Mom_def_Im_tr_paths[65] + (-1.41421356237309505) * Mom_def_Im_tr_paths[14] + (1.41421356237309505) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-1. - I * 1.) * Mom_def_Im_tr_paths[34] + (1. + I * 1.) * Mom_def_Im_tr_paths[51] + (1. - I * 1.) * Mom_def_Im_tr_paths[52] + (-1. + I * 1.) * Mom_def_Im_tr_paths[39] + (-1. - I * 1.) * Mom_def_Im_tr_paths[53] + (1. - I * 1.) * Mom_def_Im_tr_paths[37] + (1. + I * 1.) * Mom_def_Im_tr_paths[48] + (-1. + I * 1.) * Mom_def_Im_tr_paths[35] + (-1. + I * 1.) * Mom_def_Im_tr_paths[54] + (1. - I * 1.) * Mom_def_Im_tr_paths[19] + (1. + I * 1.) * Mom_def_Im_tr_paths[9] + (-1. - I * 1.) * Mom_def_Im_tr_paths[55] + (-1. + I * 1.) * Mom_def_Im_tr_paths[15] + (1. + I * 1.) * Mom_def_Im_tr_paths[56] + (1. - I * 1.) * Mom_def_Im_tr_paths[10] + (-1. - I * 1.) * Mom_def_Im_tr_paths[49] + (1. + I * 1.) * Mom_def_Im_tr_paths[57] + (-1. - I * 1.) * Mom_def_Im_tr_paths[31] + (1. - I * 1.) * Mom_def_Im_tr_paths[25] + (-1. + I * 1.) * Mom_def_Im_tr_paths[58] + (1. + I * 1.) * Mom_def_Im_tr_paths[27] + (1. - I * 1.) * Mom_def_Im_tr_paths[59] + (-1. - I * 1.) * Mom_def_Im_tr_paths[26] + (-1. + I * 1.) * Mom_def_Im_tr_paths[50] + (-1. + I * 1.) * Mom_def_Im_tr_paths[44] + (1. - I * 1.) * Mom_def_Im_tr_paths[20] + (-1. - I * 1.) * Mom_def_Im_tr_paths[11] + (1. + I * 1.) * Mom_def_Im_tr_paths[45] + (-1. + I * 1.) * Mom_def_Im_tr_paths[16] + (-1. - I * 1.) * Mom_def_Im_tr_paths[46] + (1. - I * 1.) * Mom_def_Im_tr_paths[12] + (1. + I * 1.) * Mom_def_Im_tr_paths[42] + (1. - I * 1.) * Mom_def_Im_tr_paths[21] + (-1. - I * 1.) * Mom_def_Im_tr_paths[60] + (1. + I * 1.) * Mom_def_Im_tr_paths[61] + (-1. + I * 1.) * Mom_def_Im_tr_paths[36] + (1. - I * 1.) * Mom_def_Im_tr_paths[22] + (1. + I * 1.) * Mom_def_Im_tr_paths[47] + (-1. - I * 1.) * Mom_def_Im_tr_paths[28] + (-1. + I * 1.) * Mom_def_Im_tr_paths[62] + (-1. + I * 1.) * Mom_def_Im_tr_paths[63] + (1. + I * 1.) * Mom_def_Im_tr_paths[40] + (-1. - I * 1.) * Mom_def_Im_tr_paths[13] + (1. - I * 1.) * Mom_def_Im_tr_paths[64] + (-1. + I * 1.) * Mom_def_Im_tr_paths[29] + (-1. - I * 1.) * Mom_def_Im_tr_paths[65] + (1. + I * 1.) * Mom_def_Im_tr_paths[14] + (1. - I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = +(+I * 4.) * Mom_def_Re_tr_paths[17] + (-I * 4.) * Mom_def_Re_tr_paths[4] + (4.) * Mom_def_Re_tr_paths[18] + (-4.) * Mom_def_Re_tr_paths[5] + (-I * 4.) * Mom_def_Re_tr_paths[38] + (-4.) * Mom_def_Re_tr_paths[30] + (+I * 4.) * Mom_def_Re_tr_paths[33] + (4.) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_1_n_2(double complex *op_out)
{
    *op_out = +(-1. + I * 1.) * Mom_def_Re_tr_paths[34] + (1. - I * 1.) * Mom_def_Re_tr_paths[51] + (1. + I * 1.) * Mom_def_Re_tr_paths[52] + (-1. - I * 1.) * Mom_def_Re_tr_paths[39] + (-1. + I * 1.) * Mom_def_Re_tr_paths[53] + (1. + I * 1.) * Mom_def_Re_tr_paths[37] + (1. - I * 1.) * Mom_def_Re_tr_paths[48] + (-1. - I * 1.) * Mom_def_Re_tr_paths[35] + (-1. - I * 1.) * Mom_def_Re_tr_paths[54] + (1. + I * 1.) * Mom_def_Re_tr_paths[19] + (1. - I * 1.) * Mom_def_Re_tr_paths[9] + (-1. + I * 1.) * Mom_def_Re_tr_paths[55] + (-1. - I * 1.) * Mom_def_Re_tr_paths[15] + (1. - I * 1.) * Mom_def_Re_tr_paths[56] + (1. + I * 1.) * Mom_def_Re_tr_paths[10] + (-1. + I * 1.) * Mom_def_Re_tr_paths[49] + (1. - I * 1.) * Mom_def_Re_tr_paths[57] + (-1. + I * 1.) * Mom_def_Re_tr_paths[31] + (1. + I * 1.) * Mom_def_Re_tr_paths[25] + (-1. - I * 1.) * Mom_def_Re_tr_paths[58] + (1. - I * 1.) * Mom_def_Re_tr_paths[27] + (1. + I * 1.) * Mom_def_Re_tr_paths[59] + (-1. + I * 1.) * Mom_def_Re_tr_paths[26] + (-1. - I * 1.) * Mom_def_Re_tr_paths[50] + (-1. - I * 1.) * Mom_def_Re_tr_paths[44] + (1. + I * 1.) * Mom_def_Re_tr_paths[20] + (-1. + I * 1.) * Mom_def_Re_tr_paths[11] + (1. - I * 1.) * Mom_def_Re_tr_paths[45] + (-1. - I * 1.) * Mom_def_Re_tr_paths[16] + (-1. + I * 1.) * Mom_def_Re_tr_paths[46] + (1. + I * 1.) * Mom_def_Re_tr_paths[12] + (1. - I * 1.) * Mom_def_Re_tr_paths[42] + (1. + I * 1.) * Mom_def_Re_tr_paths[21] + (-1. + I * 1.) * Mom_def_Re_tr_paths[60] + (1. - I * 1.) * Mom_def_Re_tr_paths[61] + (-1. - I * 1.) * Mom_def_Re_tr_paths[36] + (1. + I * 1.) * Mom_def_Re_tr_paths[22] + (1. - I * 1.) * Mom_def_Re_tr_paths[47] + (-1. + I * 1.) * Mom_def_Re_tr_paths[28] + (-1. - I * 1.) * Mom_def_Re_tr_paths[62] + (-1. - I * 1.) * Mom_def_Re_tr_paths[63] + (1. - I * 1.) * Mom_def_Re_tr_paths[40] + (-1. + I * 1.) * Mom_def_Re_tr_paths[13] + (1. + I * 1.) * Mom_def_Re_tr_paths[64] + (-1. - I * 1.) * Mom_def_Re_tr_paths[29] + (-1. + I * 1.) * Mom_def_Re_tr_paths[65] + (1. - I * 1.) * Mom_def_Re_tr_paths[14] + (1. + I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_1_n_3(double complex *op_out)
{
    *op_out = +(-5.6568542494923802) * Mom_def_Re_tr_paths[41] + (-5.6568542494923802) * Mom_def_Re_tr_paths[23] + (5.6568542494923802) * Mom_def_Re_tr_paths[3] + (5.6568542494923802) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_1_n_4(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[34] + (-1.41421356237309505) * Mom_def_Re_tr_paths[51] + (1.41421356237309505) * Mom_def_Re_tr_paths[52] + (1.41421356237309505) * Mom_def_Re_tr_paths[39] + (-1.41421356237309505) * Mom_def_Re_tr_paths[53] + (1.41421356237309505) * Mom_def_Re_tr_paths[37] + (-1.41421356237309505) * Mom_def_Re_tr_paths[48] + (1.41421356237309505) * Mom_def_Re_tr_paths[35] + (1.41421356237309505) * Mom_def_Re_tr_paths[54] + (1.41421356237309505) * Mom_def_Re_tr_paths[19] + (-1.41421356237309505) * Mom_def_Re_tr_paths[9] + (-1.41421356237309505) * Mom_def_Re_tr_paths[55] + (1.41421356237309505) * Mom_def_Re_tr_paths[15] + (-1.41421356237309505) * Mom_def_Re_tr_paths[56] + (1.41421356237309505) * Mom_def_Re_tr_paths[10] + (-1.41421356237309505) * Mom_def_Re_tr_paths[49] + (-1.41421356237309505) * Mom_def_Re_tr_paths[57] + (-1.41421356237309505) * Mom_def_Re_tr_paths[31] + (1.41421356237309505) * Mom_def_Re_tr_paths[25] + (1.41421356237309505) * Mom_def_Re_tr_paths[58] + (-1.41421356237309505) * Mom_def_Re_tr_paths[27] + (1.41421356237309505) * Mom_def_Re_tr_paths[59] + (-1.41421356237309505) * Mom_def_Re_tr_paths[26] + (1.41421356237309505) * Mom_def_Re_tr_paths[50] + (1.41421356237309505) * Mom_def_Re_tr_paths[44] + (1.41421356237309505) * Mom_def_Re_tr_paths[20] + (-1.41421356237309505) * Mom_def_Re_tr_paths[11] + (-1.41421356237309505) * Mom_def_Re_tr_paths[45] + (1.41421356237309505) * Mom_def_Re_tr_paths[16] + (-1.41421356237309505) * Mom_def_Re_tr_paths[46] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (-1.41421356237309505) * Mom_def_Re_tr_paths[42] + (1.41421356237309505) * Mom_def_Re_tr_paths[21] + (-1.41421356237309505) * Mom_def_Re_tr_paths[60] + (-1.41421356237309505) * Mom_def_Re_tr_paths[61] + (1.41421356237309505) * Mom_def_Re_tr_paths[36] + (1.41421356237309505) * Mom_def_Re_tr_paths[22] + (-1.41421356237309505) * Mom_def_Re_tr_paths[47] + (-1.41421356237309505) * Mom_def_Re_tr_paths[28] + (1.41421356237309505) * Mom_def_Re_tr_paths[62] + (1.41421356237309505) * Mom_def_Re_tr_paths[63] + (-1.41421356237309505) * Mom_def_Re_tr_paths[40] + (-1.41421356237309505) * Mom_def_Re_tr_paths[13] + (1.41421356237309505) * Mom_def_Re_tr_paths[64] + (1.41421356237309505) * Mom_def_Re_tr_paths[29] + (-1.41421356237309505) * Mom_def_Re_tr_paths[65] + (-1.41421356237309505) * Mom_def_Re_tr_paths[14] + (1.41421356237309505) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_1_n_5(double complex *op_out)
{
    *op_out = +(-I * 4.) * Mom_def_Re_tr_paths[17] + (+I * 4.) * Mom_def_Re_tr_paths[4] + (4.) * Mom_def_Re_tr_paths[18] + (-4.) * Mom_def_Re_tr_paths[5] + (+I * 4.) * Mom_def_Re_tr_paths[38] + (-4.) * Mom_def_Re_tr_paths[30] + (-I * 4.) * Mom_def_Re_tr_paths[33] + (4.) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_5_C_1_n_6(double complex *op_out)
{
    *op_out = +(-1. - I * 1.) * Mom_def_Re_tr_paths[34] + (1. + I * 1.) * Mom_def_Re_tr_paths[51] + (1. - I * 1.) * Mom_def_Re_tr_paths[52] + (-1. + I * 1.) * Mom_def_Re_tr_paths[39] + (-1. - I * 1.) * Mom_def_Re_tr_paths[53] + (1. - I * 1.) * Mom_def_Re_tr_paths[37] + (1. + I * 1.) * Mom_def_Re_tr_paths[48] + (-1. + I * 1.) * Mom_def_Re_tr_paths[35] + (-1. + I * 1.) * Mom_def_Re_tr_paths[54] + (1. - I * 1.) * Mom_def_Re_tr_paths[19] + (1. + I * 1.) * Mom_def_Re_tr_paths[9] + (-1. - I * 1.) * Mom_def_Re_tr_paths[55] + (-1. + I * 1.) * Mom_def_Re_tr_paths[15] + (1. + I * 1.) * Mom_def_Re_tr_paths[56] + (1. - I * 1.) * Mom_def_Re_tr_paths[10] + (-1. - I * 1.) * Mom_def_Re_tr_paths[49] + (1. + I * 1.) * Mom_def_Re_tr_paths[57] + (-1. - I * 1.) * Mom_def_Re_tr_paths[31] + (1. - I * 1.) * Mom_def_Re_tr_paths[25] + (-1. + I * 1.) * Mom_def_Re_tr_paths[58] + (1. + I * 1.) * Mom_def_Re_tr_paths[27] + (1. - I * 1.) * Mom_def_Re_tr_paths[59] + (-1. - I * 1.) * Mom_def_Re_tr_paths[26] + (-1. + I * 1.) * Mom_def_Re_tr_paths[50] + (-1. + I * 1.) * Mom_def_Re_tr_paths[44] + (1. - I * 1.) * Mom_def_Re_tr_paths[20] + (-1. - I * 1.) * Mom_def_Re_tr_paths[11] + (1. + I * 1.) * Mom_def_Re_tr_paths[45] + (-1. + I * 1.) * Mom_def_Re_tr_paths[16] + (-1. - I * 1.) * Mom_def_Re_tr_paths[46] + (1. - I * 1.) * Mom_def_Re_tr_paths[12] + (1. + I * 1.) * Mom_def_Re_tr_paths[42] + (1. - I * 1.) * Mom_def_Re_tr_paths[21] + (-1. - I * 1.) * Mom_def_Re_tr_paths[60] + (1. + I * 1.) * Mom_def_Re_tr_paths[61] + (-1. + I * 1.) * Mom_def_Re_tr_paths[36] + (1. - I * 1.) * Mom_def_Re_tr_paths[22] + (1. + I * 1.) * Mom_def_Re_tr_paths[47] + (-1. - I * 1.) * Mom_def_Re_tr_paths[28] + (-1. + I * 1.) * Mom_def_Re_tr_paths[62] + (-1. + I * 1.) * Mom_def_Re_tr_paths[63] + (1. + I * 1.) * Mom_def_Re_tr_paths[40] + (-1. - I * 1.) * Mom_def_Re_tr_paths[13] + (1. - I * 1.) * Mom_def_Re_tr_paths[64] + (-1. + I * 1.) * Mom_def_Re_tr_paths[29] + (-1. - I * 1.) * Mom_def_Re_tr_paths[65] + (1. + I * 1.) * Mom_def_Re_tr_paths[14] + (1. - I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_6_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] - Mom_def_Im_tr_paths[51] - Mom_def_Im_tr_paths[52] + Mom_def_Im_tr_paths[39] - Mom_def_Im_tr_paths[53] + Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[54] + Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[9] - Mom_def_Im_tr_paths[55] + Mom_def_Im_tr_paths[15] - Mom_def_Im_tr_paths[56] - Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] - Mom_def_Im_tr_paths[57] + Mom_def_Im_tr_paths[31] + Mom_def_Im_tr_paths[25] - Mom_def_Im_tr_paths[58] + Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[59] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[44] - Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[21] - Mom_def_Im_tr_paths[60] - Mom_def_Im_tr_paths[61] + Mom_def_Im_tr_paths[36] - Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47] + Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[62] - Mom_def_Im_tr_paths[63] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[13] - Mom_def_Im_tr_paths[64] + Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[65] - Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_6_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] - Mom_def_Re_tr_paths[51] - Mom_def_Re_tr_paths[52] + Mom_def_Re_tr_paths[39] - Mom_def_Re_tr_paths[53] + Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[54] + Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[9] - Mom_def_Re_tr_paths[55] + Mom_def_Re_tr_paths[15] - Mom_def_Re_tr_paths[56] - Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] - Mom_def_Re_tr_paths[57] + Mom_def_Re_tr_paths[31] + Mom_def_Re_tr_paths[25] - Mom_def_Re_tr_paths[58] + Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[59] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[44] - Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[21] - Mom_def_Re_tr_paths[60] - Mom_def_Re_tr_paths[61] + Mom_def_Re_tr_paths[36] - Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47] + Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[62] - Mom_def_Re_tr_paths[63] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[13] - Mom_def_Re_tr_paths[64] + Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[65] - Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_7_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] - Mom_def_Im_tr_paths[51] - Mom_def_Im_tr_paths[52] + Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[53] - Mom_def_Im_tr_paths[37] - Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[54] + Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[9] - Mom_def_Im_tr_paths[55] - Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[56] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[49] + Mom_def_Im_tr_paths[57] - Mom_def_Im_tr_paths[31] - Mom_def_Im_tr_paths[25] + Mom_def_Im_tr_paths[58] + Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[59] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] - Mom_def_Im_tr_paths[44] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[21] - Mom_def_Im_tr_paths[60] - Mom_def_Im_tr_paths[61] + Mom_def_Im_tr_paths[36] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[47] - Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[62] - Mom_def_Im_tr_paths[63] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[13] - Mom_def_Im_tr_paths[64] - Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[65] + Mom_def_Im_tr_paths[14] - Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_7_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] - Mom_def_Re_tr_paths[51] - Mom_def_Re_tr_paths[52] + Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[53] - Mom_def_Re_tr_paths[37] - Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[54] + Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[9] - Mom_def_Re_tr_paths[55] - Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[56] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[49] + Mom_def_Re_tr_paths[57] - Mom_def_Re_tr_paths[31] - Mom_def_Re_tr_paths[25] + Mom_def_Re_tr_paths[58] + Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[59] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] - Mom_def_Re_tr_paths[44] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[21] - Mom_def_Re_tr_paths[60] - Mom_def_Re_tr_paths[61] + Mom_def_Re_tr_paths[36] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[47] - Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[62] - Mom_def_Re_tr_paths[63] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[13] - Mom_def_Re_tr_paths[64] - Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[65] + Mom_def_Re_tr_paths[14] - Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_8_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[34] + (2.) * Mom_def_Im_tr_paths[51] + (2.) * Mom_def_Im_tr_paths[52] + (-2.) * Mom_def_Im_tr_paths[39] - Mom_def_Im_tr_paths[53] + Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] + (2.) * Mom_def_Im_tr_paths[54] + (-2.) * Mom_def_Im_tr_paths[19] + (-2.) * Mom_def_Im_tr_paths[9] + (2.) * Mom_def_Im_tr_paths[55] + Mom_def_Im_tr_paths[15] - Mom_def_Im_tr_paths[56] - Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] + (2.) * Mom_def_Im_tr_paths[57] + (-2.) * Mom_def_Im_tr_paths[31] + (-2.) * Mom_def_Im_tr_paths[25] + (2.) * Mom_def_Im_tr_paths[58] + Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[59] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + (-2.) * Mom_def_Im_tr_paths[44] + (2.) * Mom_def_Im_tr_paths[20] + (2.) * Mom_def_Im_tr_paths[11] + (-2.) * Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[21] - Mom_def_Im_tr_paths[60] - Mom_def_Im_tr_paths[61] + Mom_def_Im_tr_paths[36] - Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47] + Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[62] - Mom_def_Im_tr_paths[63] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[13] - Mom_def_Im_tr_paths[64] + Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[65] - Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_8_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-1.73205080756887729) * Mom_def_Im_tr_paths[53] + (1.73205080756887729) * Mom_def_Im_tr_paths[37] + (1.73205080756887729) * Mom_def_Im_tr_paths[48] + (-1.73205080756887729) * Mom_def_Im_tr_paths[35] + (1.73205080756887729) * Mom_def_Im_tr_paths[15] + (-1.73205080756887729) * Mom_def_Im_tr_paths[56] + (-1.73205080756887729) * Mom_def_Im_tr_paths[10] + (1.73205080756887729) * Mom_def_Im_tr_paths[49] + (-1.73205080756887729) * Mom_def_Im_tr_paths[27] + (1.73205080756887729) * Mom_def_Im_tr_paths[59] + (1.73205080756887729) * Mom_def_Im_tr_paths[26] + (-1.73205080756887729) * Mom_def_Im_tr_paths[50] + (1.73205080756887729) * Mom_def_Im_tr_paths[16] + (-1.73205080756887729) * Mom_def_Im_tr_paths[46] + (-1.73205080756887729) * Mom_def_Im_tr_paths[12] + (1.73205080756887729) * Mom_def_Im_tr_paths[42] + (1.73205080756887729) * Mom_def_Im_tr_paths[21] + (-1.73205080756887729) * Mom_def_Im_tr_paths[60] + (-1.73205080756887729) * Mom_def_Im_tr_paths[61] + (1.73205080756887729) * Mom_def_Im_tr_paths[36] + (1.73205080756887729) * Mom_def_Im_tr_paths[22] + (-1.73205080756887729) * Mom_def_Im_tr_paths[47] + (-1.73205080756887729) * Mom_def_Im_tr_paths[28] + (1.73205080756887729) * Mom_def_Im_tr_paths[62] + (-1.73205080756887729) * Mom_def_Im_tr_paths[63] + (1.73205080756887729) * Mom_def_Im_tr_paths[40] + (1.73205080756887729) * Mom_def_Im_tr_paths[13] + (-1.73205080756887729) * Mom_def_Im_tr_paths[64] + (-1.73205080756887729) * Mom_def_Im_tr_paths[29] + (1.73205080756887729) * Mom_def_Im_tr_paths[65] + (1.73205080756887729) * Mom_def_Im_tr_paths[14] + (-1.73205080756887729) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_8_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[34] + (2.) * Mom_def_Re_tr_paths[51] + (2.) * Mom_def_Re_tr_paths[52] + (-2.) * Mom_def_Re_tr_paths[39] - Mom_def_Re_tr_paths[53] + Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] + (2.) * Mom_def_Re_tr_paths[54] + (-2.) * Mom_def_Re_tr_paths[19] + (-2.) * Mom_def_Re_tr_paths[9] + (2.) * Mom_def_Re_tr_paths[55] + Mom_def_Re_tr_paths[15] - Mom_def_Re_tr_paths[56] - Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] + (2.) * Mom_def_Re_tr_paths[57] + (-2.) * Mom_def_Re_tr_paths[31] + (-2.) * Mom_def_Re_tr_paths[25] + (2.) * Mom_def_Re_tr_paths[58] + Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[59] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + (-2.) * Mom_def_Re_tr_paths[44] + (2.) * Mom_def_Re_tr_paths[20] + (2.) * Mom_def_Re_tr_paths[11] + (-2.) * Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[21] - Mom_def_Re_tr_paths[60] - Mom_def_Re_tr_paths[61] + Mom_def_Re_tr_paths[36] - Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47] + Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[62] - Mom_def_Re_tr_paths[63] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[13] - Mom_def_Re_tr_paths[64] + Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[65] - Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_8_C_1_n_2(double complex *op_out)
{
    *op_out = +(-1.73205080756887729) * Mom_def_Re_tr_paths[53] + (1.73205080756887729) * Mom_def_Re_tr_paths[37] + (1.73205080756887729) * Mom_def_Re_tr_paths[48] + (-1.73205080756887729) * Mom_def_Re_tr_paths[35] + (1.73205080756887729) * Mom_def_Re_tr_paths[15] + (-1.73205080756887729) * Mom_def_Re_tr_paths[56] + (-1.73205080756887729) * Mom_def_Re_tr_paths[10] + (1.73205080756887729) * Mom_def_Re_tr_paths[49] + (-1.73205080756887729) * Mom_def_Re_tr_paths[27] + (1.73205080756887729) * Mom_def_Re_tr_paths[59] + (1.73205080756887729) * Mom_def_Re_tr_paths[26] + (-1.73205080756887729) * Mom_def_Re_tr_paths[50] + (1.73205080756887729) * Mom_def_Re_tr_paths[16] + (-1.73205080756887729) * Mom_def_Re_tr_paths[46] + (-1.73205080756887729) * Mom_def_Re_tr_paths[12] + (1.73205080756887729) * Mom_def_Re_tr_paths[42] + (1.73205080756887729) * Mom_def_Re_tr_paths[21] + (-1.73205080756887729) * Mom_def_Re_tr_paths[60] + (-1.73205080756887729) * Mom_def_Re_tr_paths[61] + (1.73205080756887729) * Mom_def_Re_tr_paths[36] + (1.73205080756887729) * Mom_def_Re_tr_paths[22] + (-1.73205080756887729) * Mom_def_Re_tr_paths[47] + (-1.73205080756887729) * Mom_def_Re_tr_paths[28] + (1.73205080756887729) * Mom_def_Re_tr_paths[62] + (-1.73205080756887729) * Mom_def_Re_tr_paths[63] + (1.73205080756887729) * Mom_def_Re_tr_paths[40] + (1.73205080756887729) * Mom_def_Re_tr_paths[13] + (-1.73205080756887729) * Mom_def_Re_tr_paths[64] + (-1.73205080756887729) * Mom_def_Re_tr_paths[29] + (1.73205080756887729) * Mom_def_Re_tr_paths[65] + (1.73205080756887729) * Mom_def_Re_tr_paths[14] + (-1.73205080756887729) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_m1_n_1(double complex *op_out)
{
    *op_out = +(1. - I * 1.) * Mom_def_Im_tr_paths[34] + (1. - I * 1.) * Mom_def_Im_tr_paths[51] + (1. + I * 1.) * Mom_def_Im_tr_paths[52] + (1. + I * 1.) * Mom_def_Im_tr_paths[39] + (1. + I * 1.) * Mom_def_Im_tr_paths[53] + (1. - I * 1.) * Mom_def_Im_tr_paths[37] + (1. + I * 1.) * Mom_def_Im_tr_paths[48] + (1. - I * 1.) * Mom_def_Im_tr_paths[35] + (-1. - I * 1.) * Mom_def_Im_tr_paths[54] + (-1. - I * 1.) * Mom_def_Im_tr_paths[19] + (-1. + I * 1.) * Mom_def_Im_tr_paths[9] + (-1. + I * 1.) * Mom_def_Im_tr_paths[55] + (-1. + I * 1.) * Mom_def_Im_tr_paths[15] + (-1. - I * 1.) * Mom_def_Im_tr_paths[56] + (-1. + I * 1.) * Mom_def_Im_tr_paths[10] + (-1. - I * 1.) * Mom_def_Im_tr_paths[49] + (1. - I * 1.) * Mom_def_Im_tr_paths[57] + (1. - I * 1.) * Mom_def_Im_tr_paths[31] + (-1. - I * 1.) * Mom_def_Im_tr_paths[25] + (-1. - I * 1.) * Mom_def_Im_tr_paths[58] + (-1. - I * 1.) * Mom_def_Im_tr_paths[27] + (1. - I * 1.) * Mom_def_Im_tr_paths[59] + (-1. - I * 1.) * Mom_def_Im_tr_paths[26] + (1. - I * 1.) * Mom_def_Im_tr_paths[50] + (1. + I * 1.) * Mom_def_Im_tr_paths[44] + (1. + I * 1.) * Mom_def_Im_tr_paths[20] + (-1. + I * 1.) * Mom_def_Im_tr_paths[11] + (-1. + I * 1.) * Mom_def_Im_tr_paths[45] + (-1. + I * 1.) * Mom_def_Im_tr_paths[16] + (1. + I * 1.) * Mom_def_Im_tr_paths[46] + (-1. + I * 1.) * Mom_def_Im_tr_paths[12] + (1. + I * 1.) * Mom_def_Im_tr_paths[42] + (1. + I * 1.) * Mom_def_Im_tr_paths[21] + (1. - I * 1.) * Mom_def_Im_tr_paths[60] + (-1. + I * 1.) * Mom_def_Im_tr_paths[61] + (-1. - I * 1.) * Mom_def_Im_tr_paths[36] + (-1. - I * 1.) * Mom_def_Im_tr_paths[22] + (1. - I * 1.) * Mom_def_Im_tr_paths[47] + (-1. + I * 1.) * Mom_def_Im_tr_paths[28] + (1. + I * 1.) * Mom_def_Im_tr_paths[62] + (1. + I * 1.) * Mom_def_Im_tr_paths[63] + (1. - I * 1.) * Mom_def_Im_tr_paths[40] + (-1. + I * 1.) * Mom_def_Im_tr_paths[13] + (-1. - I * 1.) * Mom_def_Im_tr_paths[64] + (-1. - I * 1.) * Mom_def_Im_tr_paths[29] + (1. - I * 1.) * Mom_def_Im_tr_paths[65] + (-1. + I * 1.) * Mom_def_Im_tr_paths[14] + (1. + I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_m1_n_2(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Im_tr_paths[34] + (-1.41421356237309505) * Mom_def_Im_tr_paths[51] + (1.41421356237309505) * Mom_def_Im_tr_paths[52] + (-1.41421356237309505) * Mom_def_Im_tr_paths[39] + (-1.41421356237309505) * Mom_def_Im_tr_paths[53] + (-1.41421356237309505) * Mom_def_Im_tr_paths[37] + (1.41421356237309505) * Mom_def_Im_tr_paths[48] + (1.41421356237309505) * Mom_def_Im_tr_paths[35] + (1.41421356237309505) * Mom_def_Im_tr_paths[54] + (-1.41421356237309505) * Mom_def_Im_tr_paths[19] + (1.41421356237309505) * Mom_def_Im_tr_paths[9] + (-1.41421356237309505) * Mom_def_Im_tr_paths[55] + (-1.41421356237309505) * Mom_def_Im_tr_paths[15] + (-1.41421356237309505) * Mom_def_Im_tr_paths[56] + (1.41421356237309505) * Mom_def_Im_tr_paths[10] + (1.41421356237309505) * Mom_def_Im_tr_paths[49] + (1.41421356237309505) * Mom_def_Im_tr_paths[57] + (-1.41421356237309505) * Mom_def_Im_tr_paths[31] + (1.41421356237309505) * Mom_def_Im_tr_paths[25] + (-1.41421356237309505) * Mom_def_Im_tr_paths[58] + (-1.41421356237309505) * Mom_def_Im_tr_paths[27] + (-1.41421356237309505) * Mom_def_Im_tr_paths[59] + (1.41421356237309505) * Mom_def_Im_tr_paths[26] + (1.41421356237309505) * Mom_def_Im_tr_paths[50] + (1.41421356237309505) * Mom_def_Im_tr_paths[44] + (-1.41421356237309505) * Mom_def_Im_tr_paths[20] + (1.41421356237309505) * Mom_def_Im_tr_paths[11] + (-1.41421356237309505) * Mom_def_Im_tr_paths[45] + (-1.41421356237309505) * Mom_def_Im_tr_paths[16] + (-1.41421356237309505) * Mom_def_Im_tr_paths[46] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (1.41421356237309505) * Mom_def_Im_tr_paths[42] + (-1.41421356237309505) * Mom_def_Im_tr_paths[21] + (-1.41421356237309505) * Mom_def_Im_tr_paths[60] + (-1.41421356237309505) * Mom_def_Im_tr_paths[61] + (-1.41421356237309505) * Mom_def_Im_tr_paths[36] + (-1.41421356237309505) * Mom_def_Im_tr_paths[22] + (-1.41421356237309505) * Mom_def_Im_tr_paths[47] + (-1.41421356237309505) * Mom_def_Im_tr_paths[28] + (-1.41421356237309505) * Mom_def_Im_tr_paths[62] + (1.41421356237309505) * Mom_def_Im_tr_paths[63] + (1.41421356237309505) * Mom_def_Im_tr_paths[40] + (1.41421356237309505) * Mom_def_Im_tr_paths[13] + (1.41421356237309505) * Mom_def_Im_tr_paths[64] + (1.41421356237309505) * Mom_def_Im_tr_paths[29] + (1.41421356237309505) * Mom_def_Im_tr_paths[65] + (1.41421356237309505) * Mom_def_Im_tr_paths[14] + (1.41421356237309505) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-1. - I * 1.) * Mom_def_Im_tr_paths[34] + (-1. - I * 1.) * Mom_def_Im_tr_paths[51] + (-1. + I * 1.) * Mom_def_Im_tr_paths[52] + (-1. + I * 1.) * Mom_def_Im_tr_paths[39] + (-1. + I * 1.) * Mom_def_Im_tr_paths[53] + (-1. - I * 1.) * Mom_def_Im_tr_paths[37] + (-1. + I * 1.) * Mom_def_Im_tr_paths[48] + (-1. - I * 1.) * Mom_def_Im_tr_paths[35] + (1. - I * 1.) * Mom_def_Im_tr_paths[54] + (1. - I * 1.) * Mom_def_Im_tr_paths[19] + (1. + I * 1.) * Mom_def_Im_tr_paths[9] + (1. + I * 1.) * Mom_def_Im_tr_paths[55] + (1. + I * 1.) * Mom_def_Im_tr_paths[15] + (1. - I * 1.) * Mom_def_Im_tr_paths[56] + (1. + I * 1.) * Mom_def_Im_tr_paths[10] + (1. - I * 1.) * Mom_def_Im_tr_paths[49] + (-1. - I * 1.) * Mom_def_Im_tr_paths[57] + (-1. - I * 1.) * Mom_def_Im_tr_paths[31] + (1. - I * 1.) * Mom_def_Im_tr_paths[25] + (1. - I * 1.) * Mom_def_Im_tr_paths[58] + (1. - I * 1.) * Mom_def_Im_tr_paths[27] + (-1. - I * 1.) * Mom_def_Im_tr_paths[59] + (1. - I * 1.) * Mom_def_Im_tr_paths[26] + (-1. - I * 1.) * Mom_def_Im_tr_paths[50] + (-1. + I * 1.) * Mom_def_Im_tr_paths[44] + (-1. + I * 1.) * Mom_def_Im_tr_paths[20] + (1. + I * 1.) * Mom_def_Im_tr_paths[11] + (1. + I * 1.) * Mom_def_Im_tr_paths[45] + (1. + I * 1.) * Mom_def_Im_tr_paths[16] + (-1. + I * 1.) * Mom_def_Im_tr_paths[46] + (1. + I * 1.) * Mom_def_Im_tr_paths[12] + (-1. + I * 1.) * Mom_def_Im_tr_paths[42] + (-1. + I * 1.) * Mom_def_Im_tr_paths[21] + (-1. - I * 1.) * Mom_def_Im_tr_paths[60] + (1. + I * 1.) * Mom_def_Im_tr_paths[61] + (1. - I * 1.) * Mom_def_Im_tr_paths[36] + (1. - I * 1.) * Mom_def_Im_tr_paths[22] + (-1. - I * 1.) * Mom_def_Im_tr_paths[47] + (1. + I * 1.) * Mom_def_Im_tr_paths[28] + (-1. + I * 1.) * Mom_def_Im_tr_paths[62] + (-1. + I * 1.) * Mom_def_Im_tr_paths[63] + (-1. - I * 1.) * Mom_def_Im_tr_paths[40] + (1. + I * 1.) * Mom_def_Im_tr_paths[13] + (1. - I * 1.) * Mom_def_Im_tr_paths[64] + (1. - I * 1.) * Mom_def_Im_tr_paths[29] + (-1. - I * 1.) * Mom_def_Im_tr_paths[65] + (1. + I * 1.) * Mom_def_Im_tr_paths[14] + (-1. + I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_1_n_1(double complex *op_out)
{
    *op_out = +(4. + I * 4.) * Mom_def_Re_tr_paths[41] + (-4. - I * 4.) * Mom_def_Re_tr_paths[23] + (-4. + I * 4.) * Mom_def_Re_tr_paths[3] + (-4.) * Mom_def_Re_tr_paths[17] + (-4.) * Mom_def_Re_tr_paths[4] + (4. - I * 4.) * Mom_def_Re_tr_paths[32] + (+I * 4.) * Mom_def_Re_tr_paths[18] + (+I * 4.) * Mom_def_Re_tr_paths[5] + (4.) * Mom_def_Re_tr_paths[38] + (-I * 4.) * Mom_def_Re_tr_paths[30] + (4.) * Mom_def_Re_tr_paths[33] + (-I * 4.) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_1_n_2(double complex *op_out)
{
    *op_out = +(1. - I * 1.) * Mom_def_Re_tr_paths[34] + (1. - I * 1.) * Mom_def_Re_tr_paths[51] + (1. + I * 1.) * Mom_def_Re_tr_paths[52] + (1. + I * 1.) * Mom_def_Re_tr_paths[39] + (1. + I * 1.) * Mom_def_Re_tr_paths[53] + (1. - I * 1.) * Mom_def_Re_tr_paths[37] + (1. + I * 1.) * Mom_def_Re_tr_paths[48] + (1. - I * 1.) * Mom_def_Re_tr_paths[35] + (-1. - I * 1.) * Mom_def_Re_tr_paths[54] + (-1. - I * 1.) * Mom_def_Re_tr_paths[19] + (-1. + I * 1.) * Mom_def_Re_tr_paths[9] + (-1. + I * 1.) * Mom_def_Re_tr_paths[55] + (-1. + I * 1.) * Mom_def_Re_tr_paths[15] + (-1. - I * 1.) * Mom_def_Re_tr_paths[56] + (-1. + I * 1.) * Mom_def_Re_tr_paths[10] + (-1. - I * 1.) * Mom_def_Re_tr_paths[49] + (1. - I * 1.) * Mom_def_Re_tr_paths[57] + (1. - I * 1.) * Mom_def_Re_tr_paths[31] + (-1. - I * 1.) * Mom_def_Re_tr_paths[25] + (-1. - I * 1.) * Mom_def_Re_tr_paths[58] + (-1. - I * 1.) * Mom_def_Re_tr_paths[27] + (1. - I * 1.) * Mom_def_Re_tr_paths[59] + (-1. - I * 1.) * Mom_def_Re_tr_paths[26] + (1. - I * 1.) * Mom_def_Re_tr_paths[50] + (1. + I * 1.) * Mom_def_Re_tr_paths[44] + (1. + I * 1.) * Mom_def_Re_tr_paths[20] + (-1. + I * 1.) * Mom_def_Re_tr_paths[11] + (-1. + I * 1.) * Mom_def_Re_tr_paths[45] + (-1. + I * 1.) * Mom_def_Re_tr_paths[16] + (1. + I * 1.) * Mom_def_Re_tr_paths[46] + (-1. + I * 1.) * Mom_def_Re_tr_paths[12] + (1. + I * 1.) * Mom_def_Re_tr_paths[42] + (1. + I * 1.) * Mom_def_Re_tr_paths[21] + (1. - I * 1.) * Mom_def_Re_tr_paths[60] + (-1. + I * 1.) * Mom_def_Re_tr_paths[61] + (-1. - I * 1.) * Mom_def_Re_tr_paths[36] + (-1. - I * 1.) * Mom_def_Re_tr_paths[22] + (1. - I * 1.) * Mom_def_Re_tr_paths[47] + (-1. + I * 1.) * Mom_def_Re_tr_paths[28] + (1. + I * 1.) * Mom_def_Re_tr_paths[62] + (1. + I * 1.) * Mom_def_Re_tr_paths[63] + (1. - I * 1.) * Mom_def_Re_tr_paths[40] + (-1. + I * 1.) * Mom_def_Re_tr_paths[13] + (-1. - I * 1.) * Mom_def_Re_tr_paths[64] + (-1. - I * 1.) * Mom_def_Re_tr_paths[29] + (1. - I * 1.) * Mom_def_Re_tr_paths[65] + (-1. + I * 1.) * Mom_def_Re_tr_paths[14] + (1. + I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_1_n_3(double complex *op_out)
{
    *op_out = +(-5.6568542494923802) * Mom_def_Re_tr_paths[17] + (5.6568542494923802) * Mom_def_Re_tr_paths[4] + (-5.6568542494923802) * Mom_def_Re_tr_paths[18] + (5.6568542494923802) * Mom_def_Re_tr_paths[5] + (-5.6568542494923802) * Mom_def_Re_tr_paths[38] + (-5.6568542494923802) * Mom_def_Re_tr_paths[30] + (5.6568542494923802) * Mom_def_Re_tr_paths[33] + (5.6568542494923802) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_1_n_4(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Re_tr_paths[34] + (-1.41421356237309505) * Mom_def_Re_tr_paths[51] + (1.41421356237309505) * Mom_def_Re_tr_paths[52] + (-1.41421356237309505) * Mom_def_Re_tr_paths[39] + (-1.41421356237309505) * Mom_def_Re_tr_paths[53] + (-1.41421356237309505) * Mom_def_Re_tr_paths[37] + (1.41421356237309505) * Mom_def_Re_tr_paths[48] + (1.41421356237309505) * Mom_def_Re_tr_paths[35] + (1.41421356237309505) * Mom_def_Re_tr_paths[54] + (-1.41421356237309505) * Mom_def_Re_tr_paths[19] + (1.41421356237309505) * Mom_def_Re_tr_paths[9] + (-1.41421356237309505) * Mom_def_Re_tr_paths[55] + (-1.41421356237309505) * Mom_def_Re_tr_paths[15] + (-1.41421356237309505) * Mom_def_Re_tr_paths[56] + (1.41421356237309505) * Mom_def_Re_tr_paths[10] + (1.41421356237309505) * Mom_def_Re_tr_paths[49] + (1.41421356237309505) * Mom_def_Re_tr_paths[57] + (-1.41421356237309505) * Mom_def_Re_tr_paths[31] + (1.41421356237309505) * Mom_def_Re_tr_paths[25] + (-1.41421356237309505) * Mom_def_Re_tr_paths[58] + (-1.41421356237309505) * Mom_def_Re_tr_paths[27] + (-1.41421356237309505) * Mom_def_Re_tr_paths[59] + (1.41421356237309505) * Mom_def_Re_tr_paths[26] + (1.41421356237309505) * Mom_def_Re_tr_paths[50] + (1.41421356237309505) * Mom_def_Re_tr_paths[44] + (-1.41421356237309505) * Mom_def_Re_tr_paths[20] + (1.41421356237309505) * Mom_def_Re_tr_paths[11] + (-1.41421356237309505) * Mom_def_Re_tr_paths[45] + (-1.41421356237309505) * Mom_def_Re_tr_paths[16] + (-1.41421356237309505) * Mom_def_Re_tr_paths[46] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (1.41421356237309505) * Mom_def_Re_tr_paths[42] + (-1.41421356237309505) * Mom_def_Re_tr_paths[21] + (-1.41421356237309505) * Mom_def_Re_tr_paths[60] + (-1.41421356237309505) * Mom_def_Re_tr_paths[61] + (-1.41421356237309505) * Mom_def_Re_tr_paths[36] + (-1.41421356237309505) * Mom_def_Re_tr_paths[22] + (-1.41421356237309505) * Mom_def_Re_tr_paths[47] + (-1.41421356237309505) * Mom_def_Re_tr_paths[28] + (-1.41421356237309505) * Mom_def_Re_tr_paths[62] + (1.41421356237309505) * Mom_def_Re_tr_paths[63] + (1.41421356237309505) * Mom_def_Re_tr_paths[40] + (1.41421356237309505) * Mom_def_Re_tr_paths[13] + (1.41421356237309505) * Mom_def_Re_tr_paths[64] + (1.41421356237309505) * Mom_def_Re_tr_paths[29] + (1.41421356237309505) * Mom_def_Re_tr_paths[65] + (1.41421356237309505) * Mom_def_Re_tr_paths[14] + (1.41421356237309505) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_1_n_5(double complex *op_out)
{
    *op_out = +(-4. + I * 4.) * Mom_def_Re_tr_paths[41] + (4. - I * 4.) * Mom_def_Re_tr_paths[23] + (4. + I * 4.) * Mom_def_Re_tr_paths[3] + (4.) * Mom_def_Re_tr_paths[17] + (4.) * Mom_def_Re_tr_paths[4] + (-4. - I * 4.) * Mom_def_Re_tr_paths[32] + (+I * 4.) * Mom_def_Re_tr_paths[18] + (+I * 4.) * Mom_def_Re_tr_paths[5] + (-4.) * Mom_def_Re_tr_paths[38] + (-I * 4.) * Mom_def_Re_tr_paths[30] + (-4.) * Mom_def_Re_tr_paths[33] + (-I * 4.) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_9_C_1_n_6(double complex *op_out)
{
    *op_out = +(-1. - I * 1.) * Mom_def_Re_tr_paths[34] + (-1. - I * 1.) * Mom_def_Re_tr_paths[51] + (-1. + I * 1.) * Mom_def_Re_tr_paths[52] + (-1. + I * 1.) * Mom_def_Re_tr_paths[39] + (-1. + I * 1.) * Mom_def_Re_tr_paths[53] + (-1. - I * 1.) * Mom_def_Re_tr_paths[37] + (-1. + I * 1.) * Mom_def_Re_tr_paths[48] + (-1. - I * 1.) * Mom_def_Re_tr_paths[35] + (1. - I * 1.) * Mom_def_Re_tr_paths[54] + (1. - I * 1.) * Mom_def_Re_tr_paths[19] + (1. + I * 1.) * Mom_def_Re_tr_paths[9] + (1. + I * 1.) * Mom_def_Re_tr_paths[55] + (1. + I * 1.) * Mom_def_Re_tr_paths[15] + (1. - I * 1.) * Mom_def_Re_tr_paths[56] + (1. + I * 1.) * Mom_def_Re_tr_paths[10] + (1. - I * 1.) * Mom_def_Re_tr_paths[49] + (-1. - I * 1.) * Mom_def_Re_tr_paths[57] + (-1. - I * 1.) * Mom_def_Re_tr_paths[31] + (1. - I * 1.) * Mom_def_Re_tr_paths[25] + (1. - I * 1.) * Mom_def_Re_tr_paths[58] + (1. - I * 1.) * Mom_def_Re_tr_paths[27] + (-1. - I * 1.) * Mom_def_Re_tr_paths[59] + (1. - I * 1.) * Mom_def_Re_tr_paths[26] + (-1. - I * 1.) * Mom_def_Re_tr_paths[50] + (-1. + I * 1.) * Mom_def_Re_tr_paths[44] + (-1. + I * 1.) * Mom_def_Re_tr_paths[20] + (1. + I * 1.) * Mom_def_Re_tr_paths[11] + (1. + I * 1.) * Mom_def_Re_tr_paths[45] + (1. + I * 1.) * Mom_def_Re_tr_paths[16] + (-1. + I * 1.) * Mom_def_Re_tr_paths[46] + (1. + I * 1.) * Mom_def_Re_tr_paths[12] + (-1. + I * 1.) * Mom_def_Re_tr_paths[42] + (-1. + I * 1.) * Mom_def_Re_tr_paths[21] + (-1. - I * 1.) * Mom_def_Re_tr_paths[60] + (1. + I * 1.) * Mom_def_Re_tr_paths[61] + (1. - I * 1.) * Mom_def_Re_tr_paths[36] + (1. - I * 1.) * Mom_def_Re_tr_paths[22] + (-1. - I * 1.) * Mom_def_Re_tr_paths[47] + (1. + I * 1.) * Mom_def_Re_tr_paths[28] + (-1. + I * 1.) * Mom_def_Re_tr_paths[62] + (-1. + I * 1.) * Mom_def_Re_tr_paths[63] + (-1. - I * 1.) * Mom_def_Re_tr_paths[40] + (1. + I * 1.) * Mom_def_Re_tr_paths[13] + (1. - I * 1.) * Mom_def_Re_tr_paths[64] + (1. - I * 1.) * Mom_def_Re_tr_paths[29] + (-1. - I * 1.) * Mom_def_Re_tr_paths[65] + (1. + I * 1.) * Mom_def_Re_tr_paths[14] + (-1. + I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 5.6568542494923802) * Mom_def_Im_tr_paths[41] + (-I * 5.6568542494923802) * Mom_def_Im_tr_paths[23] + (5.6568542494923802) * Mom_def_Im_tr_paths[3] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[17] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[4] + (5.6568542494923802) * Mom_def_Im_tr_paths[32] + (2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[18] + (2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[5] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[38] + (2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[30] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[33] + (2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_m1_n_2(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Im_tr_paths[34] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[51] + (-1.41421356237309505) * Mom_def_Im_tr_paths[52] + (-1.41421356237309505) * Mom_def_Im_tr_paths[39] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[53] + (1.41421356237309505) * Mom_def_Im_tr_paths[37] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[48] + (1.41421356237309505) * Mom_def_Im_tr_paths[35] + (1.41421356237309505) * Mom_def_Im_tr_paths[54] + (1.41421356237309505) * Mom_def_Im_tr_paths[19] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[9] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[55] + (-1.41421356237309505) * Mom_def_Im_tr_paths[15] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[56] + (-1.41421356237309505) * Mom_def_Im_tr_paths[10] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[49] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[57] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[31] + (1.41421356237309505) * Mom_def_Im_tr_paths[25] + (1.41421356237309505) * Mom_def_Im_tr_paths[58] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[27] + (-1.41421356237309505) * Mom_def_Im_tr_paths[59] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[26] + (-1.41421356237309505) * Mom_def_Im_tr_paths[50] + (-1.41421356237309505) * Mom_def_Im_tr_paths[44] + (-1.41421356237309505) * Mom_def_Im_tr_paths[20] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[11] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[45] + (1.41421356237309505) * Mom_def_Im_tr_paths[16] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[46] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[42] + (1.41421356237309505) * Mom_def_Im_tr_paths[21] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[60] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[61] + (-1.41421356237309505) * Mom_def_Im_tr_paths[36] + (-1.41421356237309505) * Mom_def_Im_tr_paths[22] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[47] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[28] + (1.41421356237309505) * Mom_def_Im_tr_paths[62] + (1.41421356237309505) * Mom_def_Im_tr_paths[63] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[40] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[13] + (-1.41421356237309505) * Mom_def_Im_tr_paths[64] + (-1.41421356237309505) * Mom_def_Im_tr_paths[29] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[65] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[14] + (1.41421356237309505) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-4. + I * 4.) * Mom_def_Im_tr_paths[17] + (4. - I * 4.) * Mom_def_Im_tr_paths[4] + (4. - I * 4.) * Mom_def_Im_tr_paths[18] + (-4. + I * 4.) * Mom_def_Im_tr_paths[5] + (4. - I * 4.) * Mom_def_Im_tr_paths[38] + (-4. + I * 4.) * Mom_def_Im_tr_paths[30] + (-4. + I * 4.) * Mom_def_Im_tr_paths[33] + (4. - I * 4.) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_m1_n_4(double complex *op_out)
{
    *op_out = +(-1. + I * 1.) * Mom_def_Im_tr_paths[34] + (1. - I * 1.) * Mom_def_Im_tr_paths[51] + (-1. + I * 1.) * Mom_def_Im_tr_paths[52] + (1. - I * 1.) * Mom_def_Im_tr_paths[39] + (1. - I * 1.) * Mom_def_Im_tr_paths[53] + (1. - I * 1.) * Mom_def_Im_tr_paths[37] + (-1. + I * 1.) * Mom_def_Im_tr_paths[48] + (-1. + I * 1.) * Mom_def_Im_tr_paths[35] + (-1. + I * 1.) * Mom_def_Im_tr_paths[54] + (1. - I * 1.) * Mom_def_Im_tr_paths[19] + (-1. + I * 1.) * Mom_def_Im_tr_paths[9] + (1. - I * 1.) * Mom_def_Im_tr_paths[55] + (1. - I * 1.) * Mom_def_Im_tr_paths[15] + (1. - I * 1.) * Mom_def_Im_tr_paths[56] + (-1. + I * 1.) * Mom_def_Im_tr_paths[10] + (-1. + I * 1.) * Mom_def_Im_tr_paths[49] + (1. - I * 1.) * Mom_def_Im_tr_paths[57] + (-1. + I * 1.) * Mom_def_Im_tr_paths[31] + (1. - I * 1.) * Mom_def_Im_tr_paths[25] + (-1. + I * 1.) * Mom_def_Im_tr_paths[58] + (-1. + I * 1.) * Mom_def_Im_tr_paths[27] + (-1. + I * 1.) * Mom_def_Im_tr_paths[59] + (1. - I * 1.) * Mom_def_Im_tr_paths[26] + (1. - I * 1.) * Mom_def_Im_tr_paths[50] + (1. - I * 1.) * Mom_def_Im_tr_paths[44] + (-1. + I * 1.) * Mom_def_Im_tr_paths[20] + (1. - I * 1.) * Mom_def_Im_tr_paths[11] + (-1. + I * 1.) * Mom_def_Im_tr_paths[45] + (-1. + I * 1.) * Mom_def_Im_tr_paths[16] + (-1. + I * 1.) * Mom_def_Im_tr_paths[46] + (1. - I * 1.) * Mom_def_Im_tr_paths[12] + (1. - I * 1.) * Mom_def_Im_tr_paths[42] + (1. - I * 1.) * Mom_def_Im_tr_paths[21] + (1. - I * 1.) * Mom_def_Im_tr_paths[60] + (1. - I * 1.) * Mom_def_Im_tr_paths[61] + (1. - I * 1.) * Mom_def_Im_tr_paths[36] + (-1. + I * 1.) * Mom_def_Im_tr_paths[22] + (-1. + I * 1.) * Mom_def_Im_tr_paths[47] + (-1. + I * 1.) * Mom_def_Im_tr_paths[28] + (-1. + I * 1.) * Mom_def_Im_tr_paths[62] + (-1. + I * 1.) * Mom_def_Im_tr_paths[63] + (-1. + I * 1.) * Mom_def_Im_tr_paths[40] + (-1. + I * 1.) * Mom_def_Im_tr_paths[13] + (-1. + I * 1.) * Mom_def_Im_tr_paths[64] + (1. - I * 1.) * Mom_def_Im_tr_paths[29] + (1. - I * 1.) * Mom_def_Im_tr_paths[65] + (1. - I * 1.) * Mom_def_Im_tr_paths[14] + (1. - I * 1.) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_m1_n_5(double complex *op_out)
{
    *op_out = +(-5.6568542494923802) * Mom_def_Im_tr_paths[41] + (5.6568542494923802) * Mom_def_Im_tr_paths[23] + (-I * 5.6568542494923802) * Mom_def_Im_tr_paths[3] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[17] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[4] + (-I * 5.6568542494923802) * Mom_def_Im_tr_paths[32] + (-2.8284271247461901 - I * 2.8284271247461901) * Mom_def_Im_tr_paths[18] + (-2.8284271247461901 - I * 2.8284271247461901) * Mom_def_Im_tr_paths[5] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[38] + (-2.8284271247461901 - I * 2.8284271247461901) * Mom_def_Im_tr_paths[30] + (-2.8284271247461901 + I * 2.8284271247461901) * Mom_def_Im_tr_paths[33] + (-2.8284271247461901 - I * 2.8284271247461901) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_m1_n_6(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[34] + (-1.41421356237309505) * Mom_def_Im_tr_paths[51] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[52] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[39] + (1.41421356237309505) * Mom_def_Im_tr_paths[53] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[37] + (1.41421356237309505) * Mom_def_Im_tr_paths[48] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[35] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[54] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[19] + (1.41421356237309505) * Mom_def_Im_tr_paths[9] + (1.41421356237309505) * Mom_def_Im_tr_paths[55] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[15] + (-1.41421356237309505) * Mom_def_Im_tr_paths[56] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[10] + (-1.41421356237309505) * Mom_def_Im_tr_paths[49] + (-1.41421356237309505) * Mom_def_Im_tr_paths[57] + (-1.41421356237309505) * Mom_def_Im_tr_paths[31] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[25] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[58] + (1.41421356237309505) * Mom_def_Im_tr_paths[27] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[59] + (1.41421356237309505) * Mom_def_Im_tr_paths[26] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[50] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[44] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[20] + (1.41421356237309505) * Mom_def_Im_tr_paths[11] + (1.41421356237309505) * Mom_def_Im_tr_paths[45] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[16] + (-1.41421356237309505) * Mom_def_Im_tr_paths[46] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[12] + (-1.41421356237309505) * Mom_def_Im_tr_paths[42] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[21] + (1.41421356237309505) * Mom_def_Im_tr_paths[60] + (-1.41421356237309505) * Mom_def_Im_tr_paths[61] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[36] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[22] + (1.41421356237309505) * Mom_def_Im_tr_paths[47] + (-1.41421356237309505) * Mom_def_Im_tr_paths[28] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[62] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[63] + (1.41421356237309505) * Mom_def_Im_tr_paths[40] + (-1.41421356237309505) * Mom_def_Im_tr_paths[13] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[64] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[29] + (1.41421356237309505) * Mom_def_Im_tr_paths[65] + (-1.41421356237309505) * Mom_def_Im_tr_paths[14] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Re_tr_paths[34] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[51] + (-1.41421356237309505) * Mom_def_Re_tr_paths[52] + (-1.41421356237309505) * Mom_def_Re_tr_paths[39] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[53] + (1.41421356237309505) * Mom_def_Re_tr_paths[37] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[48] + (1.41421356237309505) * Mom_def_Re_tr_paths[35] + (1.41421356237309505) * Mom_def_Re_tr_paths[54] + (1.41421356237309505) * Mom_def_Re_tr_paths[19] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[9] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[55] + (-1.41421356237309505) * Mom_def_Re_tr_paths[15] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[56] + (-1.41421356237309505) * Mom_def_Re_tr_paths[10] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[49] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[57] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[31] + (1.41421356237309505) * Mom_def_Re_tr_paths[25] + (1.41421356237309505) * Mom_def_Re_tr_paths[58] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[27] + (-1.41421356237309505) * Mom_def_Re_tr_paths[59] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[26] + (-1.41421356237309505) * Mom_def_Re_tr_paths[50] + (-1.41421356237309505) * Mom_def_Re_tr_paths[44] + (-1.41421356237309505) * Mom_def_Re_tr_paths[20] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[11] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[45] + (1.41421356237309505) * Mom_def_Re_tr_paths[16] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[46] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[42] + (1.41421356237309505) * Mom_def_Re_tr_paths[21] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[60] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[61] + (-1.41421356237309505) * Mom_def_Re_tr_paths[36] + (-1.41421356237309505) * Mom_def_Re_tr_paths[22] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[47] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[28] + (1.41421356237309505) * Mom_def_Re_tr_paths[62] + (1.41421356237309505) * Mom_def_Re_tr_paths[63] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[40] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[13] + (-1.41421356237309505) * Mom_def_Re_tr_paths[64] + (-1.41421356237309505) * Mom_def_Re_tr_paths[29] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[65] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[14] + (1.41421356237309505) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_1_n_2(double complex *op_out)
{
    *op_out = +(-1. + I * 1.) * Mom_def_Re_tr_paths[34] + (1. - I * 1.) * Mom_def_Re_tr_paths[51] + (-1. + I * 1.) * Mom_def_Re_tr_paths[52] + (1. - I * 1.) * Mom_def_Re_tr_paths[39] + (1. - I * 1.) * Mom_def_Re_tr_paths[53] + (1. - I * 1.) * Mom_def_Re_tr_paths[37] + (-1. + I * 1.) * Mom_def_Re_tr_paths[48] + (-1. + I * 1.) * Mom_def_Re_tr_paths[35] + (-1. + I * 1.) * Mom_def_Re_tr_paths[54] + (1. - I * 1.) * Mom_def_Re_tr_paths[19] + (-1. + I * 1.) * Mom_def_Re_tr_paths[9] + (1. - I * 1.) * Mom_def_Re_tr_paths[55] + (1. - I * 1.) * Mom_def_Re_tr_paths[15] + (1. - I * 1.) * Mom_def_Re_tr_paths[56] + (-1. + I * 1.) * Mom_def_Re_tr_paths[10] + (-1. + I * 1.) * Mom_def_Re_tr_paths[49] + (1. - I * 1.) * Mom_def_Re_tr_paths[57] + (-1. + I * 1.) * Mom_def_Re_tr_paths[31] + (1. - I * 1.) * Mom_def_Re_tr_paths[25] + (-1. + I * 1.) * Mom_def_Re_tr_paths[58] + (-1. + I * 1.) * Mom_def_Re_tr_paths[27] + (-1. + I * 1.) * Mom_def_Re_tr_paths[59] + (1. - I * 1.) * Mom_def_Re_tr_paths[26] + (1. - I * 1.) * Mom_def_Re_tr_paths[50] + (1. - I * 1.) * Mom_def_Re_tr_paths[44] + (-1. + I * 1.) * Mom_def_Re_tr_paths[20] + (1. - I * 1.) * Mom_def_Re_tr_paths[11] + (-1. + I * 1.) * Mom_def_Re_tr_paths[45] + (-1. + I * 1.) * Mom_def_Re_tr_paths[16] + (-1. + I * 1.) * Mom_def_Re_tr_paths[46] + (1. - I * 1.) * Mom_def_Re_tr_paths[12] + (1. - I * 1.) * Mom_def_Re_tr_paths[42] + (1. - I * 1.) * Mom_def_Re_tr_paths[21] + (1. - I * 1.) * Mom_def_Re_tr_paths[60] + (1. - I * 1.) * Mom_def_Re_tr_paths[61] + (1. - I * 1.) * Mom_def_Re_tr_paths[36] + (-1. + I * 1.) * Mom_def_Re_tr_paths[22] + (-1. + I * 1.) * Mom_def_Re_tr_paths[47] + (-1. + I * 1.) * Mom_def_Re_tr_paths[28] + (-1. + I * 1.) * Mom_def_Re_tr_paths[62] + (-1. + I * 1.) * Mom_def_Re_tr_paths[63] + (-1. + I * 1.) * Mom_def_Re_tr_paths[40] + (-1. + I * 1.) * Mom_def_Re_tr_paths[13] + (-1. + I * 1.) * Mom_def_Re_tr_paths[64] + (1. - I * 1.) * Mom_def_Re_tr_paths[29] + (1. - I * 1.) * Mom_def_Re_tr_paths[65] + (1. - I * 1.) * Mom_def_Re_tr_paths[14] + (1. - I * 1.) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_0_Ir_10_C_1_n_3(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[34] + (-1.41421356237309505) * Mom_def_Re_tr_paths[51] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[52] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[39] + (1.41421356237309505) * Mom_def_Re_tr_paths[53] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[37] + (1.41421356237309505) * Mom_def_Re_tr_paths[48] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[35] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[54] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[19] + (1.41421356237309505) * Mom_def_Re_tr_paths[9] + (1.41421356237309505) * Mom_def_Re_tr_paths[55] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[15] + (-1.41421356237309505) * Mom_def_Re_tr_paths[56] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[10] + (-1.41421356237309505) * Mom_def_Re_tr_paths[49] + (-1.41421356237309505) * Mom_def_Re_tr_paths[57] + (-1.41421356237309505) * Mom_def_Re_tr_paths[31] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[25] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[58] + (1.41421356237309505) * Mom_def_Re_tr_paths[27] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[59] + (1.41421356237309505) * Mom_def_Re_tr_paths[26] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[50] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[44] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[20] + (1.41421356237309505) * Mom_def_Re_tr_paths[11] + (1.41421356237309505) * Mom_def_Re_tr_paths[45] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[16] + (-1.41421356237309505) * Mom_def_Re_tr_paths[46] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[12] + (-1.41421356237309505) * Mom_def_Re_tr_paths[42] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[21] + (1.41421356237309505) * Mom_def_Re_tr_paths[60] + (-1.41421356237309505) * Mom_def_Re_tr_paths[61] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[36] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[22] + (1.41421356237309505) * Mom_def_Re_tr_paths[47] + (-1.41421356237309505) * Mom_def_Re_tr_paths[28] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[62] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[63] + (1.41421356237309505) * Mom_def_Re_tr_paths[40] + (-1.41421356237309505) * Mom_def_Re_tr_paths[13] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[64] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[29] + (1.41421356237309505) * Mom_def_Re_tr_paths[65] + (-1.41421356237309505) * Mom_def_Re_tr_paths[14] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[41] + (2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_0_0_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_0_0_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[41] + (-2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(8.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_0_0_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[49] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[49] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Im_tr_paths[41] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[23] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[3] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Im_tr_paths[48] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[35] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[10] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[49];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[41] + (1.41421356237309505) * Mom_def_Im_tr_paths[23] + (1.41421356237309505) * Mom_def_Im_tr_paths[3] + (1.41421356237309505) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Im_tr_paths[26] + (-1.41421356237309505) * Mom_def_Im_tr_paths[50] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (-1.41421356237309505) * Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Re_tr_paths[41] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[23] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[3] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(+I * 1.41421356237309505) * Mom_def_Re_tr_paths[48] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[35] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[10] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[49];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[41] + (1.41421356237309505) * Mom_def_Re_tr_paths[23] + (1.41421356237309505) * Mom_def_Re_tr_paths[3] + (-1.41421356237309505) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(1.41421356237309505) * Mom_def_Re_tr_paths[26] + (-1.41421356237309505) * Mom_def_Re_tr_paths[50] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (-1.41421356237309505) * Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[41] + (2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[48] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[49] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[48] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[49] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[48] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[49] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[50] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_0_1_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[41] + (-2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_0_0_1_Ir_5_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[48] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[49] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[50] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] - Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * c2 * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * c3 * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (-2.) * c2 * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (-2.) * c3 * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] - Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] - Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * c2 * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_1_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * c3 * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_1_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (-2.) * c2 * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (-2.) * c3 * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[47];
}

static void OP_oneTr_p_0_1_m1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] - Mom_def_Re_tr_paths[38];
}

static void OP_oneTr_p_0_1_m1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[47];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[41] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[18] + (-2.) * Mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[44] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (4.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[41] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[18] + (2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (4.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[44] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[44] - Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[44] - Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Im_tr_paths[44] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[20] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[11] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[45];
}

static void OP_oneTr_p_0_1_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[16] + (-1.41421356237309505) * Mom_def_Im_tr_paths[46] + (1.41421356237309505) * Mom_def_Im_tr_paths[12] + (1.41421356237309505) * Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Re_tr_paths[44] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[20] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[11] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[45];
}

static void OP_oneTr_p_0_1_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[16] + (-1.41421356237309505) * Mom_def_Re_tr_paths[46] + (1.41421356237309505) * Mom_def_Re_tr_paths[12] + (1.41421356237309505) * Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[41] + (2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[18] + (2.) * Mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_0_1_0_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[44] - Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[45] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (-4.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_0_1_0_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[41] + (2.) * Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[18] + (-2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_0_1_0_Ir_4_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (-4.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_0_1_0_Ir_4_C_1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[44] - Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[45] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[44] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[45] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[46] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42];
}

static void OP_oneTr_p_0_1_0_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[44] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[45] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[46] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] - Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_1_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_1_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] - Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] - Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] - Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] - Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] - Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] + Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_0_1_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_0_1_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] + Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_0_1_1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[41] + Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_0_1_1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[42] + Mom_def_Im_tr_paths[14] - Mom_def_Im_tr_paths[43];
}

static void OP_oneTr_p_0_1_1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[41] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] - Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_0_1_1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[42] + Mom_def_Re_tr_paths[14] - Mom_def_Re_tr_paths[43];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] - Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[1] - c2 * Mom_def_Im_tr_paths[2] + c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[38] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] - Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[7] - c3 * Mom_def_Im_tr_paths[8] + c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[40] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[1] + c2 * Mom_def_Re_tr_paths[2] + c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[38] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[7] + c3 * Mom_def_Re_tr_paths[8] + c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[40] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] - Mom_def_Im_tr_paths[1] - Mom_def_Im_tr_paths[1] - c2 * Mom_def_Im_tr_paths[2] - c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[4] - Mom_def_Im_tr_paths[32] + Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[38] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] - Mom_def_Im_tr_paths[7] - Mom_def_Im_tr_paths[7] - c3 * Mom_def_Im_tr_paths[8] - c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[39] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[40] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[1] - Mom_def_Re_tr_paths[1] + c2 * Mom_def_Re_tr_paths[2] - c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] - Mom_def_Re_tr_paths[4] - Mom_def_Re_tr_paths[32] - Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[38] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[7] - Mom_def_Re_tr_paths[7] + c3 * Mom_def_Re_tr_paths[8] - c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[39] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[40] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (-I * 1.73205080756887729) * c2 * Mom_def_Im_tr_paths[2] + (+I * 1.73205080756887729) * c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[32] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[18] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (-I * 1.73205080756887729) * c3 * Mom_def_Im_tr_paths[8] + (+I * 1.73205080756887729) * c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[40];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_5(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[1] + (-2.) * Mom_def_Im_tr_paths[1] - c2 * Mom_def_Im_tr_paths[2] + c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[4] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[18] + (-2.) * Mom_def_Im_tr_paths[38] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_7(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[7] + (-2.) * Mom_def_Im_tr_paths[7] - c3 * Mom_def_Im_tr_paths[8] + c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[39] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[40] + (-2.) * Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (+I * 1.73205080756887729) * c2 * Mom_def_Re_tr_paths[2] + (+I * 1.73205080756887729) * c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[32] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[18] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (+I * 1.73205080756887729) * c3 * Mom_def_Re_tr_paths[8] + (+I * 1.73205080756887729) * c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[40];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + (-2.) * Mom_def_Re_tr_paths[1] + (-2.) * Mom_def_Re_tr_paths[1] + c2 * Mom_def_Re_tr_paths[2] + c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[4] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[18] + (-2.) * Mom_def_Re_tr_paths[38] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + (-2.) * Mom_def_Re_tr_paths[7] + (-2.) * Mom_def_Re_tr_paths[7] + c3 * Mom_def_Re_tr_paths[8] + c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[39] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[40] + (-2.) * Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (2.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (2.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_1_m1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[32];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (-2.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[32];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (-2.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_1_m1_0_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[37] + Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_m1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[37] + Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + c0 * Mom_def_Im_tr_paths[1] - c0 * Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[2] - Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[30] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + c1 * Mom_def_Im_tr_paths[7] - c1 * Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[8] - Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[36] + Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + c0 * Mom_def_Re_tr_paths[1] + c0 * Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[2] + Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[30] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + c1 * Mom_def_Re_tr_paths[7] + c1 * Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[8] + Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[36] + Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + c0 * Mom_def_Im_tr_paths[1] + c0 * Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[2] + Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[17] - Mom_def_Im_tr_paths[32] + Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[30] + Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + c1 * Mom_def_Im_tr_paths[7] + c1 * Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[8] + Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[34] - Mom_def_Im_tr_paths[35] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[36] - Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + c0 * Mom_def_Re_tr_paths[1] - c0 * Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[2] - Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] - Mom_def_Re_tr_paths[17] - Mom_def_Re_tr_paths[32] - Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[30] + Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_1_n_3(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + c1 * Mom_def_Re_tr_paths[7] - c1 * Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[8] - Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_2_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[34] - Mom_def_Re_tr_paths[35] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[36] - Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[0] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[2] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[32] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[5] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(+I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[6] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[8] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[36];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_5(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[0] + Mom_def_Im_tr_paths[0] + (-2.) * c0 * Mom_def_Im_tr_paths[1] + (2.) * c0 * Mom_def_Im_tr_paths[1] + Mom_def_Im_tr_paths[2] - Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[17] + Mom_def_Im_tr_paths[32] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[30] + (-2.) * Mom_def_Im_tr_paths[33];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_7(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[6] + Mom_def_Im_tr_paths[6] + (-2.) * c1 * Mom_def_Im_tr_paths[7] + (2.) * c1 * Mom_def_Im_tr_paths[7] + Mom_def_Im_tr_paths[8] - Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[34] + Mom_def_Im_tr_paths[35] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[36] + (-2.) * Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[0] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[2] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[3] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[32] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[5] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[6] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[8] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[35] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[36];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[0] + Mom_def_Re_tr_paths[0] + (-2.) * c0 * Mom_def_Re_tr_paths[1] + (-2.) * c0 * Mom_def_Re_tr_paths[1] + Mom_def_Re_tr_paths[2] + Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[17] + Mom_def_Re_tr_paths[32] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[30] + (-2.) * Mom_def_Re_tr_paths[33];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[6] + Mom_def_Re_tr_paths[6] + (-2.) * c1 * Mom_def_Re_tr_paths[7] + (-2.) * c1 * Mom_def_Re_tr_paths[7] + Mom_def_Re_tr_paths[8] + Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_m1_1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[34] + Mom_def_Re_tr_paths[35] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[36] + (-2.) * Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[31] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_0_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[31] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (2.) * c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] - Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (2.) * c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[31] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] - Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[31] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] + Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[31] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_0_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] + Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_0_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[31] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (2.) * c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[18] - Mom_def_Im_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (2.) * c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[31] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_m1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[18] - Mom_def_Re_tr_paths[30];
}

static void OP_oneTr_p_1_0_m1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[31] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (2.) * Mom_def_Im_tr_paths[17] + (2.) * Mom_def_Im_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (4.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[17] + (2.) * Mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (4.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_1_0_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[26] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[22] + Mom_def_Im_tr_paths[28] + Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[26] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[22] + Mom_def_Re_tr_paths[28] + Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(+I * 5.6568542494923802) * Mom_def_Im_tr_paths[1];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-I * 2.8284271247461901) * Mom_def_Im_tr_paths[17] + (+I * 2.8284271247461901) * Mom_def_Im_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(+I * 5.6568542494923802) * Mom_def_Im_tr_paths[7];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Im_tr_paths[22] + (-I * 1.41421356237309505) * Mom_def_Im_tr_paths[28] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[29] + (+I * 1.41421356237309505) * Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_5(double complex *op_out)
{
    *op_out = +(5.6568542494923802) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_6(double complex *op_out)
{
    *op_out = +(-2.8284271247461901) * Mom_def_Im_tr_paths[23] + (2.8284271247461901) * Mom_def_Im_tr_paths[3];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_7(double complex *op_out)
{
    *op_out = +(5.6568542494923802) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_m1_n_8(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Im_tr_paths[27] + (-1.41421356237309505) * Mom_def_Im_tr_paths[26] + (1.41421356237309505) * Mom_def_Im_tr_paths[16] + (1.41421356237309505) * Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 2.8284271247461901) * Mom_def_Re_tr_paths[17] + (+I * 2.8284271247461901) * Mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 1.41421356237309505) * Mom_def_Re_tr_paths[22] + (-I * 1.41421356237309505) * Mom_def_Re_tr_paths[28] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[29] + (+I * 1.41421356237309505) * Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.8284271247461901) * Mom_def_Re_tr_paths[23] + (2.8284271247461901) * Mom_def_Re_tr_paths[3];
}

static void OP_oneTr_p_1_0_0_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-1.41421356237309505) * Mom_def_Re_tr_paths[27] + (-1.41421356237309505) * Mom_def_Re_tr_paths[26] + (1.41421356237309505) * Mom_def_Re_tr_paths[16] + (1.41421356237309505) * Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[23] + (2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[17] + (-2.) * Mom_def_Im_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] - Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[29] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0] + (-4.) * Mom_def_Re_tr_paths[1];
}

static void OP_oneTr_p_1_0_0_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[23] + (2.) * Mom_def_Re_tr_paths[3] + (-2.) * Mom_def_Re_tr_paths[17] + (-2.) * Mom_def_Re_tr_paths[4];
}

static void OP_oneTr_p_1_0_0_Ir_4_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6] + (-4.) * Mom_def_Re_tr_paths[7];
}

static void OP_oneTr_p_1_0_0_Ir_4_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] - Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[29] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_5_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[27] - Mom_def_Im_tr_paths[26] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[22] - Mom_def_Im_tr_paths[28] - Mom_def_Im_tr_paths[29] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_0_0_Ir_5_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[27] - Mom_def_Re_tr_paths[26] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[22] - Mom_def_Re_tr_paths[28] - Mom_def_Re_tr_paths[29] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] + Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[25] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_0_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[25] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_0_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_0_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[25] + Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] - Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_2_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[25] + Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] - Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] + Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[25] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_0_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] - Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] + Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_0_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[25] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_0_1_Ir_4_C_m1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[23] + Mom_def_Im_tr_paths[3] - Mom_def_Im_tr_paths[5] - Mom_def_Im_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_4_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_0_1_Ir_4_C_m1_n_4(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[25] - Mom_def_Im_tr_paths[26] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_0_1_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[23] + Mom_def_Re_tr_paths[3] + Mom_def_Re_tr_paths[5] - Mom_def_Re_tr_paths[24];
}

static void OP_oneTr_p_1_0_1_Ir_4_C_1_n_2(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[25] - Mom_def_Re_tr_paths[26] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_1_m1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[21] + Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_1_m1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * c2 * Mom_def_Re_tr_paths[1] + (2.) * c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_1_m1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[17] + (2.) * Mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_1_1_m1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * c3 * Mom_def_Re_tr_paths[7] + (2.) * c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_1_m1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[21] + Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (2.) * c2 * Mom_def_Im_tr_paths[1] + (-2.) * c2 * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[17] + (2.) * Mom_def_Im_tr_paths[18];
}

static void OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (2.) * c3 * Mom_def_Im_tr_paths[7] + (-2.) * c3 * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[19] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[21] - Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_1_m1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[19] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[21] - Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[21];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[19] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[20] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[21] + (-2.) * Mom_def_Im_tr_paths[22];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[0] + (+I * 3.46410161513775459) * c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[3] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[6] + (+I * 3.46410161513775459) * c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[20] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[21];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (-4.) * c2 * Mom_def_Re_tr_paths[1] + (2.) * c2 * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (-4.) * Mom_def_Re_tr_paths[17] + (2.) * Mom_def_Re_tr_paths[18];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (-4.) * c3 * Mom_def_Re_tr_paths[7] + (2.) * c3 * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_1_m1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[19] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[20] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[21] + (-2.) * Mom_def_Re_tr_paths[22];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[0];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[3];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(4.) * Mom_def_Re_tr_paths[6];
}

static void OP_oneTr_p_1_1_0_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-4.) * Mom_def_Im_tr_paths[0];
}

static void OP_oneTr_p_1_1_0_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-4.) * Mom_def_Im_tr_paths[3];
}

static void OP_oneTr_p_1_1_0_Ir_3_C_m1_n_3(double complex *op_out)
{
    *op_out = +(-4.) * Mom_def_Im_tr_paths[6];
}

static void OP_oneTr_p_1_1_0_Ir_3_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[16] - Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[16] - Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_4_C_m1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Im_tr_paths[15] + Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[16] + Mom_def_Im_tr_paths[12];
}

static void OP_oneTr_p_1_1_0_Ir_4_C_1_n_1(double complex *op_out)
{
    *op_out = -Mom_def_Re_tr_paths[15] + Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[16] + Mom_def_Re_tr_paths[12];
}

static void OP_oneTr_p_1_1_1_Ir_1_C_m1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[9] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[13] + Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_1_1_Ir_1_C_1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (2.) * Mom_def_Re_tr_paths[1] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_1_1_Ir_1_C_1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (2.) * Mom_def_Re_tr_paths[4] + (2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_1_1_1_Ir_1_C_1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (2.) * Mom_def_Re_tr_paths[7] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_1_1_Ir_1_C_1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[9] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[13] + Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_1_1_Ir_2_C_m1_n_1(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[0] + (-2.) * Mom_def_Im_tr_paths[1] + (2.) * Mom_def_Im_tr_paths[2];
}

static void OP_oneTr_p_1_1_1_Ir_2_C_m1_n_2(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[3] + (-2.) * Mom_def_Im_tr_paths[4] + (2.) * Mom_def_Im_tr_paths[5];
}

static void OP_oneTr_p_1_1_1_Ir_2_C_m1_n_3(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Im_tr_paths[6] + (-2.) * Mom_def_Im_tr_paths[7] + (2.) * Mom_def_Im_tr_paths[8];
}

static void OP_oneTr_p_1_1_1_Ir_2_C_m1_n_4(double complex *op_out)
{
    *op_out = +Mom_def_Im_tr_paths[9] - Mom_def_Im_tr_paths[10] - Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[13] - Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_1_1_Ir_2_C_1_n_1(double complex *op_out)
{
    *op_out = +Mom_def_Re_tr_paths[9] - Mom_def_Re_tr_paths[10] - Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[13] - Mom_def_Re_tr_paths[14];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_m1_n_1(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Im_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Im_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Im_tr_paths[13];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_m1_n_2(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Im_tr_paths[9] + Mom_def_Im_tr_paths[10] + Mom_def_Im_tr_paths[11] + Mom_def_Im_tr_paths[12] + Mom_def_Im_tr_paths[13] + (-2.) * Mom_def_Im_tr_paths[14];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_1(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[0] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_2(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[3] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_3(double complex *op_out)
{
    *op_out = +(-I * 3.46410161513775459) * Mom_def_Re_tr_paths[6] + (+I * 3.46410161513775459) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_4(double complex *op_out)
{
    *op_out = +(-I * 1.73205080756887729) * Mom_def_Re_tr_paths[10] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[11] + (-I * 1.73205080756887729) * Mom_def_Re_tr_paths[12] + (+I * 1.73205080756887729) * Mom_def_Re_tr_paths[13];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_5(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[0] + (-4.) * Mom_def_Re_tr_paths[1] + (2.) * Mom_def_Re_tr_paths[2];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_6(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[3] + (-4.) * Mom_def_Re_tr_paths[4] + (2.) * Mom_def_Re_tr_paths[5];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_7(double complex *op_out)
{
    *op_out = +(2.) * Mom_def_Re_tr_paths[6] + (-4.) * Mom_def_Re_tr_paths[7] + (2.) * Mom_def_Re_tr_paths[8];
}

static void OP_oneTr_p_1_1_1_Ir_3_C_1_n_8(double complex *op_out)
{
    *op_out = +(-2.) * Mom_def_Re_tr_paths[9] + Mom_def_Re_tr_paths[10] + Mom_def_Re_tr_paths[11] + Mom_def_Re_tr_paths[12] + Mom_def_Re_tr_paths[13] + (-2.) * Mom_def_Re_tr_paths[14];
}

static int last_t = -10;
void request_space_paths_evaluation() { last_t = -10; }

void eval_time_momentum_glueball_paths(int t, int px, int py, int pz)
{
    int n_x, n_y, n_z, idx, in;
    double complex ce = 0.;
    if (path_storage == NULL)
    {
        c0 = cexp(I * PI * (4. / GLB_X));
        c1 = cexp(I * PI * (2. / GLB_X));
        c2 = cexp(I * PI * (-4. / GLB_X));
        c3 = cexp(I * PI * (-2. / GLB_X));
        path_storage = malloc(npaths * X * Y * Z * sizeof(double complex));
        Mom_def_Re_tr_paths = malloc(npaths * sizeof(double complex));
        Mom_def_Im_tr_paths = malloc(npaths * sizeof(double complex));
    }

    for (in = 0; in < npaths; in++)
    {
        Mom_def_Re_tr_paths[in] = 0.;
        Mom_def_Im_tr_paths[in] = 0.;
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
                    Mom_def_Re_tr_paths[0] += ce * creal(path_storage[0 + idx]);
                    Mom_def_Im_tr_paths[0] += I * ce * cimag(path_storage[0 + idx]);
                    path_storage[1 + idx] = path1(in);
                    Mom_def_Re_tr_paths[1] += ce * creal(path_storage[1 + idx]);
                    Mom_def_Im_tr_paths[1] += I * ce * cimag(path_storage[1 + idx]);
                    path_storage[2 + idx] = path2(in);
                    Mom_def_Re_tr_paths[2] += ce * creal(path_storage[2 + idx]);
                    Mom_def_Im_tr_paths[2] += I * ce * cimag(path_storage[2 + idx]);
                    path_storage[3 + idx] = path3(in);
                    Mom_def_Re_tr_paths[3] += ce * creal(path_storage[3 + idx]);
                    Mom_def_Im_tr_paths[3] += I * ce * cimag(path_storage[3 + idx]);
                    path_storage[4 + idx] = path4(in);
                    Mom_def_Re_tr_paths[4] += ce * creal(path_storage[4 + idx]);
                    Mom_def_Im_tr_paths[4] += I * ce * cimag(path_storage[4 + idx]);
                    path_storage[5 + idx] = path5(in);
                    Mom_def_Re_tr_paths[5] += ce * creal(path_storage[5 + idx]);
                    Mom_def_Im_tr_paths[5] += I * ce * cimag(path_storage[5 + idx]);
                    path_storage[6 + idx] = path6(in);
                    Mom_def_Re_tr_paths[6] += ce * creal(path_storage[6 + idx]);
                    Mom_def_Im_tr_paths[6] += I * ce * cimag(path_storage[6 + idx]);
                    path_storage[7 + idx] = path7(in);
                    Mom_def_Re_tr_paths[7] += ce * creal(path_storage[7 + idx]);
                    Mom_def_Im_tr_paths[7] += I * ce * cimag(path_storage[7 + idx]);
                    path_storage[8 + idx] = path8(in);
                    Mom_def_Re_tr_paths[8] += ce * creal(path_storage[8 + idx]);
                    Mom_def_Im_tr_paths[8] += I * ce * cimag(path_storage[8 + idx]);
                    path_storage[9 + idx] = path9(in);
                    Mom_def_Re_tr_paths[9] += ce * creal(path_storage[9 + idx]);
                    Mom_def_Im_tr_paths[9] += I * ce * cimag(path_storage[9 + idx]);
                    path_storage[10 + idx] = path10(in);
                    Mom_def_Re_tr_paths[10] += ce * creal(path_storage[10 + idx]);
                    Mom_def_Im_tr_paths[10] += I * ce * cimag(path_storage[10 + idx]);
                    path_storage[11 + idx] = path11(in);
                    Mom_def_Re_tr_paths[11] += ce * creal(path_storage[11 + idx]);
                    Mom_def_Im_tr_paths[11] += I * ce * cimag(path_storage[11 + idx]);
                    path_storage[12 + idx] = path12(in);
                    Mom_def_Re_tr_paths[12] += ce * creal(path_storage[12 + idx]);
                    Mom_def_Im_tr_paths[12] += I * ce * cimag(path_storage[12 + idx]);
                    path_storage[13 + idx] = path13(in);
                    Mom_def_Re_tr_paths[13] += ce * creal(path_storage[13 + idx]);
                    Mom_def_Im_tr_paths[13] += I * ce * cimag(path_storage[13 + idx]);
                    path_storage[14 + idx] = path14(in);
                    Mom_def_Re_tr_paths[14] += ce * creal(path_storage[14 + idx]);
                    Mom_def_Im_tr_paths[14] += I * ce * cimag(path_storage[14 + idx]);
                    path_storage[15 + idx] = path15(in);
                    Mom_def_Re_tr_paths[15] += ce * creal(path_storage[15 + idx]);
                    Mom_def_Im_tr_paths[15] += I * ce * cimag(path_storage[15 + idx]);
                    path_storage[16 + idx] = path16(in);
                    Mom_def_Re_tr_paths[16] += ce * creal(path_storage[16 + idx]);
                    Mom_def_Im_tr_paths[16] += I * ce * cimag(path_storage[16 + idx]);
                    path_storage[17 + idx] = path17(in);
                    Mom_def_Re_tr_paths[17] += ce * creal(path_storage[17 + idx]);
                    Mom_def_Im_tr_paths[17] += I * ce * cimag(path_storage[17 + idx]);
                    path_storage[18 + idx] = path18(in);
                    Mom_def_Re_tr_paths[18] += ce * creal(path_storage[18 + idx]);
                    Mom_def_Im_tr_paths[18] += I * ce * cimag(path_storage[18 + idx]);
                    path_storage[19 + idx] = path19(in);
                    Mom_def_Re_tr_paths[19] += ce * creal(path_storage[19 + idx]);
                    Mom_def_Im_tr_paths[19] += I * ce * cimag(path_storage[19 + idx]);
                    path_storage[20 + idx] = path20(in);
                    Mom_def_Re_tr_paths[20] += ce * creal(path_storage[20 + idx]);
                    Mom_def_Im_tr_paths[20] += I * ce * cimag(path_storage[20 + idx]);
                    path_storage[21 + idx] = path21(in);
                    Mom_def_Re_tr_paths[21] += ce * creal(path_storage[21 + idx]);
                    Mom_def_Im_tr_paths[21] += I * ce * cimag(path_storage[21 + idx]);
                    path_storage[22 + idx] = path22(in);
                    Mom_def_Re_tr_paths[22] += ce * creal(path_storage[22 + idx]);
                    Mom_def_Im_tr_paths[22] += I * ce * cimag(path_storage[22 + idx]);
                    path_storage[23 + idx] = path23(in);
                    Mom_def_Re_tr_paths[23] += ce * creal(path_storage[23 + idx]);
                    Mom_def_Im_tr_paths[23] += I * ce * cimag(path_storage[23 + idx]);
                    path_storage[24 + idx] = path24(in);
                    Mom_def_Re_tr_paths[24] += ce * creal(path_storage[24 + idx]);
                    Mom_def_Im_tr_paths[24] += I * ce * cimag(path_storage[24 + idx]);
                    path_storage[25 + idx] = path25(in);
                    Mom_def_Re_tr_paths[25] += ce * creal(path_storage[25 + idx]);
                    Mom_def_Im_tr_paths[25] += I * ce * cimag(path_storage[25 + idx]);
                    path_storage[26 + idx] = path26(in);
                    Mom_def_Re_tr_paths[26] += ce * creal(path_storage[26 + idx]);
                    Mom_def_Im_tr_paths[26] += I * ce * cimag(path_storage[26 + idx]);
                    path_storage[27 + idx] = path27(in);
                    Mom_def_Re_tr_paths[27] += ce * creal(path_storage[27 + idx]);
                    Mom_def_Im_tr_paths[27] += I * ce * cimag(path_storage[27 + idx]);
                    path_storage[28 + idx] = path28(in);
                    Mom_def_Re_tr_paths[28] += ce * creal(path_storage[28 + idx]);
                    Mom_def_Im_tr_paths[28] += I * ce * cimag(path_storage[28 + idx]);
                    path_storage[29 + idx] = path29(in);
                    Mom_def_Re_tr_paths[29] += ce * creal(path_storage[29 + idx]);
                    Mom_def_Im_tr_paths[29] += I * ce * cimag(path_storage[29 + idx]);
                    path_storage[30 + idx] = path30(in);
                    Mom_def_Re_tr_paths[30] += ce * creal(path_storage[30 + idx]);
                    Mom_def_Im_tr_paths[30] += I * ce * cimag(path_storage[30 + idx]);
                    path_storage[31 + idx] = path31(in);
                    Mom_def_Re_tr_paths[31] += ce * creal(path_storage[31 + idx]);
                    Mom_def_Im_tr_paths[31] += I * ce * cimag(path_storage[31 + idx]);
                    path_storage[32 + idx] = path32(in);
                    Mom_def_Re_tr_paths[32] += ce * creal(path_storage[32 + idx]);
                    Mom_def_Im_tr_paths[32] += I * ce * cimag(path_storage[32 + idx]);
                    path_storage[33 + idx] = path33(in);
                    Mom_def_Re_tr_paths[33] += ce * creal(path_storage[33 + idx]);
                    Mom_def_Im_tr_paths[33] += I * ce * cimag(path_storage[33 + idx]);
                    path_storage[34 + idx] = path34(in);
                    Mom_def_Re_tr_paths[34] += ce * creal(path_storage[34 + idx]);
                    Mom_def_Im_tr_paths[34] += I * ce * cimag(path_storage[34 + idx]);
                    path_storage[35 + idx] = path35(in);
                    Mom_def_Re_tr_paths[35] += ce * creal(path_storage[35 + idx]);
                    Mom_def_Im_tr_paths[35] += I * ce * cimag(path_storage[35 + idx]);
                    path_storage[36 + idx] = path36(in);
                    Mom_def_Re_tr_paths[36] += ce * creal(path_storage[36 + idx]);
                    Mom_def_Im_tr_paths[36] += I * ce * cimag(path_storage[36 + idx]);
                    path_storage[37 + idx] = path37(in);
                    Mom_def_Re_tr_paths[37] += ce * creal(path_storage[37 + idx]);
                    Mom_def_Im_tr_paths[37] += I * ce * cimag(path_storage[37 + idx]);
                    path_storage[38 + idx] = path38(in);
                    Mom_def_Re_tr_paths[38] += ce * creal(path_storage[38 + idx]);
                    Mom_def_Im_tr_paths[38] += I * ce * cimag(path_storage[38 + idx]);
                    path_storage[39 + idx] = path39(in);
                    Mom_def_Re_tr_paths[39] += ce * creal(path_storage[39 + idx]);
                    Mom_def_Im_tr_paths[39] += I * ce * cimag(path_storage[39 + idx]);
                    path_storage[40 + idx] = path40(in);
                    Mom_def_Re_tr_paths[40] += ce * creal(path_storage[40 + idx]);
                    Mom_def_Im_tr_paths[40] += I * ce * cimag(path_storage[40 + idx]);
                    path_storage[41 + idx] = path41(in);
                    Mom_def_Re_tr_paths[41] += ce * creal(path_storage[41 + idx]);
                    Mom_def_Im_tr_paths[41] += I * ce * cimag(path_storage[41 + idx]);
                    path_storage[42 + idx] = path42(in);
                    Mom_def_Re_tr_paths[42] += ce * creal(path_storage[42 + idx]);
                    Mom_def_Im_tr_paths[42] += I * ce * cimag(path_storage[42 + idx]);
                    path_storage[43 + idx] = path43(in);
                    Mom_def_Re_tr_paths[43] += ce * creal(path_storage[43 + idx]);
                    Mom_def_Im_tr_paths[43] += I * ce * cimag(path_storage[43 + idx]);
                    path_storage[44 + idx] = path44(in);
                    Mom_def_Re_tr_paths[44] += ce * creal(path_storage[44 + idx]);
                    Mom_def_Im_tr_paths[44] += I * ce * cimag(path_storage[44 + idx]);
                    path_storage[45 + idx] = path45(in);
                    Mom_def_Re_tr_paths[45] += ce * creal(path_storage[45 + idx]);
                    Mom_def_Im_tr_paths[45] += I * ce * cimag(path_storage[45 + idx]);
                    path_storage[46 + idx] = path46(in);
                    Mom_def_Re_tr_paths[46] += ce * creal(path_storage[46 + idx]);
                    Mom_def_Im_tr_paths[46] += I * ce * cimag(path_storage[46 + idx]);
                    path_storage[47 + idx] = path47(in);
                    Mom_def_Re_tr_paths[47] += ce * creal(path_storage[47 + idx]);
                    Mom_def_Im_tr_paths[47] += I * ce * cimag(path_storage[47 + idx]);
                    path_storage[48 + idx] = path48(in);
                    Mom_def_Re_tr_paths[48] += ce * creal(path_storage[48 + idx]);
                    Mom_def_Im_tr_paths[48] += I * ce * cimag(path_storage[48 + idx]);
                    path_storage[49 + idx] = path49(in);
                    Mom_def_Re_tr_paths[49] += ce * creal(path_storage[49 + idx]);
                    Mom_def_Im_tr_paths[49] += I * ce * cimag(path_storage[49 + idx]);
                    path_storage[50 + idx] = path50(in);
                    Mom_def_Re_tr_paths[50] += ce * creal(path_storage[50 + idx]);
                    Mom_def_Im_tr_paths[50] += I * ce * cimag(path_storage[50 + idx]);
                    path_storage[51 + idx] = path51(in);
                    Mom_def_Re_tr_paths[51] += ce * creal(path_storage[51 + idx]);
                    Mom_def_Im_tr_paths[51] += I * ce * cimag(path_storage[51 + idx]);
                    path_storage[52 + idx] = path52(in);
                    Mom_def_Re_tr_paths[52] += ce * creal(path_storage[52 + idx]);
                    Mom_def_Im_tr_paths[52] += I * ce * cimag(path_storage[52 + idx]);
                    path_storage[53 + idx] = path53(in);
                    Mom_def_Re_tr_paths[53] += ce * creal(path_storage[53 + idx]);
                    Mom_def_Im_tr_paths[53] += I * ce * cimag(path_storage[53 + idx]);
                    path_storage[54 + idx] = path54(in);
                    Mom_def_Re_tr_paths[54] += ce * creal(path_storage[54 + idx]);
                    Mom_def_Im_tr_paths[54] += I * ce * cimag(path_storage[54 + idx]);
                    path_storage[55 + idx] = path55(in);
                    Mom_def_Re_tr_paths[55] += ce * creal(path_storage[55 + idx]);
                    Mom_def_Im_tr_paths[55] += I * ce * cimag(path_storage[55 + idx]);
                    path_storage[56 + idx] = path56(in);
                    Mom_def_Re_tr_paths[56] += ce * creal(path_storage[56 + idx]);
                    Mom_def_Im_tr_paths[56] += I * ce * cimag(path_storage[56 + idx]);
                    path_storage[57 + idx] = path57(in);
                    Mom_def_Re_tr_paths[57] += ce * creal(path_storage[57 + idx]);
                    Mom_def_Im_tr_paths[57] += I * ce * cimag(path_storage[57 + idx]);
                    path_storage[58 + idx] = path58(in);
                    Mom_def_Re_tr_paths[58] += ce * creal(path_storage[58 + idx]);
                    Mom_def_Im_tr_paths[58] += I * ce * cimag(path_storage[58 + idx]);
                    path_storage[59 + idx] = path59(in);
                    Mom_def_Re_tr_paths[59] += ce * creal(path_storage[59 + idx]);
                    Mom_def_Im_tr_paths[59] += I * ce * cimag(path_storage[59 + idx]);
                    path_storage[60 + idx] = path60(in);
                    Mom_def_Re_tr_paths[60] += ce * creal(path_storage[60 + idx]);
                    Mom_def_Im_tr_paths[60] += I * ce * cimag(path_storage[60 + idx]);
                    path_storage[61 + idx] = path61(in);
                    Mom_def_Re_tr_paths[61] += ce * creal(path_storage[61 + idx]);
                    Mom_def_Im_tr_paths[61] += I * ce * cimag(path_storage[61 + idx]);
                    path_storage[62 + idx] = path62(in);
                    Mom_def_Re_tr_paths[62] += ce * creal(path_storage[62 + idx]);
                    Mom_def_Im_tr_paths[62] += I * ce * cimag(path_storage[62 + idx]);
                    path_storage[63 + idx] = path63(in);
                    Mom_def_Re_tr_paths[63] += ce * creal(path_storage[63 + idx]);
                    Mom_def_Im_tr_paths[63] += I * ce * cimag(path_storage[63 + idx]);
                    path_storage[64 + idx] = path64(in);
                    Mom_def_Re_tr_paths[64] += ce * creal(path_storage[64 + idx]);
                    Mom_def_Im_tr_paths[64] += I * ce * cimag(path_storage[64 + idx]);
                    path_storage[65 + idx] = path65(in);
                    Mom_def_Re_tr_paths[65] += ce * creal(path_storage[65 + idx]);
                    Mom_def_Im_tr_paths[65] += I * ce * cimag(path_storage[65 + idx]);
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
                    Mom_def_Re_tr_paths[0] += ce * creal(path_storage[0 + idx]);
                    Mom_def_Im_tr_paths[0] += I * ce * cimag(path_storage[0 + idx]);
                    Mom_def_Re_tr_paths[1] += ce * creal(path_storage[1 + idx]);
                    Mom_def_Im_tr_paths[1] += I * ce * cimag(path_storage[1 + idx]);
                    Mom_def_Re_tr_paths[2] += ce * creal(path_storage[2 + idx]);
                    Mom_def_Im_tr_paths[2] += I * ce * cimag(path_storage[2 + idx]);
                    Mom_def_Re_tr_paths[3] += ce * creal(path_storage[3 + idx]);
                    Mom_def_Im_tr_paths[3] += I * ce * cimag(path_storage[3 + idx]);
                    Mom_def_Re_tr_paths[4] += ce * creal(path_storage[4 + idx]);
                    Mom_def_Im_tr_paths[4] += I * ce * cimag(path_storage[4 + idx]);
                    Mom_def_Re_tr_paths[5] += ce * creal(path_storage[5 + idx]);
                    Mom_def_Im_tr_paths[5] += I * ce * cimag(path_storage[5 + idx]);
                    Mom_def_Re_tr_paths[6] += ce * creal(path_storage[6 + idx]);
                    Mom_def_Im_tr_paths[6] += I * ce * cimag(path_storage[6 + idx]);
                    Mom_def_Re_tr_paths[7] += ce * creal(path_storage[7 + idx]);
                    Mom_def_Im_tr_paths[7] += I * ce * cimag(path_storage[7 + idx]);
                    Mom_def_Re_tr_paths[8] += ce * creal(path_storage[8 + idx]);
                    Mom_def_Im_tr_paths[8] += I * ce * cimag(path_storage[8 + idx]);
                    Mom_def_Re_tr_paths[9] += ce * creal(path_storage[9 + idx]);
                    Mom_def_Im_tr_paths[9] += I * ce * cimag(path_storage[9 + idx]);
                    Mom_def_Re_tr_paths[10] += ce * creal(path_storage[10 + idx]);
                    Mom_def_Im_tr_paths[10] += I * ce * cimag(path_storage[10 + idx]);
                    Mom_def_Re_tr_paths[11] += ce * creal(path_storage[11 + idx]);
                    Mom_def_Im_tr_paths[11] += I * ce * cimag(path_storage[11 + idx]);
                    Mom_def_Re_tr_paths[12] += ce * creal(path_storage[12 + idx]);
                    Mom_def_Im_tr_paths[12] += I * ce * cimag(path_storage[12 + idx]);
                    Mom_def_Re_tr_paths[13] += ce * creal(path_storage[13 + idx]);
                    Mom_def_Im_tr_paths[13] += I * ce * cimag(path_storage[13 + idx]);
                    Mom_def_Re_tr_paths[14] += ce * creal(path_storage[14 + idx]);
                    Mom_def_Im_tr_paths[14] += I * ce * cimag(path_storage[14 + idx]);
                    Mom_def_Re_tr_paths[15] += ce * creal(path_storage[15 + idx]);
                    Mom_def_Im_tr_paths[15] += I * ce * cimag(path_storage[15 + idx]);
                    Mom_def_Re_tr_paths[16] += ce * creal(path_storage[16 + idx]);
                    Mom_def_Im_tr_paths[16] += I * ce * cimag(path_storage[16 + idx]);
                    Mom_def_Re_tr_paths[17] += ce * creal(path_storage[17 + idx]);
                    Mom_def_Im_tr_paths[17] += I * ce * cimag(path_storage[17 + idx]);
                    Mom_def_Re_tr_paths[18] += ce * creal(path_storage[18 + idx]);
                    Mom_def_Im_tr_paths[18] += I * ce * cimag(path_storage[18 + idx]);
                    Mom_def_Re_tr_paths[19] += ce * creal(path_storage[19 + idx]);
                    Mom_def_Im_tr_paths[19] += I * ce * cimag(path_storage[19 + idx]);
                    Mom_def_Re_tr_paths[20] += ce * creal(path_storage[20 + idx]);
                    Mom_def_Im_tr_paths[20] += I * ce * cimag(path_storage[20 + idx]);
                    Mom_def_Re_tr_paths[21] += ce * creal(path_storage[21 + idx]);
                    Mom_def_Im_tr_paths[21] += I * ce * cimag(path_storage[21 + idx]);
                    Mom_def_Re_tr_paths[22] += ce * creal(path_storage[22 + idx]);
                    Mom_def_Im_tr_paths[22] += I * ce * cimag(path_storage[22 + idx]);
                    Mom_def_Re_tr_paths[23] += ce * creal(path_storage[23 + idx]);
                    Mom_def_Im_tr_paths[23] += I * ce * cimag(path_storage[23 + idx]);
                    Mom_def_Re_tr_paths[24] += ce * creal(path_storage[24 + idx]);
                    Mom_def_Im_tr_paths[24] += I * ce * cimag(path_storage[24 + idx]);
                    Mom_def_Re_tr_paths[25] += ce * creal(path_storage[25 + idx]);
                    Mom_def_Im_tr_paths[25] += I * ce * cimag(path_storage[25 + idx]);
                    Mom_def_Re_tr_paths[26] += ce * creal(path_storage[26 + idx]);
                    Mom_def_Im_tr_paths[26] += I * ce * cimag(path_storage[26 + idx]);
                    Mom_def_Re_tr_paths[27] += ce * creal(path_storage[27 + idx]);
                    Mom_def_Im_tr_paths[27] += I * ce * cimag(path_storage[27 + idx]);
                    Mom_def_Re_tr_paths[28] += ce * creal(path_storage[28 + idx]);
                    Mom_def_Im_tr_paths[28] += I * ce * cimag(path_storage[28 + idx]);
                    Mom_def_Re_tr_paths[29] += ce * creal(path_storage[29 + idx]);
                    Mom_def_Im_tr_paths[29] += I * ce * cimag(path_storage[29 + idx]);
                    Mom_def_Re_tr_paths[30] += ce * creal(path_storage[30 + idx]);
                    Mom_def_Im_tr_paths[30] += I * ce * cimag(path_storage[30 + idx]);
                    Mom_def_Re_tr_paths[31] += ce * creal(path_storage[31 + idx]);
                    Mom_def_Im_tr_paths[31] += I * ce * cimag(path_storage[31 + idx]);
                    Mom_def_Re_tr_paths[32] += ce * creal(path_storage[32 + idx]);
                    Mom_def_Im_tr_paths[32] += I * ce * cimag(path_storage[32 + idx]);
                    Mom_def_Re_tr_paths[33] += ce * creal(path_storage[33 + idx]);
                    Mom_def_Im_tr_paths[33] += I * ce * cimag(path_storage[33 + idx]);
                    Mom_def_Re_tr_paths[34] += ce * creal(path_storage[34 + idx]);
                    Mom_def_Im_tr_paths[34] += I * ce * cimag(path_storage[34 + idx]);
                    Mom_def_Re_tr_paths[35] += ce * creal(path_storage[35 + idx]);
                    Mom_def_Im_tr_paths[35] += I * ce * cimag(path_storage[35 + idx]);
                    Mom_def_Re_tr_paths[36] += ce * creal(path_storage[36 + idx]);
                    Mom_def_Im_tr_paths[36] += I * ce * cimag(path_storage[36 + idx]);
                    Mom_def_Re_tr_paths[37] += ce * creal(path_storage[37 + idx]);
                    Mom_def_Im_tr_paths[37] += I * ce * cimag(path_storage[37 + idx]);
                    Mom_def_Re_tr_paths[38] += ce * creal(path_storage[38 + idx]);
                    Mom_def_Im_tr_paths[38] += I * ce * cimag(path_storage[38 + idx]);
                    Mom_def_Re_tr_paths[39] += ce * creal(path_storage[39 + idx]);
                    Mom_def_Im_tr_paths[39] += I * ce * cimag(path_storage[39 + idx]);
                    Mom_def_Re_tr_paths[40] += ce * creal(path_storage[40 + idx]);
                    Mom_def_Im_tr_paths[40] += I * ce * cimag(path_storage[40 + idx]);
                    Mom_def_Re_tr_paths[41] += ce * creal(path_storage[41 + idx]);
                    Mom_def_Im_tr_paths[41] += I * ce * cimag(path_storage[41 + idx]);
                    Mom_def_Re_tr_paths[42] += ce * creal(path_storage[42 + idx]);
                    Mom_def_Im_tr_paths[42] += I * ce * cimag(path_storage[42 + idx]);
                    Mom_def_Re_tr_paths[43] += ce * creal(path_storage[43 + idx]);
                    Mom_def_Im_tr_paths[43] += I * ce * cimag(path_storage[43 + idx]);
                    Mom_def_Re_tr_paths[44] += ce * creal(path_storage[44 + idx]);
                    Mom_def_Im_tr_paths[44] += I * ce * cimag(path_storage[44 + idx]);
                    Mom_def_Re_tr_paths[45] += ce * creal(path_storage[45 + idx]);
                    Mom_def_Im_tr_paths[45] += I * ce * cimag(path_storage[45 + idx]);
                    Mom_def_Re_tr_paths[46] += ce * creal(path_storage[46 + idx]);
                    Mom_def_Im_tr_paths[46] += I * ce * cimag(path_storage[46 + idx]);
                    Mom_def_Re_tr_paths[47] += ce * creal(path_storage[47 + idx]);
                    Mom_def_Im_tr_paths[47] += I * ce * cimag(path_storage[47 + idx]);
                    Mom_def_Re_tr_paths[48] += ce * creal(path_storage[48 + idx]);
                    Mom_def_Im_tr_paths[48] += I * ce * cimag(path_storage[48 + idx]);
                    Mom_def_Re_tr_paths[49] += ce * creal(path_storage[49 + idx]);
                    Mom_def_Im_tr_paths[49] += I * ce * cimag(path_storage[49 + idx]);
                    Mom_def_Re_tr_paths[50] += ce * creal(path_storage[50 + idx]);
                    Mom_def_Im_tr_paths[50] += I * ce * cimag(path_storage[50 + idx]);
                    Mom_def_Re_tr_paths[51] += ce * creal(path_storage[51 + idx]);
                    Mom_def_Im_tr_paths[51] += I * ce * cimag(path_storage[51 + idx]);
                    Mom_def_Re_tr_paths[52] += ce * creal(path_storage[52 + idx]);
                    Mom_def_Im_tr_paths[52] += I * ce * cimag(path_storage[52 + idx]);
                    Mom_def_Re_tr_paths[53] += ce * creal(path_storage[53 + idx]);
                    Mom_def_Im_tr_paths[53] += I * ce * cimag(path_storage[53 + idx]);
                    Mom_def_Re_tr_paths[54] += ce * creal(path_storage[54 + idx]);
                    Mom_def_Im_tr_paths[54] += I * ce * cimag(path_storage[54 + idx]);
                    Mom_def_Re_tr_paths[55] += ce * creal(path_storage[55 + idx]);
                    Mom_def_Im_tr_paths[55] += I * ce * cimag(path_storage[55 + idx]);
                    Mom_def_Re_tr_paths[56] += ce * creal(path_storage[56 + idx]);
                    Mom_def_Im_tr_paths[56] += I * ce * cimag(path_storage[56 + idx]);
                    Mom_def_Re_tr_paths[57] += ce * creal(path_storage[57 + idx]);
                    Mom_def_Im_tr_paths[57] += I * ce * cimag(path_storage[57 + idx]);
                    Mom_def_Re_tr_paths[58] += ce * creal(path_storage[58 + idx]);
                    Mom_def_Im_tr_paths[58] += I * ce * cimag(path_storage[58 + idx]);
                    Mom_def_Re_tr_paths[59] += ce * creal(path_storage[59 + idx]);
                    Mom_def_Im_tr_paths[59] += I * ce * cimag(path_storage[59 + idx]);
                    Mom_def_Re_tr_paths[60] += ce * creal(path_storage[60 + idx]);
                    Mom_def_Im_tr_paths[60] += I * ce * cimag(path_storage[60 + idx]);
                    Mom_def_Re_tr_paths[61] += ce * creal(path_storage[61 + idx]);
                    Mom_def_Im_tr_paths[61] += I * ce * cimag(path_storage[61 + idx]);
                    Mom_def_Re_tr_paths[62] += ce * creal(path_storage[62 + idx]);
                    Mom_def_Im_tr_paths[62] += I * ce * cimag(path_storage[62 + idx]);
                    Mom_def_Re_tr_paths[63] += ce * creal(path_storage[63 + idx]);
                    Mom_def_Im_tr_paths[63] += I * ce * cimag(path_storage[63 + idx]);
                    Mom_def_Re_tr_paths[64] += ce * creal(path_storage[64 + idx]);
                    Mom_def_Im_tr_paths[64] += I * ce * cimag(path_storage[64 + idx]);
                    Mom_def_Re_tr_paths[65] += ce * creal(path_storage[65 + idx]);
                    Mom_def_Im_tr_paths[65] += I * ce * cimag(path_storage[65 + idx]);
                }
    }
};
void eval_all_glueball_ops(int t, double complex *numerical_op)
{
    eval_time_momentum_glueball_paths(t, -1, -1, -1);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_1(numerical_op + 1);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_2(numerical_op + 2);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_3(numerical_op + 3);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_4(numerical_op + 4);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_1(numerical_op + 5);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_2(numerical_op + 6);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_3(numerical_op + 7);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_4(numerical_op + 8);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_1_n_1(numerical_op + 9);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_m1_n_1(numerical_op + 10);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_m1_n_2(numerical_op + 11);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_1(numerical_op + 12);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_2(numerical_op + 13);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_3(numerical_op + 14);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_4(numerical_op + 15);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_5(numerical_op + 16);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_6(numerical_op + 17);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_7(numerical_op + 18);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_8(numerical_op + 19);
    eval_time_momentum_glueball_paths(t, -1, -1, 0);
    OP_oneTr_p_m1_m1_0_Ir_1_C_m1_n_1(numerical_op + 20);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_1(numerical_op + 21);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_2(numerical_op + 22);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_3(numerical_op + 23);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_4(numerical_op + 24);
    OP_oneTr_p_m1_m1_0_Ir_2_C_m1_n_1(numerical_op + 25);
    OP_oneTr_p_m1_m1_0_Ir_2_C_1_n_1(numerical_op + 26);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_1(numerical_op + 27);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_2(numerical_op + 28);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_3(numerical_op + 29);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_4(numerical_op + 30);
    OP_oneTr_p_m1_m1_0_Ir_3_C_1_n_1(numerical_op + 31);
    OP_oneTr_p_m1_m1_0_Ir_4_C_m1_n_1(numerical_op + 32);
    OP_oneTr_p_m1_m1_0_Ir_4_C_1_n_1(numerical_op + 33);
    eval_time_momentum_glueball_paths(t, -1, -1, 1);
    OP_oneTr_p_m1_m1_1_Ir_1_C_m1_n_1(numerical_op + 34);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_1(numerical_op + 35);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_2(numerical_op + 36);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_3(numerical_op + 37);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_4(numerical_op + 38);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_1(numerical_op + 39);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_2(numerical_op + 40);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_3(numerical_op + 41);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_4(numerical_op + 42);
    OP_oneTr_p_m1_m1_1_Ir_2_C_1_n_1(numerical_op + 43);
    OP_oneTr_p_m1_m1_1_Ir_3_C_m1_n_1(numerical_op + 44);
    OP_oneTr_p_m1_m1_1_Ir_3_C_m1_n_2(numerical_op + 45);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_1(numerical_op + 46);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_2(numerical_op + 47);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_3(numerical_op + 48);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_4(numerical_op + 49);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_5(numerical_op + 50);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_6(numerical_op + 51);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_7(numerical_op + 52);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_8(numerical_op + 53);
    eval_time_momentum_glueball_paths(t, -1, 0, -1);
    OP_oneTr_p_m1_0_m1_Ir_1_C_m1_n_1(numerical_op + 54);
    OP_oneTr_p_m1_0_m1_Ir_1_C_m1_n_2(numerical_op + 55);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_1(numerical_op + 56);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_2(numerical_op + 57);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_3(numerical_op + 58);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_4(numerical_op + 59);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_1(numerical_op + 60);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_2(numerical_op + 61);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_3(numerical_op + 62);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_4(numerical_op + 63);
    OP_oneTr_p_m1_0_m1_Ir_2_C_1_n_1(numerical_op + 64);
    OP_oneTr_p_m1_0_m1_Ir_2_C_1_n_2(numerical_op + 65);
    OP_oneTr_p_m1_0_m1_Ir_3_C_m1_n_1(numerical_op + 66);
    OP_oneTr_p_m1_0_m1_Ir_3_C_m1_n_2(numerical_op + 67);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_1(numerical_op + 68);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_2(numerical_op + 69);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_3(numerical_op + 70);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_4(numerical_op + 71);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_1(numerical_op + 72);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_2(numerical_op + 73);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_3(numerical_op + 74);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_4(numerical_op + 75);
    OP_oneTr_p_m1_0_m1_Ir_4_C_1_n_1(numerical_op + 76);
    OP_oneTr_p_m1_0_m1_Ir_4_C_1_n_2(numerical_op + 77);
    eval_time_momentum_glueball_paths(t, -1, 0, 0);
    OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_1(numerical_op + 78);
    OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_2(numerical_op + 79);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_1(numerical_op + 80);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_2(numerical_op + 81);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_3(numerical_op + 82);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_4(numerical_op + 83);
    OP_oneTr_p_m1_0_0_Ir_2_C_m1_n_1(numerical_op + 84);
    OP_oneTr_p_m1_0_0_Ir_2_C_1_n_1(numerical_op + 85);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_1(numerical_op + 86);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_2(numerical_op + 87);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_3(numerical_op + 88);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_4(numerical_op + 89);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_5(numerical_op + 90);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_6(numerical_op + 91);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_7(numerical_op + 92);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_8(numerical_op + 93);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_1(numerical_op + 94);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_2(numerical_op + 95);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_3(numerical_op + 96);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_4(numerical_op + 97);
    OP_oneTr_p_m1_0_0_Ir_4_C_m1_n_1(numerical_op + 98);
    OP_oneTr_p_m1_0_0_Ir_4_C_m1_n_2(numerical_op + 99);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_1(numerical_op + 100);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_2(numerical_op + 101);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_3(numerical_op + 102);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_4(numerical_op + 103);
    OP_oneTr_p_m1_0_0_Ir_5_C_m1_n_1(numerical_op + 104);
    OP_oneTr_p_m1_0_0_Ir_5_C_1_n_1(numerical_op + 105);
    eval_time_momentum_glueball_paths(t, -1, 0, 1);
    OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_1(numerical_op + 106);
    OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_2(numerical_op + 107);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_1(numerical_op + 108);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_2(numerical_op + 109);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_3(numerical_op + 110);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_4(numerical_op + 111);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_1(numerical_op + 112);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_2(numerical_op + 113);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_3(numerical_op + 114);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_4(numerical_op + 115);
    OP_oneTr_p_m1_0_1_Ir_2_C_1_n_1(numerical_op + 116);
    OP_oneTr_p_m1_0_1_Ir_2_C_1_n_2(numerical_op + 117);
    OP_oneTr_p_m1_0_1_Ir_3_C_m1_n_1(numerical_op + 118);
    OP_oneTr_p_m1_0_1_Ir_3_C_m1_n_2(numerical_op + 119);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_1(numerical_op + 120);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_2(numerical_op + 121);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_3(numerical_op + 122);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_4(numerical_op + 123);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_1(numerical_op + 124);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_2(numerical_op + 125);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_3(numerical_op + 126);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_4(numerical_op + 127);
    OP_oneTr_p_m1_0_1_Ir_4_C_1_n_1(numerical_op + 128);
    OP_oneTr_p_m1_0_1_Ir_4_C_1_n_2(numerical_op + 129);
    eval_time_momentum_glueball_paths(t, -1, 1, -1);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_1(numerical_op + 130);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_2(numerical_op + 131);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_3(numerical_op + 132);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_4(numerical_op + 133);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_1(numerical_op + 134);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_2(numerical_op + 135);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_3(numerical_op + 136);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_4(numerical_op + 137);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_1(numerical_op + 138);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_2(numerical_op + 139);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_3(numerical_op + 140);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_4(numerical_op + 141);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_1(numerical_op + 142);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_2(numerical_op + 143);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_3(numerical_op + 144);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_4(numerical_op + 145);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_1(numerical_op + 146);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_2(numerical_op + 147);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_3(numerical_op + 148);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_4(numerical_op + 149);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_5(numerical_op + 150);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_6(numerical_op + 151);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_7(numerical_op + 152);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_8(numerical_op + 153);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_1(numerical_op + 154);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_2(numerical_op + 155);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_3(numerical_op + 156);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_4(numerical_op + 157);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_5(numerical_op + 158);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_6(numerical_op + 159);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_7(numerical_op + 160);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_8(numerical_op + 161);
    eval_time_momentum_glueball_paths(t, -1, 1, 0);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_1(numerical_op + 162);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_2(numerical_op + 163);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_3(numerical_op + 164);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_4(numerical_op + 165);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_1(numerical_op + 166);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_2(numerical_op + 167);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_3(numerical_op + 168);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_4(numerical_op + 169);
    OP_oneTr_p_m1_1_0_Ir_2_C_m1_n_1(numerical_op + 170);
    OP_oneTr_p_m1_1_0_Ir_2_C_1_n_1(numerical_op + 171);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_1(numerical_op + 172);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_2(numerical_op + 173);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_3(numerical_op + 174);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_4(numerical_op + 175);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_1(numerical_op + 176);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_2(numerical_op + 177);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_3(numerical_op + 178);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_4(numerical_op + 179);
    OP_oneTr_p_m1_1_0_Ir_4_C_m1_n_1(numerical_op + 180);
    OP_oneTr_p_m1_1_0_Ir_4_C_1_n_1(numerical_op + 181);
    eval_time_momentum_glueball_paths(t, -1, 1, 1);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_1(numerical_op + 182);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_2(numerical_op + 183);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_3(numerical_op + 184);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_4(numerical_op + 185);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_1(numerical_op + 186);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_2(numerical_op + 187);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_3(numerical_op + 188);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_4(numerical_op + 189);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_1(numerical_op + 190);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_2(numerical_op + 191);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_3(numerical_op + 192);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_4(numerical_op + 193);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_1(numerical_op + 194);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_2(numerical_op + 195);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_3(numerical_op + 196);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_4(numerical_op + 197);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_1(numerical_op + 198);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_2(numerical_op + 199);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_3(numerical_op + 200);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_4(numerical_op + 201);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_5(numerical_op + 202);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_6(numerical_op + 203);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_7(numerical_op + 204);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_8(numerical_op + 205);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_1(numerical_op + 206);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_2(numerical_op + 207);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_3(numerical_op + 208);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_4(numerical_op + 209);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_5(numerical_op + 210);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_6(numerical_op + 211);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_7(numerical_op + 212);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_8(numerical_op + 213);
    eval_time_momentum_glueball_paths(t, 0, -1, -1);
    OP_oneTr_p_0_m1_m1_Ir_1_C_m1_n_1(numerical_op + 214);
    OP_oneTr_p_0_m1_m1_Ir_1_C_m1_n_2(numerical_op + 215);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_1(numerical_op + 216);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_2(numerical_op + 217);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_3(numerical_op + 218);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_4(numerical_op + 219);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_1(numerical_op + 220);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_2(numerical_op + 221);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_3(numerical_op + 222);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_4(numerical_op + 223);
    OP_oneTr_p_0_m1_m1_Ir_2_C_1_n_1(numerical_op + 224);
    OP_oneTr_p_0_m1_m1_Ir_2_C_1_n_2(numerical_op + 225);
    OP_oneTr_p_0_m1_m1_Ir_3_C_m1_n_1(numerical_op + 226);
    OP_oneTr_p_0_m1_m1_Ir_3_C_m1_n_2(numerical_op + 227);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_1(numerical_op + 228);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_2(numerical_op + 229);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_3(numerical_op + 230);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_4(numerical_op + 231);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_1(numerical_op + 232);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_2(numerical_op + 233);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_3(numerical_op + 234);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_4(numerical_op + 235);
    OP_oneTr_p_0_m1_m1_Ir_4_C_1_n_1(numerical_op + 236);
    OP_oneTr_p_0_m1_m1_Ir_4_C_1_n_2(numerical_op + 237);
    eval_time_momentum_glueball_paths(t, 0, -1, 0);
    OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_1(numerical_op + 238);
    OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_2(numerical_op + 239);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_1(numerical_op + 240);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_2(numerical_op + 241);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_3(numerical_op + 242);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_4(numerical_op + 243);
    OP_oneTr_p_0_m1_0_Ir_2_C_m1_n_1(numerical_op + 244);
    OP_oneTr_p_0_m1_0_Ir_2_C_1_n_1(numerical_op + 245);
    OP_oneTr_p_0_m1_0_Ir_3_C_m1_n_1(numerical_op + 246);
    OP_oneTr_p_0_m1_0_Ir_3_C_m1_n_2(numerical_op + 247);
    OP_oneTr_p_0_m1_0_Ir_3_C_1_n_1(numerical_op + 248);
    OP_oneTr_p_0_m1_0_Ir_3_C_1_n_2(numerical_op + 249);
    OP_oneTr_p_0_m1_0_Ir_4_C_m1_n_1(numerical_op + 250);
    OP_oneTr_p_0_m1_0_Ir_4_C_m1_n_2(numerical_op + 251);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_1(numerical_op + 252);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_2(numerical_op + 253);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_3(numerical_op + 254);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_4(numerical_op + 255);
    OP_oneTr_p_0_m1_0_Ir_5_C_m1_n_1(numerical_op + 256);
    OP_oneTr_p_0_m1_0_Ir_5_C_1_n_1(numerical_op + 257);
    eval_time_momentum_glueball_paths(t, 0, -1, 1);
    OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_1(numerical_op + 258);
    OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_2(numerical_op + 259);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_1(numerical_op + 260);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_2(numerical_op + 261);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_3(numerical_op + 262);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_4(numerical_op + 263);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_1(numerical_op + 264);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_2(numerical_op + 265);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_3(numerical_op + 266);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_4(numerical_op + 267);
    OP_oneTr_p_0_m1_1_Ir_2_C_1_n_1(numerical_op + 268);
    OP_oneTr_p_0_m1_1_Ir_2_C_1_n_2(numerical_op + 269);
    OP_oneTr_p_0_m1_1_Ir_3_C_m1_n_1(numerical_op + 270);
    OP_oneTr_p_0_m1_1_Ir_3_C_m1_n_2(numerical_op + 271);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_1(numerical_op + 272);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_2(numerical_op + 273);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_3(numerical_op + 274);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_4(numerical_op + 275);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_1(numerical_op + 276);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_2(numerical_op + 277);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_3(numerical_op + 278);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_4(numerical_op + 279);
    OP_oneTr_p_0_m1_1_Ir_4_C_1_n_1(numerical_op + 280);
    OP_oneTr_p_0_m1_1_Ir_4_C_1_n_2(numerical_op + 281);
    eval_time_momentum_glueball_paths(t, 0, 0, -1);
    OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_1(numerical_op + 282);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_1(numerical_op + 283);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_2(numerical_op + 284);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_3(numerical_op + 285);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_4(numerical_op + 286);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_1(numerical_op + 287);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_2(numerical_op + 288);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_3(numerical_op + 289);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_4(numerical_op + 290);
    OP_oneTr_p_0_0_m1_Ir_2_C_1_n_1(numerical_op + 291);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_1(numerical_op + 292);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_2(numerical_op + 293);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_3(numerical_op + 294);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_4(numerical_op + 295);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_1(numerical_op + 296);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_2(numerical_op + 297);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_3(numerical_op + 298);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_4(numerical_op + 299);
    OP_oneTr_p_0_0_m1_Ir_4_C_m1_n_1(numerical_op + 300);
    OP_oneTr_p_0_0_m1_Ir_4_C_m1_n_2(numerical_op + 301);
    OP_oneTr_p_0_0_m1_Ir_4_C_1_n_1(numerical_op + 302);
    OP_oneTr_p_0_0_m1_Ir_5_C_m1_n_1(numerical_op + 303);
    OP_oneTr_p_0_0_m1_Ir_5_C_1_n_1(numerical_op + 304);
    OP_oneTr_p_0_0_m1_Ir_5_C_1_n_2(numerical_op + 305);
    eval_time_momentum_glueball_paths(t, 0, 0, 0);
    OP_oneTr_p_0_0_0_Ir_1_C_m1_n_1(numerical_op + 306);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_1(numerical_op + 307);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_2(numerical_op + 308);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_3(numerical_op + 309);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_4(numerical_op + 310);
    OP_oneTr_p_0_0_0_Ir_2_C_m1_n_1(numerical_op + 311);
    OP_oneTr_p_0_0_0_Ir_2_C_m1_n_2(numerical_op + 312);
    OP_oneTr_p_0_0_0_Ir_2_C_1_n_1(numerical_op + 313);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_1(numerical_op + 314);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_2(numerical_op + 315);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_3(numerical_op + 316);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_4(numerical_op + 317);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_1(numerical_op + 318);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_2(numerical_op + 319);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_3(numerical_op + 320);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_4(numerical_op + 321);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_5(numerical_op + 322);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_6(numerical_op + 323);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_7(numerical_op + 324);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_8(numerical_op + 325);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_1(numerical_op + 326);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_2(numerical_op + 327);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_3(numerical_op + 328);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_4(numerical_op + 329);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_5(numerical_op + 330);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_6(numerical_op + 331);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_7(numerical_op + 332);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_8(numerical_op + 333);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_9(numerical_op + 334);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_10(numerical_op + 335);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_11(numerical_op + 336);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_12(numerical_op + 337);
    OP_oneTr_p_0_0_0_Ir_4_C_1_n_1(numerical_op + 338);
    OP_oneTr_p_0_0_0_Ir_4_C_1_n_2(numerical_op + 339);
    OP_oneTr_p_0_0_0_Ir_4_C_1_n_3(numerical_op + 340);
    OP_oneTr_p_0_0_0_Ir_5_C_m1_n_1(numerical_op + 341);
    OP_oneTr_p_0_0_0_Ir_5_C_m1_n_2(numerical_op + 342);
    OP_oneTr_p_0_0_0_Ir_5_C_m1_n_3(numerical_op + 343);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_1(numerical_op + 344);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_2(numerical_op + 345);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_3(numerical_op + 346);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_4(numerical_op + 347);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_5(numerical_op + 348);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_6(numerical_op + 349);
    OP_oneTr_p_0_0_0_Ir_6_C_m1_n_1(numerical_op + 350);
    OP_oneTr_p_0_0_0_Ir_6_C_1_n_1(numerical_op + 351);
    OP_oneTr_p_0_0_0_Ir_7_C_m1_n_1(numerical_op + 352);
    OP_oneTr_p_0_0_0_Ir_7_C_1_n_1(numerical_op + 353);
    OP_oneTr_p_0_0_0_Ir_8_C_m1_n_1(numerical_op + 354);
    OP_oneTr_p_0_0_0_Ir_8_C_m1_n_2(numerical_op + 355);
    OP_oneTr_p_0_0_0_Ir_8_C_1_n_1(numerical_op + 356);
    OP_oneTr_p_0_0_0_Ir_8_C_1_n_2(numerical_op + 357);
    OP_oneTr_p_0_0_0_Ir_9_C_m1_n_1(numerical_op + 358);
    OP_oneTr_p_0_0_0_Ir_9_C_m1_n_2(numerical_op + 359);
    OP_oneTr_p_0_0_0_Ir_9_C_m1_n_3(numerical_op + 360);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_1(numerical_op + 361);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_2(numerical_op + 362);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_3(numerical_op + 363);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_4(numerical_op + 364);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_5(numerical_op + 365);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_6(numerical_op + 366);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_1(numerical_op + 367);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_2(numerical_op + 368);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_3(numerical_op + 369);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_4(numerical_op + 370);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_5(numerical_op + 371);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_6(numerical_op + 372);
    OP_oneTr_p_0_0_0_Ir_10_C_1_n_1(numerical_op + 373);
    OP_oneTr_p_0_0_0_Ir_10_C_1_n_2(numerical_op + 374);
    OP_oneTr_p_0_0_0_Ir_10_C_1_n_3(numerical_op + 375);
    eval_time_momentum_glueball_paths(t, 0, 0, 1);
    OP_oneTr_p_0_0_1_Ir_1_C_m1_n_1(numerical_op + 376);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_1(numerical_op + 377);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_2(numerical_op + 378);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_3(numerical_op + 379);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_4(numerical_op + 380);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_1(numerical_op + 381);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_2(numerical_op + 382);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_3(numerical_op + 383);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_4(numerical_op + 384);
    OP_oneTr_p_0_0_1_Ir_2_C_1_n_1(numerical_op + 385);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_1(numerical_op + 386);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_2(numerical_op + 387);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_3(numerical_op + 388);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_4(numerical_op + 389);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_1(numerical_op + 390);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_2(numerical_op + 391);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_3(numerical_op + 392);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_4(numerical_op + 393);
    OP_oneTr_p_0_0_1_Ir_4_C_m1_n_1(numerical_op + 394);
    OP_oneTr_p_0_0_1_Ir_4_C_m1_n_2(numerical_op + 395);
    OP_oneTr_p_0_0_1_Ir_4_C_1_n_1(numerical_op + 396);
    OP_oneTr_p_0_0_1_Ir_5_C_m1_n_1(numerical_op + 397);
    OP_oneTr_p_0_0_1_Ir_5_C_1_n_1(numerical_op + 398);
    OP_oneTr_p_0_0_1_Ir_5_C_1_n_2(numerical_op + 399);
    eval_time_momentum_glueball_paths(t, 0, 1, -1);
    OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_1(numerical_op + 400);
    OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_2(numerical_op + 401);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_1(numerical_op + 402);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_2(numerical_op + 403);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_3(numerical_op + 404);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_4(numerical_op + 405);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_1(numerical_op + 406);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_2(numerical_op + 407);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_3(numerical_op + 408);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_4(numerical_op + 409);
    OP_oneTr_p_0_1_m1_Ir_2_C_1_n_1(numerical_op + 410);
    OP_oneTr_p_0_1_m1_Ir_2_C_1_n_2(numerical_op + 411);
    OP_oneTr_p_0_1_m1_Ir_3_C_m1_n_1(numerical_op + 412);
    OP_oneTr_p_0_1_m1_Ir_3_C_m1_n_2(numerical_op + 413);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_1(numerical_op + 414);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_2(numerical_op + 415);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_3(numerical_op + 416);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_4(numerical_op + 417);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_1(numerical_op + 418);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_2(numerical_op + 419);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_3(numerical_op + 420);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_4(numerical_op + 421);
    OP_oneTr_p_0_1_m1_Ir_4_C_1_n_1(numerical_op + 422);
    OP_oneTr_p_0_1_m1_Ir_4_C_1_n_2(numerical_op + 423);
    eval_time_momentum_glueball_paths(t, 0, 1, 0);
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(numerical_op + 424);
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_2(numerical_op + 425);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_1(numerical_op + 426);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_2(numerical_op + 427);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_3(numerical_op + 428);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_4(numerical_op + 429);
    OP_oneTr_p_0_1_0_Ir_2_C_m1_n_1(numerical_op + 430);
    OP_oneTr_p_0_1_0_Ir_2_C_1_n_1(numerical_op + 431);
    OP_oneTr_p_0_1_0_Ir_3_C_m1_n_1(numerical_op + 432);
    OP_oneTr_p_0_1_0_Ir_3_C_m1_n_2(numerical_op + 433);
    OP_oneTr_p_0_1_0_Ir_3_C_1_n_1(numerical_op + 434);
    OP_oneTr_p_0_1_0_Ir_3_C_1_n_2(numerical_op + 435);
    OP_oneTr_p_0_1_0_Ir_4_C_m1_n_1(numerical_op + 436);
    OP_oneTr_p_0_1_0_Ir_4_C_m1_n_2(numerical_op + 437);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_1(numerical_op + 438);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_2(numerical_op + 439);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_3(numerical_op + 440);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_4(numerical_op + 441);
    OP_oneTr_p_0_1_0_Ir_5_C_m1_n_1(numerical_op + 442);
    OP_oneTr_p_0_1_0_Ir_5_C_1_n_1(numerical_op + 443);
    eval_time_momentum_glueball_paths(t, 0, 1, 1);
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_1(numerical_op + 444);
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_2(numerical_op + 445);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_1(numerical_op + 446);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_2(numerical_op + 447);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_3(numerical_op + 448);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_4(numerical_op + 449);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_1(numerical_op + 450);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_2(numerical_op + 451);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_3(numerical_op + 452);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_4(numerical_op + 453);
    OP_oneTr_p_0_1_1_Ir_2_C_1_n_1(numerical_op + 454);
    OP_oneTr_p_0_1_1_Ir_2_C_1_n_2(numerical_op + 455);
    OP_oneTr_p_0_1_1_Ir_3_C_m1_n_1(numerical_op + 456);
    OP_oneTr_p_0_1_1_Ir_3_C_m1_n_2(numerical_op + 457);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_1(numerical_op + 458);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_2(numerical_op + 459);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_3(numerical_op + 460);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_4(numerical_op + 461);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_1(numerical_op + 462);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_2(numerical_op + 463);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_3(numerical_op + 464);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_4(numerical_op + 465);
    OP_oneTr_p_0_1_1_Ir_4_C_1_n_1(numerical_op + 466);
    OP_oneTr_p_0_1_1_Ir_4_C_1_n_2(numerical_op + 467);
    eval_time_momentum_glueball_paths(t, 1, -1, -1);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_1(numerical_op + 468);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_2(numerical_op + 469);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_3(numerical_op + 470);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_4(numerical_op + 471);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_1(numerical_op + 472);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_2(numerical_op + 473);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_3(numerical_op + 474);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_4(numerical_op + 475);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_1(numerical_op + 476);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_2(numerical_op + 477);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_3(numerical_op + 478);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_4(numerical_op + 479);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_1(numerical_op + 480);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_2(numerical_op + 481);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_3(numerical_op + 482);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_4(numerical_op + 483);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_1(numerical_op + 484);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_2(numerical_op + 485);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_3(numerical_op + 486);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_4(numerical_op + 487);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_5(numerical_op + 488);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_6(numerical_op + 489);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_7(numerical_op + 490);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_8(numerical_op + 491);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_1(numerical_op + 492);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_2(numerical_op + 493);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_3(numerical_op + 494);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_4(numerical_op + 495);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_5(numerical_op + 496);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_6(numerical_op + 497);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_7(numerical_op + 498);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_8(numerical_op + 499);
    eval_time_momentum_glueball_paths(t, 1, -1, 0);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_1(numerical_op + 500);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_2(numerical_op + 501);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_3(numerical_op + 502);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_4(numerical_op + 503);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_1(numerical_op + 504);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_2(numerical_op + 505);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_3(numerical_op + 506);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_4(numerical_op + 507);
    OP_oneTr_p_1_m1_0_Ir_2_C_m1_n_1(numerical_op + 508);
    OP_oneTr_p_1_m1_0_Ir_2_C_1_n_1(numerical_op + 509);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_1(numerical_op + 510);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_2(numerical_op + 511);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_3(numerical_op + 512);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_4(numerical_op + 513);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_1(numerical_op + 514);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_2(numerical_op + 515);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_3(numerical_op + 516);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_4(numerical_op + 517);
    OP_oneTr_p_1_m1_0_Ir_4_C_m1_n_1(numerical_op + 518);
    OP_oneTr_p_1_m1_0_Ir_4_C_1_n_1(numerical_op + 519);
    eval_time_momentum_glueball_paths(t, 1, -1, 1);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_1(numerical_op + 520);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_2(numerical_op + 521);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_3(numerical_op + 522);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_4(numerical_op + 523);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_1(numerical_op + 524);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_2(numerical_op + 525);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_3(numerical_op + 526);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_4(numerical_op + 527);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_1(numerical_op + 528);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_2(numerical_op + 529);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_3(numerical_op + 530);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_4(numerical_op + 531);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_1(numerical_op + 532);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_2(numerical_op + 533);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_3(numerical_op + 534);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_4(numerical_op + 535);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_1(numerical_op + 536);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_2(numerical_op + 537);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_3(numerical_op + 538);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_4(numerical_op + 539);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_5(numerical_op + 540);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_6(numerical_op + 541);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_7(numerical_op + 542);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_8(numerical_op + 543);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_1(numerical_op + 544);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_2(numerical_op + 545);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_3(numerical_op + 546);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_4(numerical_op + 547);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_5(numerical_op + 548);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_6(numerical_op + 549);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_7(numerical_op + 550);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_8(numerical_op + 551);
    eval_time_momentum_glueball_paths(t, 1, 0, -1);
    OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_1(numerical_op + 552);
    OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_2(numerical_op + 553);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_1(numerical_op + 554);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_2(numerical_op + 555);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_3(numerical_op + 556);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_4(numerical_op + 557);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_1(numerical_op + 558);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_2(numerical_op + 559);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_3(numerical_op + 560);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_4(numerical_op + 561);
    OP_oneTr_p_1_0_m1_Ir_2_C_1_n_1(numerical_op + 562);
    OP_oneTr_p_1_0_m1_Ir_2_C_1_n_2(numerical_op + 563);
    OP_oneTr_p_1_0_m1_Ir_3_C_m1_n_1(numerical_op + 564);
    OP_oneTr_p_1_0_m1_Ir_3_C_m1_n_2(numerical_op + 565);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_1(numerical_op + 566);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_2(numerical_op + 567);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_3(numerical_op + 568);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_4(numerical_op + 569);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_1(numerical_op + 570);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_2(numerical_op + 571);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_3(numerical_op + 572);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_4(numerical_op + 573);
    OP_oneTr_p_1_0_m1_Ir_4_C_1_n_1(numerical_op + 574);
    OP_oneTr_p_1_0_m1_Ir_4_C_1_n_2(numerical_op + 575);
    eval_time_momentum_glueball_paths(t, 1, 0, 0);
    OP_oneTr_p_1_0_0_Ir_1_C_m1_n_1(numerical_op + 576);
    OP_oneTr_p_1_0_0_Ir_1_C_m1_n_2(numerical_op + 577);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_1(numerical_op + 578);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_2(numerical_op + 579);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_3(numerical_op + 580);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_4(numerical_op + 581);
    OP_oneTr_p_1_0_0_Ir_2_C_m1_n_1(numerical_op + 582);
    OP_oneTr_p_1_0_0_Ir_2_C_1_n_1(numerical_op + 583);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_1(numerical_op + 584);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_2(numerical_op + 585);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_3(numerical_op + 586);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_4(numerical_op + 587);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_5(numerical_op + 588);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_6(numerical_op + 589);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_7(numerical_op + 590);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_8(numerical_op + 591);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_1(numerical_op + 592);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_2(numerical_op + 593);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_3(numerical_op + 594);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_4(numerical_op + 595);
    OP_oneTr_p_1_0_0_Ir_4_C_m1_n_1(numerical_op + 596);
    OP_oneTr_p_1_0_0_Ir_4_C_m1_n_2(numerical_op + 597);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_1(numerical_op + 598);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_2(numerical_op + 599);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_3(numerical_op + 600);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_4(numerical_op + 601);
    OP_oneTr_p_1_0_0_Ir_5_C_m1_n_1(numerical_op + 602);
    OP_oneTr_p_1_0_0_Ir_5_C_1_n_1(numerical_op + 603);
    eval_time_momentum_glueball_paths(t, 1, 0, 1);
    OP_oneTr_p_1_0_1_Ir_1_C_m1_n_1(numerical_op + 604);
    OP_oneTr_p_1_0_1_Ir_1_C_m1_n_2(numerical_op + 605);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_1(numerical_op + 606);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_2(numerical_op + 607);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_3(numerical_op + 608);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_4(numerical_op + 609);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_1(numerical_op + 610);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_2(numerical_op + 611);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_3(numerical_op + 612);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_4(numerical_op + 613);
    OP_oneTr_p_1_0_1_Ir_2_C_1_n_1(numerical_op + 614);
    OP_oneTr_p_1_0_1_Ir_2_C_1_n_2(numerical_op + 615);
    OP_oneTr_p_1_0_1_Ir_3_C_m1_n_1(numerical_op + 616);
    OP_oneTr_p_1_0_1_Ir_3_C_m1_n_2(numerical_op + 617);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_1(numerical_op + 618);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_2(numerical_op + 619);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_3(numerical_op + 620);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_4(numerical_op + 621);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_1(numerical_op + 622);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_2(numerical_op + 623);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_3(numerical_op + 624);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_4(numerical_op + 625);
    OP_oneTr_p_1_0_1_Ir_4_C_1_n_1(numerical_op + 626);
    OP_oneTr_p_1_0_1_Ir_4_C_1_n_2(numerical_op + 627);
    eval_time_momentum_glueball_paths(t, 1, 1, -1);
    OP_oneTr_p_1_1_m1_Ir_1_C_m1_n_1(numerical_op + 628);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_1(numerical_op + 629);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_2(numerical_op + 630);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_3(numerical_op + 631);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_4(numerical_op + 632);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_1(numerical_op + 633);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_2(numerical_op + 634);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_3(numerical_op + 635);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_4(numerical_op + 636);
    OP_oneTr_p_1_1_m1_Ir_2_C_1_n_1(numerical_op + 637);
    OP_oneTr_p_1_1_m1_Ir_3_C_m1_n_1(numerical_op + 638);
    OP_oneTr_p_1_1_m1_Ir_3_C_m1_n_2(numerical_op + 639);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_1(numerical_op + 640);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_2(numerical_op + 641);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_3(numerical_op + 642);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_4(numerical_op + 643);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_5(numerical_op + 644);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_6(numerical_op + 645);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_7(numerical_op + 646);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_8(numerical_op + 647);
    eval_time_momentum_glueball_paths(t, 1, 1, 0);
    OP_oneTr_p_1_1_0_Ir_1_C_m1_n_1(numerical_op + 648);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_1(numerical_op + 649);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_2(numerical_op + 650);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_3(numerical_op + 651);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_4(numerical_op + 652);
    OP_oneTr_p_1_1_0_Ir_2_C_m1_n_1(numerical_op + 653);
    OP_oneTr_p_1_1_0_Ir_2_C_1_n_1(numerical_op + 654);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_1(numerical_op + 655);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_2(numerical_op + 656);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_3(numerical_op + 657);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_4(numerical_op + 658);
    OP_oneTr_p_1_1_0_Ir_3_C_1_n_1(numerical_op + 659);
    OP_oneTr_p_1_1_0_Ir_4_C_m1_n_1(numerical_op + 660);
    OP_oneTr_p_1_1_0_Ir_4_C_1_n_1(numerical_op + 661);
    eval_time_momentum_glueball_paths(t, 1, 1, 1);
    OP_oneTr_p_1_1_1_Ir_1_C_m1_n_1(numerical_op + 662);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_1(numerical_op + 663);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_2(numerical_op + 664);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_3(numerical_op + 665);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_4(numerical_op + 666);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_1(numerical_op + 667);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_2(numerical_op + 668);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_3(numerical_op + 669);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_4(numerical_op + 670);
    OP_oneTr_p_1_1_1_Ir_2_C_1_n_1(numerical_op + 671);
    OP_oneTr_p_1_1_1_Ir_3_C_m1_n_1(numerical_op + 672);
    OP_oneTr_p_1_1_1_Ir_3_C_m1_n_2(numerical_op + 673);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_1(numerical_op + 674);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_2(numerical_op + 675);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_3(numerical_op + 676);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_4(numerical_op + 677);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_5(numerical_op + 678);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_6(numerical_op + 679);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_7(numerical_op + 680);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_8(numerical_op + 681);
}

void eval_all_glueball_ops_p_m1_m1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_m1_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_m1_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_m1_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_m1_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_m1_m1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_0_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_3_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_4_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_0_Ir_4_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_m1_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_m1_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_1_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_m1_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_m1_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_m1_1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_m1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_m1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_m1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_m1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_m1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_0_0_Ir_3_C_m1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_0_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_4_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_0_Ir_4_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_5_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_0_0_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_0_Ir_5_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_0_1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_0_1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_0_1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_0_1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_m1_1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_m1_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_m1_Ir_2_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_1_m1_Ir_3_C_m1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_m1_1_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_1_m1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_0_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_0_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_0_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_4_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_1_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_0_Ir_4_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_m1_1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_1_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_1_Ir_2_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_m1_1_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_1_1_Ir_3_C_m1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_m1_1_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_m1_1_1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_m1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_m1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_m1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_m1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_m1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_0_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_0_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_0_Ir_3_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_0_Ir_4_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_0_Ir_4_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_5_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_m1_0_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_0_Ir_5_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_m1_1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_m1_1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_m1_1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_m1_1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_m1_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_m1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_m1_Ir_4_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_4_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_5_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_m1_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_m1_Ir_5_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_m1_Ir_5_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_2_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_0_0_0_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_0_0_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_8(numerical_op + 7);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_9(numerical_op + 8);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_10(numerical_op + 9);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_11(numerical_op + 10);
    OP_oneTr_p_0_0_0_Ir_4_C_m1_n_12(numerical_op + 11);
}

void eval_all_glueball_ops_p_0_0_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_4_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_4_C_1_n_3(numerical_op + 2);
}

void eval_all_glueball_ops_p_0_0_0_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_5_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_5_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_5_C_m1_n_3(numerical_op + 2);
}

void eval_all_glueball_ops_p_0_0_0_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_0_0_0_Ir_5_C_1_n_6(numerical_op + 5);
}

void eval_all_glueball_ops_p_0_0_0_Ir_6_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_6_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_0_Ir_6_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_6_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_0_Ir_7_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_7_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_0_Ir_7_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_7_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_0_Ir_8_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_8_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_8_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_0_Ir_8_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_8_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_8_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_0_Ir_9_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_9_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_9_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_9_C_m1_n_3(numerical_op + 2);
}

void eval_all_glueball_ops_p_0_0_0_Ir_9_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_0_0_0_Ir_9_C_1_n_6(numerical_op + 5);
}

void eval_all_glueball_ops_p_0_0_0_Ir_10_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_0_0_0_Ir_10_C_m1_n_6(numerical_op + 5);
}

void eval_all_glueball_ops_p_0_0_0_Ir_10_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_0_Ir_10_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_0_Ir_10_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_0_Ir_10_C_1_n_3(numerical_op + 2);
}

void eval_all_glueball_ops_p_0_0_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_1_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_0_1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_0_1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_4_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_0_1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_4_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_1_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_5_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_0_1_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_0_1_Ir_5_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_0_1_Ir_5_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_m1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_m1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_m1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_m1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_m1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_1_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_1_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_3_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_4_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_0_Ir_4_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_0_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_5_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_1_0_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_0_Ir_5_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_0_1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_0_1_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_0_1_1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_0_1_1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_0_1_1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_0_1_1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_m1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_m1_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_m1_Ir_2_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_1_m1_m1_Ir_3_C_m1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_1_m1_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_1_m1_m1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_0_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_0_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_0_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_4_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_m1_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_0_Ir_4_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_m1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_1_Ir_1_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_1_Ir_2_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_m1_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_1_m1_1_Ir_3_C_m1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_1_m1_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_1_m1_1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_m1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_m1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_m1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_m1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_m1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_0_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_0_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_4(numerical_op + 3);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_5(numerical_op + 4);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_6(numerical_op + 5);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_7(numerical_op + 6);
    OP_oneTr_p_1_0_0_Ir_3_C_m1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_1_0_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_0_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_4_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_0_Ir_4_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_0_Ir_5_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_5_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_0_0_Ir_5_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_0_Ir_5_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_0_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_1_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_1_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_2_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_2_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_0_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_1_Ir_3_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_1_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_0_1_Ir_4_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_0_1_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_0_1_Ir_4_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_0_1_Ir_4_C_1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_1_m1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_m1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_m1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_m1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_1_m1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_m1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_1_m1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_m1_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_m1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_m1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_m1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_1_m1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_1_1_m1_Ir_3_C_1_n_8(numerical_op + 7);
}

void eval_all_glueball_ops_p_1_1_0_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_0_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_0_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_1_0_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_2_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_0_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_0_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_0_Ir_3_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_1_0_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_3_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_0_Ir_4_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_4_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_0_Ir_4_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_0_Ir_4_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_1_Ir_1_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_1_Ir_1_C_m1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_1_Ir_1_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_1_Ir_1_C_1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_1_1_Ir_2_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_1_Ir_2_C_m1_n_4(numerical_op + 3);
}

void eval_all_glueball_ops_p_1_1_1_Ir_2_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_1_Ir_2_C_1_n_1(numerical_op + 0);
}

void eval_all_glueball_ops_p_1_1_1_Ir_3_C_m1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_1_Ir_3_C_m1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_1_Ir_3_C_m1_n_2(numerical_op + 1);
}

void eval_all_glueball_ops_p_1_1_1_Ir_3_C_1(double complex *numerical_op)
{
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_1(numerical_op + 0);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_2(numerical_op + 1);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_3(numerical_op + 2);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_4(numerical_op + 3);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_5(numerical_op + 4);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_6(numerical_op + 5);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_7(numerical_op + 6);
    OP_oneTr_p_1_1_1_Ir_3_C_1_n_8(numerical_op + 7);
}
