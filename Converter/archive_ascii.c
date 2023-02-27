/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "libhr.h"

void read_gauge_field_ascii(char filename[]) {
    int counter = 0, counter0 = 0, pointcounter = 0;
    int Vdone[GLB_T][GLB_X][GLB_Y][GLB_Z][4];

    Timer clock;
    timer_set(&clock);

    FILE *fp = fopen(filename, "r");
    error(fp == NULL, 1, "read_gauge_field_ascii", "Failed to open file for reading");

    for (int g0 = 0; g0 < GLB_T; g0++) {
        for (int g1 = 0; g1 < GLB_X; g1++) {
            for (int g2 = 0; g2 < GLB_Y; g2++) {
                for (int g3 = 0; g3 < GLB_Z; g3++) {
                    for (int mu = 0; mu < 4; mu++) {
                        Vdone[g0][g1][g2][g3][mu] = 0;
                    }
                }
            }
        }
    }

    while (1) {
        /* u( row , col , x , y , z , t , dir) = (re,im) */
        int g0, g1, g2, g3;
        int hm = fscanf(fp, " %d %d %d %d\n", &g0, &g1, &g2, &g3);
        if (hm != 4) { break; }
        pointcounter++;
        for (int mu = 0; mu < 4; mu++) {
            suNg tmpmat;

            for (int gamma = 0; gamma < NG; gamma++) {
                for (int alpha = 0; alpha < NG; alpha++) {
                    double re, im;
                    hm = fscanf(fp, " %lf %lf\n", &re, &im);
                    if (hm != 2) { error(0, 1, "read_gauge_field_ascii", "Bad number of element in the gauge field\n"); }
                    tmpmat.c[gamma * NG + alpha] = re + I * im;
                    counter++;
                }
            }
            *pu_gauge(ipt(g0 + 1, g1 - 1, g2 - 1, g3 - 1), mu) = tmpmat;
            Vdone[g0 + 1][g1 - 1][g2 - 1][g3 - 1][mu] = 1;
        }
    }

    for (int g0 = 0; g0 < GLB_T; g0++) {
        for (int g1 = 0; g1 < GLB_X; g1++) {
            for (int g2 = 0; g2 < GLB_Y; g2++) {
                for (int g3 = 0; g3 < GLB_Z; g3++) {
                    for (int mu = 0; mu < 4; mu++) {
                        if (Vdone[g0][g1][g2][g3][mu] == 0) {
                            counter0 += NG * NG;
                            _suNg_zero(*pu_gauge(ipt(g0, g1, g2, g3), mu));
                        }
                    }
                }
            }
        }
    }

    lprintf("IO", 0, "Read %d lines\n", counter + pointcounter);
    error(counter + counter0 != NG * NG * 4 * GLB_T * GLB_X * GLB_Y * GLB_Z, 1, "read_gauge_field_ascii " __FILE__,
          "Bad number of lines in file");

    fclose(fp);

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Configuration [%s] read [%lf sec]\n", filename, elapsed_sec);
}
