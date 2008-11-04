/* This is part of the ranlux 3.3 distribution by Martin Luesher */

/*******************************************************************************
*
* File testlx.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* This program checks that ranlxs and ranlxd work correctly
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Include/ranlux.h"

#define NXS 204
#define NXD 99

int main(void)
{
   int k,test1,test2;
   int *state1,*state2;
   float sbase;
   float xs[NXS],ys[NXS],xsn[96];
   double base;
   double xd[NXD],yd[NXD],xdn[48];

   sbase=(float)(ldexp(1.0,24));
   base=ldexp(1.0,48);
   state1=malloc(rlxs_size()*sizeof(int));
   state2=malloc(rlxd_size()*sizeof(int));
   
   rlxs_init(0,32767);
   rlxd_init(1,32767);


/*******************************************************************************
*
* Check that the correct sequences of random numbers are obtained
*
*******************************************************************************/

   for (k=0;k<20;k++)
   {
      ranlxs(xs,NXS);
      ranlxd(xd,NXD);
   }

   xsn[0]=13257445.0f;
   xsn[1]=15738482.0f;
   xsn[2]=5448599.0f;
   xsn[3]=9610459.0f;
   xsn[4]=1046025.0f;
   xsn[5]=2811360.0f;
   xsn[6]=14923726.0f;
   xsn[7]=2287739.0f;
   xsn[8]=16133204.0f;
   xsn[9]=16328320.0f;
   xsn[10]=12980218.0f;
   xsn[11]=9256959.0f;
   xsn[12]=5633754.0f;
   xsn[13]=7422961.0f;
   xsn[14]=6032411.0f;
   xsn[15]=14970828.0f;
   xsn[16]=10717272.0f;
   xsn[17]=2520878.0f;
   xsn[18]=8906135.0f;
   xsn[19]=8507426.0f;
   xsn[20]=11925022.0f;
   xsn[21]=12042827.0f;
   xsn[22]=12263021.0f;
   xsn[23]=4828801.0f;
   xsn[24]=5300508.0f;
   xsn[25]=13346776.0f;
   xsn[26]=10869790.0f;
   xsn[27]=8520207.0f;
   xsn[28]=11213953.0f;
   xsn[29]=14439320.0f;
   xsn[30]=5716476.0f;
   xsn[31]=13600448.0f;
   xsn[32]=12545579.0f;
   xsn[33]=3466523.0f;
   xsn[34]=113906.0f;
   xsn[35]=10407879.0f;
   xsn[36]=12058596.0f;
   xsn[37]=4390921.0f;
   xsn[38]=1634350.0f;
   xsn[39]=9823280.0f;
   xsn[40]=12569690.0f;
   xsn[41]=8267856.0f;
   xsn[42]=5869501.0f;
   xsn[43]=7210219.0f;
   xsn[44]=1362361.0f;
   xsn[45]=2956909.0f;
   xsn[46]=504465.0f;
   xsn[47]=6664636.0f;
   xsn[48]=6048963.0f;
   xsn[49]=1098525.0f;
   xsn[50]=1261330.0f;
   xsn[51]=2401071.0f;
   xsn[52]=8087317.0f;
   xsn[53]=1293933.0f;
   xsn[54]=555494.0f;
   xsn[55]=14872475.0f;
   xsn[56]=11261534.0f;
   xsn[57]=166813.0f;
   xsn[58]=13424516.0f;
   xsn[59]=15280818.0f;
   xsn[60]=4644497.0f;
   xsn[61]=6333595.0f;
   xsn[62]=10012569.0f;
   xsn[63]=6878028.0f;
   xsn[64]=9176136.0f;
   xsn[65]=8379433.0f;
   xsn[66]=11073957.0f;
   xsn[67]=2465529.0f;
   xsn[68]=13633550.0f;
   xsn[69]=12721649.0f;
   xsn[70]=569725.0f;
   xsn[71]=6375015.0f;
   xsn[72]=2164250.0f;
   xsn[73]=6725885.0f;
   xsn[74]=7223108.0f;
   xsn[75]=4890858.0f;
   xsn[76]=11298261.0f;
   xsn[77]=12086020.0f;
   xsn[78]=4447706.0f;
   xsn[79]=1164782.0f;
   xsn[80]=1904399.0f;
   xsn[81]=16669839.0f;
   xsn[82]=2586766.0f;
   xsn[83]=3605708.0f;
   xsn[84]=15761082.0f;
   xsn[85]=14937769.0f;
   xsn[86]=13965017.0f;
   xsn[87]=2175021.0f;
   xsn[88]=16668997.0f;
   xsn[89]=13996602.0f;
   xsn[90]=6313099.0f;
   xsn[91]=15646036.0f;
   xsn[92]=9746447.0f;
   xsn[93]=9596781.0f;
   xsn[94]=9244169.0f;
   xsn[95]=4731726.0f;

   xdn[0]=135665102723086.0;
   xdn[1]=259840970195871.0;
   xdn[2]=110726726657103.0;
   xdn[3]=53972500363809.0;
   xdn[4]=199301297412157.0;
   xdn[5]=63744794353870.0;
   xdn[6]=178745978725904.0;
   xdn[7]=243549380863176.0;
   xdn[8]=244796821836177.0;
   xdn[9]=223788809121855.0;
   xdn[10]=113720856430443.0;
   xdn[11]=124607822268499.0;
   xdn[12]=25705458431399.0;
   xdn[13]=155476863764950.0;
   xdn[14]=195602097736933.0;
   xdn[15]=183038707238950.0;
   xdn[16]=62268883953527.0;
   xdn[17]=157047615112119.0;
   xdn[18]=58134973897037.0;
   xdn[19]=26908869337679.0;
   xdn[20]=259927185454290.0;
   xdn[21]=130534606773507.0;
   xdn[22]=205295065526788.0;
   xdn[23]=40201323262686.0;
   xdn[24]=193822255723177.0;
   xdn[25]=239720285097881.0;
   xdn[26]=54433631586673.0;
   xdn[27]=31313178820772.0;
   xdn[28]=152904879618865.0;
   xdn[29]=256187025780734.0;
   xdn[30]=110292144635528.0;
   xdn[31]=26555117184469.0;
   xdn[32]=228913371644996.0;
   xdn[33]=126837665590799.0;
   xdn[34]=141069100232139.0;
   xdn[35]=96171028602910.0;
   xdn[36]=259271018918511.0;
   xdn[37]=65257892816619.0;
   xdn[38]=14254344610711.0;
   xdn[39]=137794868158301.0;
   xdn[40]=269703238916504.0;
   xdn[41]=35782602710520.0;
   xdn[42]=51447305327263.0;
   xdn[43]=247852246697199.0;
   xdn[44]=65072958134912.0;
   xdn[45]=273325640150591.0;
   xdn[46]=2768714666444.0;
   xdn[47]=173907458721736.0;
   
   test1=0;
   test2=0;

   for (k=0;k<96;k++)
   {
      if (xsn[k]!=(xs[k+60]*sbase))
         test1=1;
   }

   for (k=0;k<48;k++)
   {
      if (xdn[k]!=(xd[k+39]*base))
         test2=1;
   }

   if (test1==1)
   {
      printf("\n");
      printf("Test failed: ranlxs gives incorrect results\n");
      printf("=> do not use ranlxs on this machine\n");
      printf("\n");
   }

   if (test2==1)
   {
      printf("\n");
      printf("Test failed: ranlxd gives incorrect results\n");
      printf("=> do not use ranlxd on this machine\n");
      printf("\n");
   }


/*******************************************************************************
*
* Check of the I/O routines
*
*******************************************************************************/

   rlxs_get(state1);
   rlxd_get(state2);

   for (k=0;k<10;k++)
   {
      ranlxs(xs,NXS);
      ranlxd(xd,NXD);
   }

   rlxs_reset(state1);
   rlxd_reset(state2);

   for (k=0;k<10;k++)
   {
      ranlxs(ys,NXS);
      ranlxd(yd,NXD);
   }

   for (k=0;k<NXS;k++)
   {
      if (xs[k]!=ys[k])
         test1=2;
   }

   for (k=0;k<NXD;k++)
   {
      if (xd[k]!=yd[k])
         test2=2;
   }

   if (test1==2)
   {
      printf("\n");
      printf("Test failed: I/O routines for ranlxs do not work properly\n");
      printf("=> do not use ranlxs on this machine\n");
      printf("\n");
   }

   if (test2==2)
   {
      printf("\n");
      printf("Test failed: I/O routines for ranlxd do not work properly\n");
      printf("=> do not use ranlxd on this machine\n");
      printf("\n");
   }


/*******************************************************************************
*
* Success messages
*
*******************************************************************************/

   if ((test1==0)&&(test2==0))
   {
      printf("\n");
      printf("All tests passed\n");
      printf("=> ranlxs and ranlxd work correctly on this machine\n");
      printf("\n");
   }
   else if (test1==0)
   {
      printf("\n");
      printf("All tests on ranlxs passed\n");
      printf("=> ranlxs works correctly on this machine\n");
      printf("\n");
   }
   else if (test2==0)
   {
      printf("\n");
      printf("All tests on ranlxd passed\n");
      printf("=> ranlxd works correctly on this machine\n");
      printf("\n");
   }
   exit(0);
}

