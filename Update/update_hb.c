/*******************************************************************************
*
* File update.c
*
* Update programs
*
*******************************************************************************/

#define PROJECT_INTERVAL 10

#include "suN.h"
#include "utils.h"
#include "global.h"
#include "update.h"

static suNg v __attribute__ ((aligned (16)));

void project_gauge_field(void)
{
   int ix,mu;

   for (ix=0;ix<VOLUME;ix++)
     for (mu=0;mu<4;mu++)
       project_to_suNg(pu_gauge(ix,mu));
}


static void update_all(float beta,int type)
{
   int ix,mu;
   static int count=PROJECT_INTERVAL;

   if (count>=PROJECT_INTERVAL) {
     project_gauge_field();
     count=0;
   }
   ++count;

   for (ix=0;ix<VOLUME;ix++) {
     for (mu=0;mu<4;mu++) {
       staples(ix,mu,&v);
       cabmar(beta,pu_gauge(ix,mu),&v,type);
     }
   }
} 


void update(float beta,int nhb,int nor)
{
   int n;

   for (n=0;n<nhb;n++)
      update_all(beta,0);

   for (n=0;n<nor;n++)
      update_all(beta,1);
}

