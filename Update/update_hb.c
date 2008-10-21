/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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
#include "communications.h"

static suNg v __attribute__ ((aligned (16)));

void project_gauge_field(void)
{
  _DECLARE_INT_ITERATOR(ix);

  _MASTER_FOR(&glattice,ix) {
    project_to_suNg(pu_gauge(ix,0));
    project_to_suNg(pu_gauge(ix,1));
    project_to_suNg(pu_gauge(ix,2));
    project_to_suNg(pu_gauge(ix,3));
  }

  start_gf_sendrecv(u_gauge);

}


static void update_all(double beta,int type)
{
  _DECLARE_INT_ITERATOR(ix);
   static int count=PROJECT_INTERVAL;

   if (count>=PROJECT_INTERVAL) {
     project_gauge_field();
     count=0;
   }
   ++count;

  _MASTER_FOR(&glattice,ix) {
    staples(ix,0,&v);
    cabmar(beta,pu_gauge(ix,0),&v,type);
    staples(ix,1,&v);
    cabmar(beta,pu_gauge(ix,1),&v,type);
    staples(ix,2,&v);
    cabmar(beta,pu_gauge(ix,2),&v,type);
    staples(ix,3,&v);
    cabmar(beta,pu_gauge(ix,3),&v,type);
   }
} 


void update(double beta,int nhb,int nor)
{
   int n;

   for (n=0;n<nhb;n++)
      update_all(beta,0);

   for (n=0;n<nor;n++)
      update_all(beta,1);
}

