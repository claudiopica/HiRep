/*******************************************************************************
*
* File ranlxs.c
*
* Random number generator "ranlxs"
*
*******************************************************************************/

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"

#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static float one_bit;
static vec_t carry;

static union
{
   dble_vec_t vec[12];
   int num[96];
} x;

#define STEP(pi,pj) \
      d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
      (*pi).c2.c1+=(d<0); \
      d+=BASE; \
      (*pi).c1.c1=d&MASK; \
      d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
      (*pi).c2.c2+=(d<0); \
      d+=BASE; \
      (*pi).c1.c2=d&MASK; \
      d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
      (*pi).c2.c3+=(d<0); \
      d+=BASE; \
      (*pi).c1.c3=d&MASK; \
      d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
      (*pi).c2.c4+=(d<0); \
      d+=BASE; \
      (*pi).c1.c4=d&MASK; \
      d=(*pj).c2.c1-(*pi).c2.c1; \
      carry.c1=(d<0); \
      d+=BASE; \
      (*pi).c2.c1=d&MASK; \
      d=(*pj).c2.c2-(*pi).c2.c2; \
      carry.c2=(d<0); \
      d+=BASE; \
      (*pi).c2.c2=d&MASK; \
      d=(*pj).c2.c3-(*pi).c2.c3; \
      carry.c3=(d<0); \
      d+=BASE; \
      (*pi).c2.c3=d&MASK; \
      d=(*pj).c2.c4-(*pi).c2.c4; \
      carry.c4=(d<0); \
      d+=BASE; \
      (*pi).c2.c4=d&MASK


static void update()
{
   int k,kmax,d;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];
      
   for (k=0;k<kmax;k++) 
   {
      STEP(pi,pj);
      pi+=1;
      pj+=1;
      if (pi==pmax)
         pi=pmin;      
      if (pj==pmax)
         pj=pmin; 
   }

   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


static void define_constants()
{
   int k;

   one_bit=(float)(ldexp(1.0,-24));

   for (k=0;k<96;k++)
      next[k]=(k+1)%96;
}


void rlxs_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   error((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24),1,
         "rlxs_init [ranlxs.c]",
         "Arithmetic on this machine is not suitable for ranlxs");

   define_constants();
   
   error((level<0)||(level>2),1,"rlxs_init [ranlxs.c]",
         "Bad choice of luxury level (should be 0,1 or 2)");   
   
   if (level==0)
      pr=109;
   else if (level==1)
      pr=202;
   else if (level==2)
      pr=397;
   
   i=seed;

   for (k=0;k<31;k++) 
   {
      xbit[k]=i%2;
      i/=2;
   }

   error((seed<=0)||(i!=0),1,"rlxs_init [ranlxs.c]",
         "Bad choice of seed (should be between 1 and 2^31-1)");

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++) 
         {
            iy=xbit[ibit];
            ix=2*ix+iy;
         
            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)==i)
            ix=16777215-ix;

         x.num[4*k+i]=ix;
      }
   }

   carry.c1=0;
   carry.c2=0;
   carry.c3=0;
   carry.c4=0;

   ir=0;
   jr=7;
   is=95;
   is_old=0;
   prm=pr%12;
   init=1;
}


void ranlxs(float r[],int n)
{
   int k;

   if (init==0)
      rlxs_init(0,1);

   for (k=0;k<n;k++) 
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=one_bit*(float)(x.num[is]);      
   }
}


int rlxs_size(void)
{
   return(105);
}


void rlxs_get(int state[])
{
   int k;

   error(init==0,1,"rlxs_get [ranlxs.c]",
         "Undefined state (ranlxs is not initialized");

   state[0]=rlxs_size();

   for (k=0;k<96;k++)
      state[k+1]=x.num[k];

   state[97]=carry.c1;
   state[98]=carry.c2;
   state[99]=carry.c3;
   state[100]=carry.c4;

   state[101]=pr;
   state[102]=ir;
   state[103]=jr;
   state[104]=is;
}


void rlxs_reset(int state[])
{
   int k;

   error((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24),1,
         "rlxs_reset [ranlxs.c]",
         "Arithmetic on this machine is not suitable for ranlxs");

   define_constants();

   error(state[0]!=rlxs_size(),1,"rlxs_reset [ranlxs.c]",
         "Unexpected input data");   

   for (k=0;k<96;k++)
   {
      error((state[k+1]<0)||(state[k+1]>=167777216),1,
            "rlxs_reset [ranlxs.c]","Unexpected input data");      

      x.num[k]=state[k+1];
   }

   error(((state[97]!=0)&&(state[97]!=1))||
         ((state[98]!=0)&&(state[98]!=1))||
         ((state[99]!=0)&&(state[99]!=1))||
         ((state[100]!=0)&&(state[100]!=1)),1,
         "rlxs_reset [ranlxs.c]","Unexpected input data");   
   
   carry.c1=state[97];
   carry.c2=state[98];
   carry.c3=state[99];
   carry.c4=state[100];

   pr=state[101];
   ir=state[102];
   jr=state[103];
   is=state[104];
   is_old=8*ir;
   prm=pr%12;
   init=1;

   error(((pr!=109)&&(pr!=202)&&(pr!=397))||
         (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
         (is<0)||(is>95),1,
         "rlxs_reset [ranlxs.c]","Unexpected input data");   
}


void rlxs_read_random(char filename[])
{
  int i,nran,*state;
  FILE *fp;

  nran=rlxs_size();
  state=malloc(nran*sizeof(int));

  error((fp=fopen(filename,"r"))==NULL,1,"read_random [main_utils.c]",
        "Failed to open file for reading state of ranlux\n");

  for (i=0;i<nran;i++)
    fscanf(fp,"%d",&state[i]);

  error(ferror(fp)!=0,1,"read_random [main_utils.c]",
        "Failed to read state of ranlux from file\n");

  fclose(fp);

  rlxs_reset(state);
  free(state);
}

void rlxs_write_random(char filename[])
{
  int i,nran,*state;
  FILE *fp;

  nran=rlxs_size();
  state=malloc(nran*sizeof(int));

  rlxs_get(state);

  error((fp=fopen(filename,"w"))==NULL,1,"write_random [main_utils.c]",
        "Failed to open file for writing state of ranlux\n");

  for (i=0;i<nran;i++)
    fprintf(fp,"%d\n",state[i]);

  error(ferror(fp)!=0,1,"write_random [main_utils.c]",
        "Failed to write state of ranlux to file\n");

  fclose(fp);
  free(state);
}


