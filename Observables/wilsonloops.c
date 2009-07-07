#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include <math.h>


#ifdef WITH_MPI
#error Please compile without MPI!
#endif


int multistep(int ix,int mu,int dx) {
  if(dx==0) return ix;
  else if(dx>0) return multistep(iup(ix,mu),mu,dx-1);
  else if(dx<0) return multistep(idn(ix,mu),mu,dx+1);
  return 0;
}

void wilsonloops(int mu, int x, suNg_field* g) {
  suNg_field* gfixed = alloc_gfield(&glattice);
  int i, y, nu;
  _DECLARE_INT_ITERATOR(ix);
  int ix1;
  suNg tmp, tmp1;
  int glb_size[4] = {GLB_T, GLB_X, GLB_Y, GLB_Z};
  complex wilson_fund;
  complex ctmp;
  double wilson_adj;

  error(g->type!=&glattice,1,"wilsonloops.c","The gauge field must be defined on the whole lattice");
  error(x<=0,1,"wilsonloops.c","The length must be > 0");
  error(mu<0 || mu>3,1,"wilsonloops.c","Not defined direction");

  suNg_field_copy(gfixed,g);

  
  _MASTER_FOR(&glattice,ix) {
    ix1=ix;
    for(i=1;i<x;i++) {
      ix1=iup(ix1,mu);
      _suNg_times_suNg(tmp,*_4FIELD_AT(gfixed,ix,mu),*_4FIELD_AT(g,ix1,mu));
      *_4FIELD_AT(gfixed,ix,mu)=tmp;
    }
  }
  
  for(nu=0;nu<4;nu++) {
    if(nu==mu) continue;
    for(y=1;y<glb_size[nu];y++) {

      wilson_fund.re=wilson_fund.im=0.0;
      wilson_adj=0.0;
      _MASTER_FOR(&glattice,ix) {
        ix1=multistep(ix,nu,y);
        _suNg_times_suNg(tmp,*_4FIELD_AT(gfixed,ix,nu),*_4FIELD_AT(gfixed,ix1,mu));
        ix1=multistep(ix,mu,x);
        _suNg_times_suNg_dagger(tmp1,tmp,*_4FIELD_AT(gfixed,ix1,nu));
        _suNg_times_suNg_dagger(tmp,tmp1,*_4FIELD_AT(gfixed,ix,mu));
        _suNg_trace_re(ctmp.re,tmp);
        _suNg_trace_im(ctmp.im,tmp);
        _complex_add_assign(wilson_fund,ctmp);
        wilson_adj += _complex_prod_re(ctmp,ctmp)-1.;
      }
      _complex_mulr(wilson_fund,1./(GLB_T*GLB_X*GLB_Y*GLB_Z*NG),wilson_fund);
      wilson_adj /= GLB_T*GLB_X*GLB_Y*GLB_Z*(NG*NG-1);
      
      lprintf("WILSON LOOPS",0,"(mu,x,nu,y,f.re,f.im,adj) = %d %d %d %d %.8e %.8e %.8e\n",mu,x,nu,y,wilson_fund.re,wilson_fund.im,wilson_adj);
      
      if(y==glb_size[nu]-1) break;

      _MASTER_FOR(&glattice,ix) {
        ix1=multistep(ix,nu,y);
        _suNg_times_suNg(tmp,*_4FIELD_AT(gfixed,ix,nu),*_4FIELD_AT(g,ix1,nu));
        *_4FIELD_AT(gfixed,ix,nu)=tmp;
      }

    }
  }
  
  free_gfield(gfixed);
}


static int perm_abc[6][3]={{0,1,2},{1,2,0},{2,0,1},{2,1,0},{1,0,2},{0,2,1}};
static int perm_aab[3][3]={{0,1,2},{1,2,0},{2,0,1}};
static int perm_abb[3][3]={{0,1,2},{1,2,0},{2,0,1}};
static int perm_aaa[1][3]={{0,1,2}};


void ara_temporalwilsonloops(const int t, const int cc[3], const int nsteps, const suNg_field* g) {
  int i, j, k, m;
  int c[3], d[3];
  _DECLARE_INT_ITERATOR(ix);
  int iy;
  int* iz;
  suNg tmp[2];
  int max_size;
  int* perm=NULL;
  int n, nperms=0;
  int rotinv;
  int nave;
  int p;
  int chi1, chi2;

  suNg *t_llinks=malloc(sizeof(suNg)*glattice.gsize);
  suNg *s_llinks=malloc(sizeof(suNg)*glattice.gsize);
  
  error(g->type!=&glattice,1,"wilsonloops.c","The gauge field must be defined on the whole lattice");
  error(t<=0,1,"wilsonloops.c","The length must be > 0");

  if(GLB_X==GLB_Y && GLB_Y==GLB_Z && bc[1]==bc[2] && bc[2]==bc[3]) {
    rotinv=(1==1);
    lprintf("ARA WILSON LOOPS",0,"Spatial rotational invariance found.\n");
  } else
    lprintf("ARA WILSON LOOPS",0,"No spatial rotational invariance found.\n");

  max_size=(GLB_X-2>GLB_Y-2)?GLB_X-2:GLB_Y-2;
  max_size=(GLB_Z-2>max_size)?GLB_Z-2:max_size;


  _MASTER_FOR(&glattice,ix) {
    t_llinks[ix]=*_4FIELD_AT(g,ix,0);
    iy=ix;
    for(i=1;i<t;i++) {
      iy=iup(iy,0);
      _suNg_times_suNg(tmp[0],t_llinks[ix],*_4FIELD_AT(g,iy,0));
      t_llinks[ix]=tmp[0];
    }
  }

  d[0]=d[1]=d[2]=0;
  c[0]=(cc[0]>=0)?cc[0]:-cc[0];
  c[1]=(cc[1]>=0)?cc[1]:-cc[1];
  c[2]=(cc[2]>=0)?cc[2]:-cc[2];

  for(i=0;i<3;i++)
  for(j=0;j<3;j++) {
    if(i==j) continue;
    k=3-i-j;
    if(c[i]>=c[j] && c[j]>=c[k]) {
      d[0]=c[i];d[1]=c[j];d[2]=c[k];
    }
  }

  lprintf("ARA WILSON LOOPS",50,"d={ %d , %d , %d }\n",d[0],d[1],d[2]);

  if(d[0]>d[1] && d[1]>d[2]) {
    perm=perm_abc[0];
    nperms=6;
  } else if(d[0]==d[1] && d[1]>d[2]) {
    perm=perm_aab[0];
    nperms=3;
  } else if(d[0]>d[1] && d[1]==d[2]) {
    perm=perm_abb[0];
    nperms=3;
  } else if(d[0]==d[1] && d[1]==d[2]) {
    perm=perm_aaa[0];
    nperms=1;
  }

  lprintf("ARA WILSON LOOPS",50,"n permutations = %d\n",nperms);
  
  iz=malloc(glattice.gsize*sizeof(int));

  double wilson_adj;
  complex wilson_fund;
  complex ctmp;
  complex wilson_fund_ave[nsteps];
  double wilson_adj_ave[nsteps];

  for(p=1;p<=nsteps;p++) {
    wilson_fund_ave[p].re=wilson_fund_ave[p].im=0.0;
    wilson_adj_ave[p]=0.0;
  }
  nave=0;

  for(n=0;n<nperms;n++) {
  
    c[perm[3*n+0]]=d[0];
    c[perm[3*n+1]]=d[1];
    c[perm[3*n+2]]=d[2];
    
    if(!rotinv) {
      for(p=1;p<=nsteps;p++) {
        wilson_fund_ave[p].re=wilson_fund_ave[p].im=0.0;
        wilson_adj_ave[p]=0.0;
      }
      nave=0;
    }

    for(m=0;m<4;m++) {

      c[perm[3*n+0]]=d[0];
      c[perm[3*n+1]]=d[1];
      c[perm[3*n+2]]=d[2];

      if(m%2==1) {
        if(c[1]==0 || c[0]==0) continue;
        else c[1]=-c[1];
      }

      if(m/2==1) {
        if(c[2]==0 || (c[0]==0 && c[1]==0)) continue;
        else c[2]=-c[2];
      }

/*      lprintf("ARA WILSON LOOPS",50,"c={ %d , %d , %d }\n",c[0],c[1],c[2]);*/

      _MASTER_FOR(&glattice,ix) {
        iz[ix]=ix;
        _suNg_unit(s_llinks[ix]);
      }

      for(p=1;p<=nsteps;p++) {
      
        chi1=2*d[2]-d[0];
        chi2=2*d[1]-d[0];
        for(i=0;i<d[0];i++) {
          
          _MASTER_FOR(&glattice,ix) {
            if(c[perm[3*n+0]]>0) {
              for(i=0;i<c[perm[3*n+0]];i++) {
                _suNg_times_suNg(tmp[0],s_llinks[ix],*_4FIELD_AT(g,iz[ix],perm[3*n+0]+1));
                s_llinks[ix]=tmp[0];
                iz[ix]=iup(iz[ix],perm[3*n+0]+1);
              }
            } else {
              for(i=0;i<-c[perm[3*n+0]];i++) {
                iz[ix]=idn(iz[ix],perm[3*n+0]+1);
                _suNg_times_suNg_dagger(tmp[0],s_llinks[ix],*_4FIELD_AT(g,iz[ix],perm[3*n+0]+1));
                s_llinks[ix]=tmp[0];
              }
            }
          }
          
          if(chi2>=0 && d[1]!=0) {
            chi2=chi2-2*d[0];
            
            _MASTER_FOR(&glattice,ix) {
              if(c[perm[3*n+1]]>0) {
                for(i=0;i<c[perm[3*n+1]];i++) {
                  _suNg_times_suNg(tmp[0],s_llinks[ix],*_4FIELD_AT(g,iz[ix],perm[3*n+1]+1));
                  s_llinks[ix]=tmp[0];
                  iz[ix]=iup(iz[ix],perm[3*n+1]+1);
                }
              } else {
                for(i=0;i<-c[perm[3*n+1]];i++) {
                  iz[ix]=idn(iz[ix],perm[3*n+1]+1);
                  _suNg_times_suNg_dagger(tmp[0],s_llinks[ix],*_4FIELD_AT(g,iz[ix],perm[3*n+1]+1));
                  s_llinks[ix]=tmp[0];
                }
              }
            }
            
          }
          
          if(chi1>=0 && d[2]!=0) {
            chi1=chi1-2*d[0];
            
            _MASTER_FOR(&glattice,ix) {
              if(c[perm[3*n+2]]>0) {
                for(i=0;i<c[perm[3*n+2]];i++) {
                  _suNg_times_suNg(tmp[0],s_llinks[ix],*_4FIELD_AT(g,iz[ix],perm[3*n+2]+1));
                  s_llinks[ix]=tmp[0];
                  iz[ix]=iup(iz[ix],perm[3*n+2]+1);
                }
              } else {
                for(i=0;i<-c[perm[3*n+2]];i++) {
                  iz[ix]=idn(iz[ix],perm[3*n+2]+1);
                  _suNg_times_suNg_dagger(tmp[0],s_llinks[ix],*_4FIELD_AT(g,iz[ix],perm[3*n+2]+1));
                  s_llinks[ix]=tmp[0];
                }
              }
            }
            
          }

        
        }
        
        /* compute the Wilson loops */
        wilson_fund.re=wilson_fund.im=0.;
        wilson_adj=0.;
        
        _MASTER_FOR(&glattice,ix) {
          _suNg_times_suNg(tmp[0],s_llinks[ix],t_llinks[iz[ix]]);
          _suNg_times_suNg_dagger(tmp[1],tmp[0],s_llinks[multistep(ix,0,t)]);
          _suNg_times_suNg_dagger(tmp[0],tmp[1],t_llinks[ix]);
          _suNg_trace_re(ctmp.re,tmp[0]);
          _suNg_trace_im(ctmp.im,tmp[0]);
          _complex_add_assign(wilson_fund,ctmp);
          wilson_adj += _complex_prod_re(ctmp,ctmp)-1.;
        }
        _complex_mulr(wilson_fund,1./(GLB_T*GLB_X*GLB_Y*GLB_Z*NG),wilson_fund);
        wilson_adj /= GLB_T*GLB_X*GLB_Y*GLB_Z*(NG*NG-1);
        
        nave++;
        wilson_fund_ave[p].re+=wilson_fund.re;
        wilson_fund_ave[p].im+=wilson_fund.im;
        wilson_adj_ave[p]+=wilson_adj;

      } /* for(p=1;p<=nsteps;p++) */
    
    } /* for(m=0;m<4;m++) */

    if(!rotinv) {
      for(p=1;p<=nsteps;p++) {
        wilson_fund_ave[p].re/=nave;
        wilson_fund_ave[p].im/=nave;
        wilson_adj_ave[p]/=nave;

        c[0]=(c[0]>0)?c[0]:-c[0];
        c[1]=(c[1]>0)?c[1]:-c[1];
        c[2]=(c[2]>0)?c[2]:-c[2];
        lprintf("ARA WILSON LOOPS",0,"(T,dx,dy,dz,R,f.re,f.im,adj) = %d %d %d %d %.8e %.8e %.8e %.8e\n",
                t,p*c[0],p*c[1],p*c[2],p*sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]),
                wilson_fund_ave[p].re,wilson_fund_ave[p].im,
                wilson_adj_ave[p]);

      }
    }

  } /* for(n=0;n<nperms;n++) */

  if(rotinv) {
    for(p=1;p<=nsteps;p++) {
      wilson_fund_ave[p].re/=nave;
      wilson_fund_ave[p].im/=nave;
      wilson_adj_ave[p]/=nave;

      c[0]=(c[0]>0)?c[0]:-c[0];
      c[1]=(c[1]>0)?c[1]:-c[1];
      c[2]=(c[2]>0)?c[2]:-c[2];
      lprintf("ARA WILSON LOOPS",0,"(T,dx,dy,dz,R,f.re,f.im,adj) = %d %d %d %d %.8e %.8e %.8e %.8e\n",
              t,p*d[0],p*d[1],p*d[2],p*sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]),
                wilson_fund_ave[p].re,wilson_fund_ave[p].im,
                wilson_adj_ave[p]);

    }
  }

  
  free(t_llinks);
  free(s_llinks);
  free(iz);
}



