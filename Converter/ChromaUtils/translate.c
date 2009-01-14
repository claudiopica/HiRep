#include <qdp-config.h>
#include "suN.h"
#include "global.h"

static int nu;

void translate(QLA_ColorMatrix *m, int coords[])
{
  int i,j;
  QLA_Complex t2;
  suNg tmpmat;
  for(i=0; i<QDP_Nc; i++) {
    for(j=0; j<QDP_Nc; j++) {
      t2 = QLA_elem_M(*m,i,j);

      tmpmat.c[i*QDP_Nc+j].re=(double)t2.real;
      tmpmat.c[i*QDP_Nc+j].im=(double)t2.imag;

      *pu_gauge(ipt(coords[3],coords[0],coords[1],coords[2]),nu)=tmpmat;
    }
  }
}
 


void chroma_to_HiRep(QDP_ColorMatrix *link[])
{
  int mu;
  printf0("\n\nParametri vari: Nc=%d\n",QDP_Nc);
  printf0("Lattice size: (T=%d X=%d Y=%d Z=%d)\n\n",QDP_coord_size(3),QDP_coord_size(0),QDP_coord_size(1),QDP_coord_size(2));

  for (mu = 0; mu < NDIM; mu++) {
    nu = (mu+1)%4;
    QDP_M_eq_func(link[mu],  translate, QDP_all);
  }
}
