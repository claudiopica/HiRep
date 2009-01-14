#include <qdp-config.h>

QLA_Real
plaquette(QDP_ColorMatrix *link[])
{
  int mu, nu;
  QLA_Real plaq, total;
  QDP_ColorMatrix *tmp1;
  QDP_ColorMatrix *tmp2;
  QDP_ColorMatrix *tmp3;
  QDP_ColorMatrix *tmp4;

  total = 0;

  tmp1 = QDP_create_M();
  tmp2 = QDP_create_M();
  tmp3 = QDP_create_M();
  tmp4 = QDP_create_M();
  

/*   QDP_extract_C(tmp2,tmp1,QDP_all); */

  for (mu = 0; mu < NDIM; mu++) {
    for (nu = mu + 1; nu < NDIM; nu++) {
      QDP_M_eq_sM(tmp1, link[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
      QDP_M_eq_sM(tmp2, link[mu], QDP_neighbor[nu], QDP_forward, QDP_all);
      QDP_M_eq_Ma_times_M(tmp3, link[nu], link[mu], QDP_all);
      QDP_M_eq_M_times_M(tmp4, tmp3, tmp1, QDP_all);
      QDP_r_eq_re_M_dot_M(&plaq, tmp2, tmp4, QDP_all);
      total += plaq;
    }
  }

  QDP_destroy_M(tmp1);
  QDP_destroy_M(tmp2);
  QDP_destroy_M(tmp3);
  QDP_destroy_M(tmp4);

  return total;
}

