#include <math.h>
#include <iostream>
#include <iomanip>

#include "type.h"
#include "fit.h"

void fitexp(par *apar, double *diagzed_corr)
{
  double *acorn = new double[apar->numbinjk * apar->ndt];
  double *acor = new double[apar->ndt];
  double *scor = new double[apar->ndt];
  double *energy = new double[apar->numbinjk];
  double *chisq = new double[apar->numbinjk];
  double *coef = new double[apar->numbinjk];

  double anorm = 0.0, snorm = 0.0, amm, amu, aml, avn, avd, ff, etry, el, eu, abest, chis, chi, ebest, cof, sen, aen, scoef, acoef, schisq, achisq;
  int numdel;
  double small = 0.00001;

  for (int i5 = 0; i5 < apar->numbinjk; i5++)
  {
    anorm += diagzed_corr[apar->ndt * i5];
    snorm += diagzed_corr[apar->ndt * i5] * diagzed_corr[apar->ndt * i5];
  }
  anorm /= (apar->numbinjk);
  snorm = (snorm - anorm * anorm * apar->numbinjk);
  if (snorm > 0.0 || abs(snorm) < 1.e-13)
    snorm = sqrt(abs(snorm));
  else
    std::cerr << "[ERROR][fitexp] Ill defined variance of the operator" << std::endl;

  for (int i4 = 0; i4 < apar->ndt; i4++)
  {
    acor[i4] = 0.0;
    scor[i4] = 0.0;
    for (int ib = 0; ib < apar->numbinjk; ib++)
      acorn[ib + apar->numbinjk * i4] = 0;
  }

  std::cout << "[INFO][fitexp] Correlator " << std::endl;

  for (std::map<int, dtcor>::iterator it_dt = apar->corrdef.begin(); it_dt != apar->corrdef.end(); ++it_dt)
  {
    int nt = it_dt->second.index;
    for (int ib = 0; ib < apar->numbinjk; ib++)
    {
      acorn[nt + apar->ndt * ib] = diagzed_corr[nt + apar->ndt * ib] / anorm;
      acor[nt] += acorn[nt + apar->ndt * ib];
      scor[nt] += acorn[nt + apar->ndt * ib] * acorn[nt + apar->ndt * ib];
    }

    acor[nt] /= apar->numbinjk;
    scor[nt] = (scor[nt] - acor[nt] * acor[nt] * apar->numbinjk);

    if (scor[nt] > 0.0 || abs(scor[nt]) < 1.e-13)
      scor[nt] = sqrt(abs(scor[nt]));
    else
      std::cerr << "[ERROR][fitexp]Ill defined correlator for distance " << it_dt->first << std::endl;

    std::cout << "\t cor[ " << it_dt->first << " ]= " << acor[nt] * anorm << " +/- " << (apar->numbinjk - 1) / (apar->numbinjk * sqrt(apar->numbinjk)) * scor[nt] << std::endl;
  }

  if (scor[0] < small)
    scor[0] = scor[1];
  for (int nt = 0; nt < apar->ndt; nt++)
    if (scor[nt] < small)
      scor[nt] = 1.0;

  std::cout << std::endl
            << "[INFO][fitexp] Mass fit " << std::endl;

  int iOBC = 0;
  if (apar->Tflag == 1)
    iOBC = 1;

  std::map<int, dtcor>::iterator it_dt = apar->corrdef.begin();
  std::map<int, dtcor>::iterator it_dtp1 = apar->corrdef.begin();
  it_dtp1++;

  for (; it_dtp1 != apar->corrdef.end(); ++it_dt, ++it_dtp1)
  {
    int itl = it_dt->second.index;
    int itlp1 = it_dtp1->second.index;
    int dist = it_dtp1->first - it_dt->first;

    if (acor[itl] > small && acor[itlp1] > small && acor[itlp1] - scor[itlp1] > small && acor[itl] - scor[itl] > small)
    {
      amm = log(acor[itl] / acor[itlp1]) / dist;
      amu = log((acor[itl] + scor[itl]) / (acor[itlp1] - scor[itlp1])) / dist;
      aml = log((acor[itl] - scor[itl]) / (acor[itlp1] + scor[itlp1])) / dist;
      el = amm - 3.0 * (amm - aml);
      eu = amm + 3.0 * (amu - amm);
      numdel = 1000;
    }
    else
    {
      if (acor[1] > small && acor[0] > small)
        eu = 2.0 * log(acor[0] / acor[1]);
      else
        eu = 5.0;
      el = 0.0;
      numdel = 10000;
    }
    
    if (eu > 5.0)
      eu = 5.0;
    if (el < 0.0)
      el = 0.0;

    for (int itu = itl + 1; itu < apar->ndt; itu++)
    {
      for (int ib = 0; ib < apar->numbinjk; ib++)
      {
        for (int idel = 0; idel < numdel; idel++)
        {
          etry = el + idel * ((eu - el) / numdel);
          avn = 0.0;
          avd = 0.0;
          for (int idt = itl; idt < itu + 1; idt++)
          {
            ff = exp(-etry * idt) + iOBC * exp(-etry * (apar->lt - idt));
            avn += ff * acorn[idt + apar->ndt * ib] / (scor[idt] * scor[idt]);
            avd += ff * ff / (scor[idt] * scor[idt]);
          }
          if (avd < small)
            avd = small;
          abest = avn / avd;
          chis = 0.0;
          for (int idt = itl; idt < itu + 1; idt++)
          {
            ff = abest * (exp(-etry * idt) + iOBC * exp(-etry * (apar->lt - idt)));
            chis += ((acorn[idt + apar->ndt * ib] - ff) / scor[idt]) * ((acorn[idt + apar->ndt * ib] - ff) / scor[idt]);
          }
          if (idel == 0)
          {
            chi = chis;
            ebest = etry;
            cof = abest;
          }
          else if (chis < chi)
          {
            chi = chis;
            ebest = etry;
            cof = abest;
          }
        }
        energy[ib] = ebest;
        chisq[ib] = chi;
        coef[ib] = cof;
      }

      sen = aen = scoef = acoef = schisq = achisq = 0.0;
      for (int ib = 0; ib < apar->numbinjk; ib++)
      {
        aen += energy[ib];
        sen += energy[ib] * energy[ib];
        acoef += coef[ib];
        scoef += coef[ib] * coef[ib];
        achisq += chisq[ib];
        schisq += chisq[ib] * chisq[ib];
      }
      aen /= apar->numbinjk;
      sen /= apar->numbinjk;
      sen = (apar->numbinjk * (sen - aen * aen));
      acoef /= apar->numbinjk;
      scoef /= apar->numbinjk;
      scoef = (apar->numbinjk * (scoef - acoef * acoef));
      achisq /= apar->numbinjk;
      schisq /= apar->numbinjk;
      schisq = (apar->numbinjk * (schisq - achisq * achisq));

      if (sen > 0)
        sen = sqrt(sen);
      if (scoef > 0)
        scoef = sqrt(scoef);
      if (schisq > 0)
        schisq = sqrt(schisq);

      std::cout << std::fixed << std::setprecision(6);
      std::cout << "dt range " << itl << "-" << itu << "\tmass= " << aen << " " << sen << " \tamplitude= " << acoef << " " << scoef << "\t\tchisq= " << achisq << " " << schisq << std::endl;
      if (achisq > 20)
        itu = apar->ndt;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  delete[] acorn;
  delete[] acor;
  delete[] scor;
  delete[] energy;
  delete[] chisq;
  delete[] coef;
}
