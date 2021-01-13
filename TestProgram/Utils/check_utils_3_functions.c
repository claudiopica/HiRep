/*This is an automatically generated function, do not edit.*/
#define Complex(a, b) ((a) + I * (b))
static int fullcheck(int rotid, double complex *rotated, double complex *unrotated)
{
#define rotfun(a) rotated[(a)]
#define unrotfun(a) unrotated[(a)]
  double complex tmp[3];
  int return_value = 0;
  if (0 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(0);
    tmp[1] = 1. * unrotfun(0);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, az, -ay, az, -ax, -ay, -az, -az][-1, 0, 0] + Pim[ax, az, az, ay, -az, ay, -ax, -az, -ay, -ay][-1, 0, 0] + Pim[ay, ay, ax, az, az, -ay, -az, -ay, -ax, -az][-1, 0, 0] + Pim[ay, az, ay, -ax, az, -ay, -ay, ax, -az, -az][0, 0, 0] + Pim[ay, az, az, ax, -ay, -ay, -az, ay, -az, -ax][-1, 0, 0] + Pim[az, ay, ay, ax, -az, -az, -ay, az, -ay, -ax][-1, 0, 0] + Pim[az, ay, az, -ax, ay, -az, -az, ax, -ay, -ay][0, 0, 0] + Pim[az, az, ax, ay, ay, -az, -ay, -az, -ax, -ay][-1, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(1);
    tmp[1] = 1. * unrotfun(1);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][-1, 0, 0] + Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][-1, 0, 0] + Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][-1, 0, 0] - Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, 0, 0] - Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][-1, 0, 0] + Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][-1, 0, 0] - Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", -1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(2);
    tmp[1] = 1. * unrotfun(2);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, -1, 0] + Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, -1, 0] + Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, -1, 0] + Pim[ay, az, az, az, -ax, -az, -ay, -az, ax, -az][0, -1, 0] + Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ax, az, ay, -az, -az, -az, -ax, az, -ay][0, -1, 0] - Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, -1, 0] - Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(3);
    tmp[1] = 1. * unrotfun(3);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ay, az, az, -ax, -az, -ax, -ay, -az][0, -1, 0] + Pim[ax, az, ax, -ay, az, -ax, -ax, ay, -az, -az][0, 0, 0] + Pim[ax, az, az, ay, -ax, -ax, -az, ax, -az, -ay][0, -1, 0] + Pim[ay, ax, ax, az, -ax, az, -ay, -ax, -az, -az][0, -1, 0] + Pim[ay, az, az, ax, -az, ax, -ay, -az, -ax, -ax][0, -1, 0] + Pim[az, ax, ax, ay, -az, -az, -ax, az, -ax, -ay][0, -1, 0] + Pim[az, ax, az, -ay, ax, -az, -az, ay, -ax, -ax][0, 0, 0] + Pim[az, az, ay, ax, ax, -az, -ax, -az, -ay, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, -1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(4);
    tmp[1] = 1. * unrotfun(4);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, -1] - Pim[ax, ay, ay, ay, -az, -ay, -ax, -ay, az, -ay][0, 0, 0] - Pim[ax, az, ax, ay, -ax, -ax, -ax, -az, ax, -ay][0, 0, -1] - Pim[ay, ax, ax, ax, -az, -ax, -ay, -ax, az, -ax][0, 0, 0] - Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, -1] - Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] - Pim[ay, az, ay, ax, -ay, -ay, -ay, -az, ay, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(5);
    tmp[1] = 1. * unrotfun(5);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, az, ay, ay, -ax, -ay, -ax, -az, -ay][0, 0, -1] + Pim[ax, ay, ax, -az, ay, -ax, -ax, az, -ay, -ay][0, 0, 0] + Pim[ax, ay, ay, az, -ax, -ax, -ay, ax, -ay, -az][0, 0, -1] + Pim[ay, ax, ax, az, -ay, -ay, -ax, ay, -ax, -az][0, 0, -1] + Pim[ay, ax, ay, -az, ax, -ay, -ay, az, -ax, -ax][0, 0, 0] + Pim[ay, ay, az, ax, ax, -ay, -ax, -ay, -az, -ax][0, 0, -1] + Pim[az, ax, ax, ay, -ax, ay, -az, -ax, -ay, -ay][0, 0, -1] + Pim[az, ay, ay, ax, -ay, ax, -az, -ay, -ax, -ax][0, 0, -1]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, -1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (1 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (1 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (1 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (2 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (2 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (2 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (3 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (3 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (3 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (4 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (4 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (4 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (5 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (5 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (5 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (6 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (6 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (6 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (7 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (7 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (7 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (8 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (8 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (8 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (9 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (9 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (9 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (10 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (10 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (10 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (11 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (11 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (11 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (12 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (12 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (12 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (13 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (13 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (13 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (14 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (14 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (14 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (24 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (24 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (24 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (31 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (31 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (31 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (32 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (32 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (32 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (33 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (33 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (33 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (34 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (34 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (34 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (35 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (35 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (35 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (36 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (36 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (36 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (37 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., -0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (37 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (37 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (38 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (38 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (38 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = Complex(0., 0.5) * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - Complex(0., 0.5) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (39 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (39 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (39 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (40 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - Complex(0., 1.) * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (40 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (40 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + Complex(0., 1.) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (41 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (41 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 0.7071067811865475 * unrotfun(6) + 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (41 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (42 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) - 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (42 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 0.7071067811865475 * unrotfun(6) - 0.7071067811865475 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (42 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = -0.5 * unrotfun(6) + 0.7071067811865475 * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (43 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (43 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - Complex(0., 0.7071067811865475) * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (43 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) - Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (44 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = -0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) + 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (44 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + Complex(0., 0.7071067811865475) * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (44 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0.5 * unrotfun(6) + Complex(0., 0.7071067811865475) * unrotfun(7) - 0.5 * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. - 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. - 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. - 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(6);
    tmp[1] = 0. + 1. * unrotfun(8);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_1_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(7);
    tmp[1] = 0. + 1. * unrotfun(7);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_2_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 2, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(8);
    tmp[1] = 0. + 1. * unrotfun(6);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=000_T1minusOhP+_3_001_A1Dic4-_xx-y-x-yz-xyy-z_001_A1Dic4-_zx-yx-zxy-x-x-x\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 0, 9, 3, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(9);
    tmp[1] = 1. * unrotfun(9);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -ay, -ax, -ay, az, -ax, ay, ay, -az][0, 0, 0] + Pim[ax, ax, ay, -ax, ay, az, -ax, -ay, -ay, -az][0, 0, 0] + Pim[ax, ay, ax, az, ay, -ax, -ax, -az, -ay, -ay][0, 0, 0] + Pim[ay, ax, ay, az, ax, -ay, -ay, -az, -ax, -ax][0, 0, 0] + Pim[ay, ay, -ax, -ay, -ax, az, -ay, ax, ax, -az][0, 0, 0] + Pim[ay, ay, ax, -ay, ax, az, -ay, -ax, -ax, -az][0, 0, 0] + Pim[az, ax, ay, ay, -az, -ax, -ax, -ay, ax, -ay][0, 0, 0] + Pim[az, ay, ax, ax, -az, -ay, -ay, -ax, ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (15 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (16 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (21 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (25 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(10);
    tmp[1] = 1. * unrotfun(10);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=-Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] - Pim[ax, ay, ay, ay, az, -ay, -ax, -ay, -az, -ay][0, 0, 0] - Pim[ay, ax, ax, ax, az, -ax, -ay, -ax, -az, -ax][0, 0, 0] - Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] - Pim[az, ax, -ay, ax, -az, ax, ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] - Pim[az, ay, -ax, ay, -az, ay, ax, -ay, -ay, -ay][0, 0, 0] - Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 0, 1, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(11);
    tmp[1] = 1. * unrotfun(11);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, -az, -ax, -az, ay, -ax, az, az, -ay][0, 0, 0] + Pim[ax, ax, az, -ax, az, ay, -ax, -az, -az, -ay][0, 0, 0] + Pim[ax, az, ax, ay, az, -ax, -ax, -ay, -az, -az][0, 0, 0] + Pim[ay, ax, az, az, -ay, -ax, -ax, -az, ax, -az][0, 0, 0] + Pim[ay, az, ax, ax, -ay, -az, -az, -ax, az, -ax][0, 0, 0] + Pim[az, ax, az, ay, ax, -az, -az, -ay, -ax, -ax][0, 0, 0] + Pim[az, az, -ax, -az, -ax, ay, -az, ax, ax, -ay][0, 0, 0] + Pim[az, az, ax, -az, ax, ay, -az, -ax, -ax, -ay][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (17 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (18 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (22 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (29 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (30 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (47 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(12);
    tmp[1] = 1. * unrotfun(12);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ax, ax, -az, -ax, ay, -ax, az, -ax, -ay][0, 0, 0] + Pim[ax, ax, ax, az, -ax, ay, -ax, -az, -ax, -ay][0, 0, 0] + Pim[ax, ay, ax, az, ax, -ay, -ax, -ax, -ax, -az][0, 0, 0] - Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, ay, ax, -az, ax, -ay, -ax, -ax, -ax][0, 0, 0] + Pim[az, ay, az, ax, az, -ay, -az, -az, -az, -ax][0, 0, 0] + Pim[az, az, az, -ax, -az, ay, -az, ax, -az, -ay][0, 0, 0] - Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 0, 1, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(13);
    tmp[1] = 1. * unrotfun(13);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ay, ax, ay, az, ay, -ax, -ay, -ay, -ay, -az][0, 0, 0] + Pim[ay, ay, ay, -az, -ay, ax, -ay, az, -ay, -ax][0, 0, 0] + Pim[ay, ay, ay, az, -ay, ax, -ay, -az, -ay, -ax][0, 0, 0] + Pim[ay, az, ax, az, -ay, az, -ax, -az, -az, -az][0, 0, 0] + Pim[az, ax, az, ay, az, -ax, -az, -az, -az, -ay][0, 0, 0] + Pim[az, ay, ax, ay, -az, ay, -ax, -ay, -ay, -ay][0, 0, 0] + Pim[az, az, az, -ay, -az, ax, -az, ay, -az, -ax][0, 0, 0] + Pim[az, az, az, ay, -az, ax, -az, -ay, -az, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (20 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (19 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (23 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (27 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (28 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (46 == rotid)
  {
    tmp[0] = rotfun(14);
    tmp[1] = 1. * unrotfun(14);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pim[ax, ay, az, az, -ax, -ay, -ay, -az, ay, -az][0, 0, 0] + Pim[ax, az, ay, ay, -ax, -az, -az, -ay, az, -ay][0, 0, 0] + Pim[ay, ay, -az, -ay, -az, ax, -ay, az, az, -ax][0, 0, 0] + Pim[ay, ay, az, -ay, az, ax, -ay, -az, -az, -ax][0, 0, 0] + Pim[ay, az, ay, ax, az, -ay, -ay, -ax, -az, -az][0, 0, 0] + Pim[az, ay, az, ax, ay, -az, -az, -ax, -ay, -ay][0, 0, 0] + Pim[az, az, -ay, -az, -ay, ax, -az, ay, ay, -ax][0, 0, 0] + Pim[az, az, ay, -az, ay, ax, -az, -ay, -ay, -ax][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 0, 0, 1, 1, -1, 2, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (0 == rotid)
  {
    tmp[0] = rotfun(15);
    tmp[1] = 1. * unrotfun(15);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pre[ax, -ay, az, -ax, ay, -az][-1, 1, 0] + Pre[ax, -ay, az, -ax, ay, -az][0, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][-1, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 1, 0, 4, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (1 == rotid)
  {
    tmp[0] = rotfun(15);
    tmp[1] = -1. * unrotfun(15);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pre[ax, -ay, az, -ax, ay, -az][-1, 1, 0] + Pre[ax, -ay, az, -ax, ay, -az][0, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][-1, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 1, 0, 4, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (26 == rotid)
  {
    tmp[0] = rotfun(15);
    tmp[1] = 1. * unrotfun(15);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pre[ax, -ay, az, -ax, ay, -az][-1, 1, 0] + Pre[ax, -ay, az, -ax, ay, -az][0, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][-1, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 1, 0, 4, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
  if (45 == rotid)
  {
    tmp[0] = rotfun(15);
    tmp[1] = -1. * unrotfun(15);
    _complex_mul_star(tmp[2], tmp[0] - tmp[1], tmp[0] - tmp[1]);
    if (sqrt(creal(tmp[2])) >= 1.e-11)
    {
      lprintf("Error", 0, " Op=Pre[ax, -ay, az, -ax, ay, -az][-1, 1, 0] + Pre[ax, -ay, az, -ax, ay, -az][0, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][-1, 0, 0] - Pre[ay, az, ax, -ay, -az, -ax][0, -1, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n", 1, 1, 0, 4, 1, 1, 1, creal(tmp[0]), cimag(tmp[0]), creal(tmp[1]), cimag(tmp[1]), sqrt(creal(tmp[2])));
      return_value++;
    }
  }
#undef unrotfunreturn
#undef rotfun
#undef Complex
  return return_value;
}
