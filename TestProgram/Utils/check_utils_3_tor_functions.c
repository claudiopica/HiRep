/*This is an automatically generated function, do not edit.*/
#define Complex(a,b) ((a)+I*(b))
static int fulltorcheck(int rotid, double complex *rotated, double complex *unrotated)
{
#define rotfun(a) rotated[(a)]
#define unrotfun(a) unrotated[(a)]
#if total_n_tor_op>0
double complex tmp[3];
#endif
int return_value=0;
if(0==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(1==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(2==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(3==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(4==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(5==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(6==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(7==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(8==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(9==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(10==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(11==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(12==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(13==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(14==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(15==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(16==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(17==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(18==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(19==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(20==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(21==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(22==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(23==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(24==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(25==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(26==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(27==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(28==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(29==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(30==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(31==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(32==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(33==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(34==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(35==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(36==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(37==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(38==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(39==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(40==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(41==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(42==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(43==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(44==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(45==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(46==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(47==rotid)
{
tmp[0]=rotfun(0);
tmp[1]=1.*unrotfun(0);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=4 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[az]][0, 0, 0] Pim[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pim[L[ax]][0, 0, 0] Pim[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pim[L[ay]][0, 0, 0] Pim[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ax, az, az, ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[ax, ay, ay, -ax, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ax, az, az, -ax, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[-ay, az, az, ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[ay, ax, ax, -ay, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[az]][0, 0, 0] Pre[ay, az, az, -ay, -2 az + L[az]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[-az, ax, ax, az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[-az, ay, ay, az, -2 ay + L[ay]][0, 0, 0] + 4 Pre[L[ax]][0, 0, 0] Pre[az, ax, ax, -az, -2 ax + L[ax]][0, 0, 0] + 4 Pre[L[ay]][0, 0, 0] Pre[az, ay, ay, -az, -2 ay + L[ay]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",0,0,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(0==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(20==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(19==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(23==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(27==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(28==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(45==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(46==rotid)
{
tmp[0]=rotfun(1);
tmp[1]=1.*unrotfun(1);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[-ay, ax, ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[ay, ax, -ay, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[-az, ax, az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] + 2 Pim[az, ax, -az, -ax + L[ax]][0, 0, 0] Pre[L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-ay, ax, ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[ay, ax, -ay, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[-az, ax, az, -ax + L[ax]][0, 0, 0] - 2 Pim[L[ax]][0, 0, 0] Pre[az, ax, -az, -ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,0,0,1,1,-1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(0==rotid)
{
tmp[0]=rotfun(4);
tmp[1]=1.*unrotfun(4);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 2 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,1,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(1==rotid)
{
tmp[0]=rotfun(4);
tmp[1]=1.*unrotfun(4);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 2 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,1,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(26==rotid)
{
tmp[0]=rotfun(4);
tmp[1]=1.*unrotfun(4);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 2 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,1,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
if(45==rotid)
{
tmp[0]=rotfun(4);
tmp[1]=1.*unrotfun(4);
_complex_mul_star(tmp[2],tmp[0]-tmp[1],tmp[0]-tmp[1]);
if(sqrt(creal(tmp[2]))>=1.e-10){
  lprintf("Error",0," Tor=2 Pim[L[ay]][0, 0, 0] Pim[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pim[L[ax]][0, 0, 0] Pim[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0] + 2 Pre[L[ay]][0, 0, 0] Pre[-ax, ay, ay, ax, -2 ay + L[ay]][0, 0, 0] + 2 Pre[L[ax]][0, 0, 0] Pre[-ay, ax, ax, ay, -2 ax + L[ax]][0, 0, 0]\n px=%d py=%d pz=%d Irrep=%d ev=%d charge=%d multiplet id=%d (%2.10e %2.10e) (%2.10e %2.10e) %2.10e \n",1,1,0,1,1,1,1,creal(tmp[0]),cimag(tmp[0]),creal(tmp[1]),cimag(tmp[1]),sqrt(creal(tmp[2])));
  return_value++;
  }
}
#undef unrotfunreturn
#undef rotfun
#undef Complex
return return_value;
}
