/*This is an automatically generated function, do not edit.*/
static int fullcheck(int rotid, double complex *rotated, double complex *unrotated)
{
#define rotfun(a) rotated[(a)]
#define unrotfun(a) unrotated[(a)]
double complex tmp[2];
int return_value=0;
if(0==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(6==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(6==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(2==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(2==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(4==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(5) - 1.*unrotfun(5);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(4==rotid){
 tmp[0]=rotfun(5) - 1.*unrotfun(5);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(5) - 1.*unrotfun(5);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(5) - 1.*unrotfun(5);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(1==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(2==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(3==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(4==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(5==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(6==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(7==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(8==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(9==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(10==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(11==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(12==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(13==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(14==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(15==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(16==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(17==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(18==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(21==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(22==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(24==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(30==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(31==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(32==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(33==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(34==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(35==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(36==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(37==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(38==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(39==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(40==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(41==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(42==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(43==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(44==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(6) - 1.*unrotfun(6);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(1==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(2==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(3==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(4==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(5==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(6==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(7==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(8==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(9==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(10==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(11==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(12==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(13==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(14==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(15==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(16==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(17==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(18==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(21==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(22==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(24==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(30==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(31==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(32==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(33==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(34==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(35==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(36==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(37==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(38==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(39==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(40==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(41==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(42==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(43==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(44==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(7) - 1.*unrotfun(7);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(15==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(16==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(21==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(12) - 1.*unrotfun(12);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(15==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(16==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(21==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(13) - 1.*unrotfun(13);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(14) - 1.*unrotfun(14);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(4==rotid){
 tmp[0]=rotfun(14) - 1.*unrotfun(14);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(14) - 1.*unrotfun(14);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(14) - 1.*unrotfun(14);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(15) - 1.*unrotfun(15);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(4==rotid){
 tmp[0]=rotfun(15) - 1.*unrotfun(15);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(15) - 1.*unrotfun(15);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(15) - 1.*unrotfun(15);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(17==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(18==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(22==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(30==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(16) - 1.*unrotfun(16);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(17==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(18==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(22==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(30==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(17) - 1.*unrotfun(17);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(18) - 1.*unrotfun(18);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(3==rotid){
 tmp[0]=rotfun(18) - 1.*unrotfun(18);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(18) - 1.*unrotfun(18);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(18) - 1.*unrotfun(18);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(19) - 1.*unrotfun(19);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(3==rotid){
 tmp[0]=rotfun(19) - 1.*unrotfun(19);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(19) - 1.*unrotfun(19);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(47==rotid){
 tmp[0]=rotfun(19) - 1.*unrotfun(19);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(20) - 1.*unrotfun(20);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(2==rotid){
 tmp[0]=rotfun(20) - 1.*unrotfun(20);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(20) - 1.*unrotfun(20);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(20) - 1.*unrotfun(20);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(21) - 1.*unrotfun(21);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(2==rotid){
 tmp[0]=rotfun(21) - 1.*unrotfun(21);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(25==rotid){
 tmp[0]=rotfun(21) - 1.*unrotfun(21);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(21) - 1.*unrotfun(21);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(22) - 1.*unrotfun(22);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(6==rotid){
 tmp[0]=rotfun(22) - 1.*unrotfun(22);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(22) - 1.*unrotfun(22);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(22) - 1.*unrotfun(22);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(23) - 1.*unrotfun(23);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(6==rotid){
 tmp[0]=rotfun(23) - 1.*unrotfun(23);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(29==rotid){
 tmp[0]=rotfun(23) - 1.*unrotfun(23);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(23) - 1.*unrotfun(23);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(24) - 1.*unrotfun(24);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(25) - 1.*unrotfun(25);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(26) - 1.*unrotfun(26);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(5==rotid){
 tmp[0]=rotfun(26) - 1.*unrotfun(26);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(30==rotid){
 tmp[0]=rotfun(26) - 1.*unrotfun(26);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(26) - 1.*unrotfun(26);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(27) - 1.*unrotfun(27);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(5==rotid){
 tmp[0]=rotfun(27) - 1.*unrotfun(27);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(30==rotid){
 tmp[0]=rotfun(27) - 1.*unrotfun(27);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(27) - 1.*unrotfun(27);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(28) - 1.*unrotfun(28);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(1==rotid){
 tmp[0]=rotfun(28) - 1.*unrotfun(28);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(28) - 1.*unrotfun(28);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(28) - 1.*unrotfun(28);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(29) - 1.*unrotfun(29);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(1==rotid){
 tmp[0]=rotfun(29) - 1.*unrotfun(29);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(29) - 1.*unrotfun(29);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(29) - 1.*unrotfun(29);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-12) {
return_value++;
 }
}
#undef unrotfunreturn
#undef rotfun
return return_value;
}
