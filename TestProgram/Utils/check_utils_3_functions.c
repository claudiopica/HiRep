/*This is an automatically generated function, do not edit.*/
static int fullcheck(int rotid, double complex *rotated, double complex *unrotated)
{
#define rotfun(a) rotated[(a)]
#define unrotfun(a) unrotated[(a)]
double complex tmp[2];
int return_value=0;
if(0==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,1)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,1));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(MapOptoCindex(1,0,0,1,-1,2)) - 1.*unrotfun(MapOptoCindex(1,0,0,1,-1,2));
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(0) - 1.*unrotfun(0);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(1) - 1.*unrotfun(1);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(2) - 1.*unrotfun(2);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(20==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(19==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(23==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(27==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(28==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(46==rotid){
 tmp[0]=rotfun(3) - 1.*unrotfun(3);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(0==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(1==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(26==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
if(45==rotid){
 tmp[0]=rotfun(4) - 1.*unrotfun(4);
 _complex_mul_star(tmp[1],tmp[0],tmp[0]);
 if(sqrt(creal(tmp[1]))>1.e-14) {
return_value++;
 }
}
#undef unrotfunreturn
#undef rotfun
return return_value;
}
