#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <list>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>
#include "type.h"
#include "utils.h"
#include "read_arg.h"
#include "matrix_algebra.h"
#include <stdlib.h>     /* srand, rand */
#include <iomanip>      // std::setprecision

using namespace std;

void print(const double * a, const int n);

int main(){
  int n=5;

  double *a=new double[n*n];
  
  double *b=new double[n*n];
  
  double *binv=new double[n*n];
  
  double *c=new double[n*n];
 
  double *wi=new double[n];
  double *wr=new double[n];
  
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) c[j+n*i]=b[j+n*i]=a[i*n+j]=rand()%10;
  
  print(a,n);
  
  int idx[n];
  double dd;

  ludcmp(a, n, idx, &dd);

  double det = 1.0;
  
  for(int i=0;i<n;i++) det *= a[i*n+i];
  
  cout << endl << "Det = " << std::setprecision(20) << det*dd << endl;

  cout << endl << "function Det = " << std::setprecision(20) << Det(b,n) << endl<< endl;


  inverse(b, binv, n);
  cout << endl<<"Evaluating the inverse matrix" << endl<< endl;
  
  print(binv,n);
  
  for(int i=0;i<n*n;i++) a[i]=0.;
  
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      for(int k=0;k<n;k++)
	a[i*n+j]+=b[i*n+k]*binv[k*n+j];
  
  
  for(int i=0;i<n;i++) a[i*(n+1)]-=1.0;
  
  double sum=0.;
  for(int i=0;i<n*n;i++) sum+=a[i]*a[i];
  
  cout << endl <<"Distance of m.m1 from the identity " << sqrt(sum) << endl << endl;
  cout << endl <<"Evaluating the Hessenberg form " << endl<< endl;
  
  elmhes(c, n);
  
  for (int j=0;j<n-2;j++)
    for (int i=j+2;i<n;i++) c[i*n+j]=0.0;
  
  print(c,n);

  cout << endl <<"Evaluating the Eigenvalues" << endl<< endl;
	
	hqr(c,n,wr,wi);
	
	for(int j=0;j<n;j++) cout << wr[j]<< " "<< wi[j] <<endl;
	


 return 0;
}

void print(const double * a, const int n)
{
	cout << "m={";
  for(int i=0;i<n;i++){
    cout << "{";
    for(int j=0;j<n;j++) {
      cout << a[i*n+j];
      if(j!=n-1) cout <<", ";
    }
    cout << "}";
    if(i!=n-1) cout<<",";
    else  cout<<"};";
    cout << endl;
  }
	
}
