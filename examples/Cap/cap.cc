#include <stdio.h>
#include <complex>
using namespace std;

int main()
{
  double PI = 3.14159265358979323846;
  double eps= 8.854187818e-12*3.9;
  double s  = 1e-6*1e-6;
  double d  = 1e-6;
  double C  = eps*s/d;
  double R  = 1000;
  double U  = 0.0026;
  double f;
  
  FILE *fp=fopen("ivac.ana.txt", "w");
  for(f=1e6;f<1e7;f*=1.1)
  {
    complex<double> Z(R,-1.0/(2*PI*f*C));
    complex<double> I=1e3*U/Z; 
    fprintf(fp,"%e %e %e\n",f/1e6,I.real(),I.imag());
  }
  fclose(fp);  

} 
