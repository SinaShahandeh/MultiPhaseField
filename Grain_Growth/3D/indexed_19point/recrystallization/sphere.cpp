#include <math.h>
// #include "sphere.h"

double *sphere(int radius)
{
double *matsphe;
matsphe= new double[(2*radius+2)*(2*radius+2)*(2*radius+2)];
int i,j,k;
for (k=0;k<2*radius+3;k++)
{
  for(j=0;j<2*radius+3;j++)
  {
    for(i=0;i<2*radius+3;i++)
    {
      if (((i-radius-1)*(i-radius-1)+(j-radius-1)*(j-radius-1)+(k-radius-1)*(k-radius-1))<radius*radius)
      {
        matsphe[i+j*2*radius+k*2*radius*2*radius]=1;
      }
      else
      {
        matsphe[i+j*2*radius+k*2*radius*2*radius]=0;
      }
    }
  }
}
return matsphe;
}
