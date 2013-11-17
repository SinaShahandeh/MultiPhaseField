#include <math.h>
// #include "sphere.h"

double *sphere(int radius)
{
double *matsphe;
matsphe= new double[2*radius*2*radius];
int i,j;

  for(j=0;j<2*radius;j++)
  {
    for(i=0;i<2*radius;i++)
    {
      if (((i-radius)*(i-radius)+(j-radius)*(j-radius))<radius*radius)
      {
        matsphe[i+j*2*radius]=1;
      }
      else
      {
        matsphe[i+j*2*radius]=0;
      }
    }
  }

return matsphe;
}
