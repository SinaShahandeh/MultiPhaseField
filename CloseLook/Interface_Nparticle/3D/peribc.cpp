

//getting an index of the domain and returning outindex considering \
periodic boundary condition

// #include "peribc.h"
int peribc(int inindex, int& size)
{
int outindex;
outindex=inindex;
if(inindex<0)
    {
     outindex=inindex+size;
    }
if (inindex>=size)
   {
    outindex=inindex-size;
   }
return (outindex);
}
