

//getting an index of the domain and returning outindex considering \
symmetric boundary condition

// #include "symbc.h"
int symbc(int inindex, int& size)
{
int outindex;
outindex=inindex;
if(inindex<0)
    {
     outindex=-inindex;
    }
if (inindex>=size)
   {
    outindex=2*size-inindex-1;
   }
return (outindex);
}
