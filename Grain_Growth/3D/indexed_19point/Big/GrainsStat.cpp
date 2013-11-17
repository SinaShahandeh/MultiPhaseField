#include <iostream>

using namespace std;
int GrainsStat(int* Inds, int* MGsize, int nuclein, int mboxsize, int nboxsize, int kboxsize)
{
  int ni, i, j, k, jn, kn;
  
  #pragma omp parallel for
  for (ni=1;ni<nuclein+1;ni++)
  {
    //search for ni in the Inds matrix
    MGsize[ni]=0;
    for (i=0;i<mboxsize;i++)
    {
      for(j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for(k=0;k<kboxsize;k++)
        {
          kn=k*mboxsize*nboxsize;
          if(ni==Inds[i+jn+kn])
          {
            MGsize[ni]=MGsize[ni]+1;
          }
        }
      }
    }
  }

  
  //count grains with larger size than a thereshold
  int nzgrains=0;
  for (ni=1;ni<nuclein+1;ni++)
  {
    if (MGsize[ni]>20)
    {
      nzgrains=nzgrains+1;
    }
  }

  //finding average grain size
  int meanG=0;
  for (ni=1;ni<nuclein+1;ni++)
  {
    if (MGsize[ni]>20)
      meanG=meanG+MGsize[ni];
  }
  meanG=meanG/nzgrains;

  /*  
  int* MGnonzero;
  MGnonzero=new int[ninz+1];
  for (ni=1;ni<ninz+1;ni++)
  {
    MGnonzero[ni]=MG[ni];
}
*/
  
  return meanG;
}



