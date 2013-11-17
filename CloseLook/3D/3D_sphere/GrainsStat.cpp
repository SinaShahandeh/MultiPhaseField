using namespace std;
int GrainsStat(int* Inds, int* MGsize, int nuclein, int mboxsize, int nboxsize, int kboxsize)
{
  int ni, i, j, k, jn, kn;
  
  /*
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
  */
  ni=9;  
  int G;
  //search for ni in the Inds matrix
  G=0;
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
          G=G+1;
        }
      }
    }
  }
  
  return G;
}



