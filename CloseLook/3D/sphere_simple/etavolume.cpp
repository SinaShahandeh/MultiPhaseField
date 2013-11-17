//calculates volume of the first order parameter in the eta matrix
double etavolume(double* eta,int mboxsize, int nboxsize, int kboxsize)
{
  int i, j, k, kn, jn;
  int size2=mboxsize*nboxsize;
  double vol=0.00;
  
  for (k=0;k<kboxsize;k++)
  {
    kn=k*size2;
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
        vol=vol+eta[i+jn+kn];
      }
    }
  }
  return vol;
}