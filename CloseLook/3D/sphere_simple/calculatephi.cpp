
double* calculatephi(double* eta, double* phi, int mboxsize,int nboxsize, int kboxsize, int p)
{
   int i, j, k, kn, jn, pnn, pind;
   int size2=mboxsize*nboxsize;
   int size3=mboxsize*nboxsize*kboxsize;
      // making the phi array
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            phi[i+jn+kn]=0;
            for (pind=0;pind<p;pind++)
            {
              pnn=pind*size3;
              phi[i+jn+kn]=phi[i+jn+kn]+eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn];
            }
          }
        }
      }
return phi;
}
