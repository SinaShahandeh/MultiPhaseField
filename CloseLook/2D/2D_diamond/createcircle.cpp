//create a dome of value 1 in the order paramter number 1.
double* createcircle(double *eta,int mboxsize, int nboxsize, int domeR)
{
  int i,j,jn;
  int i0=mboxsize/2;   //centre of the dome position
  int j0=nboxsize/2;
  int size2=nboxsize*mboxsize;
  //create a circle at (i0,j0)
  for (j=0;j<nboxsize;j++)
  {
    jn=j*mboxsize;
    for (i=0;i<mboxsize;i++)
    {
      if ((i-i0)*(i-i0)+(j-j0)*(j-j0)<domeR*domeR)
      {
        eta[i+jn]=1;
        eta[i+jn+1*size2]=0;
        eta[i+jn+2*size2]=0;
        eta[i+jn+3*size2]=0;
        eta[i+jn+4*size2]=0;
       }
    }
  }
  return eta;
}

