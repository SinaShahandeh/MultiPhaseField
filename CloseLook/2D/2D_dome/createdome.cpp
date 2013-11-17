//create a dome of value 1 in the order paramter number 1.
double* createdome(double *eta,int mboxsize, int nboxsize, int domeR)
{
  int i,j,jn;
  int i0=mboxsize/2;   //centre of the dome position
  int j0=nboxsize/2-nboxsize/8;
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
        eta[i+jn+size2]=0;
      }
      else
      {
        eta[i+jn+size2]=1;
        eta[i+jn]=0;
      }
    }
  }
  //fill under the half circle
  for (j=j0;j<nboxsize;j++)
  {
    jn=j*mboxsize;
    for(i=i0-domeR;i<i0+domeR;i++)
    {
      eta[i+jn]=1;
      eta[i+jn+size2]=0;
    }
  }
  return eta;
}

