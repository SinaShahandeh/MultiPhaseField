 
//create a rectangle
int rectgrain(double *eta,int mboxsize, int nboxsize, int pi, int x1, int x2, int y1, int y2)
{
  int i,j,jn;
  int size2=nboxsize*mboxsize;
  //create a circle at (i0,j0)
  for (j=x1;j<x2;j++)
  {
    jn=j*mboxsize;
    for (i=y1;i<y2;i++)
    {
     eta[i+jn+pi*size2]=1;
    }
  }
  return 0;
}

