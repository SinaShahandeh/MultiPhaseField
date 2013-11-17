void drawphi(*phi,tn)
{
  // make picture of phi
  if  (tn % 100 ==0)
  {
    char filename[100];
    n=sprintf (filename, "im_ %d .png", tn);
    pngwriter png(mboxsize,nboxsize,1,filename); 
    for (i=0;i<size-1;i++)
    {
      for (j=0;j<size-1;j++)
      {
        png.plot(i,j, phi[i][j], phi[i][j], phi[i][j]);
      }
    }
    png.close()
  }
}