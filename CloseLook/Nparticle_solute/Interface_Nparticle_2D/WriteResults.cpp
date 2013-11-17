
int WriteResults(double* phi, int mboxsize, int nboxsize,int p,int tn)
{
      ofstream myfile;
      myfile.open (filename);
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            if (phi[i+jn+kn]<0.90)
            {
              myfile << i <<" "<< j    <<" " << k   <<" " << phi[i+jn+kn] << endl; 
            }
          }
        }
      }
      myfile.close();
  
return 0;
}