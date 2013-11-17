
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

int WriteResults(double eta[], double Inds[], int mboxsize, int nboxsize, int kboxsize, int p,int tn)
{
    double *phi;
  phi= new double[mboxsize*nboxsize*kboxsize];
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
      // writing
      int n;
      char filename[200];
      // writing
      ofstream myfile2;
      // make a string like "result_5.txt"
      n=sprintf (filename, "/media/disk/sim_res/test/Fullres_%d.txt", tn);
      myfile2.open (filename);
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            myfile2 << phi[i+jn+kn] << "   "; 
          }
          myfile2 << endl;
        }
      }
      myfile2.close();
      // writing
      int maxind;
      double maxeta;
      char filename3[200];
      // writing
      ofstream myfile3;
      // make a string like "result_5.txt"
      n=sprintf (filename3, "/media/disk/sim_res/test/Inds_%d.txt", tn);
      myfile3.open (filename3);
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            maxeta=eta[i+jn+kn];
            maxind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(maxeta<eta[i+jn+kn+indc*size3]){maxeta=eta[i+jn+kn+indc*size3]; maxind=indc;}
            }
            myfile3 << inds[i+jn+kn+maxind*size3] << "   "; 
            
          }
          myfile3 << endl;
        }
      }
      myfile3.close();
      
  */
return(0);
}


// For making index maps:


   //--------------------------------------
   if  (tn % 5 ==0)
   {
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
     // writing
     int n;
     int maxind;
     double maxeta;
     char filename[200];
     // writing
     ofstream myfile2;
     // make a string like "result_5.txt"
     n=sprintf (filename, "results/Fullres_ %d .txt", tn);
     myfile2.open (filename);
     for (k=0;k<kboxsize;k++)
     {
       kn=k*size2;
       for (j=0;j<nboxsize;j++)
       {
         jn=j*mboxsize;
         for (i=0;i<mboxsize;i++)
         {
           maxeta=eta[i+jn+kn];
           maxind=0;
           for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
             {
               if(maxeta<eta[i+jn+kn+indc*size3]){maxeta=eta[i+jn+kn+indc*size3]; maxind=indc;}
             }
           myfile2 << inds[i+jn+kn+maxind*size3] << "   "; 
           
         }
         myfile2 << endl;
       }
     }
     myfile2.close();
   }
   
   