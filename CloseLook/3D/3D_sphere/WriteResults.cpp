void WriteResults(double* eta, double* phi, int
if  (tn % writingtimesteps ==0)
    {
      cout << "Time required for time step number [" << tn << " ] is :" << double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds. \n";
      time1=clock();
      int n;
      if (fullwrite==true)
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
        char filename[200];
        ofstream myfile2;
        // make a string like "result_5.txt"
        n=sprintf (filename, "%sFullres_%d.txt",dirstr, tn);
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
      }
      // writing indexed matrix
      int maxind;
      double maxeta;
      char filename3[200];
      ofstream myfile3;
      // make a string like "result_5.txt"
      n=sprintf (filename3, "%sInds_%d.txt",dirstr, tn);
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
            Maxindfield[i+jn+kn]=inds[i+jn+kn+maxind*size3];
            myfile3 << Maxindfield[i+jn+kn] << "   ";
            
          }
          myfile3 << endl;
        }
      }
      myfile3.close();
      int meanG;
      meanG=GrainsStat(Maxindfield, MGsize, nuclein, mboxsize, nboxsize, kboxsize);
      cout << "average grain size is : " << meanG <<endl;
      //writing grain statistics data
      char filenamestat[200];
      ofstream fileGstat;
      n=sprintf (filenamestat, "%sGrainStat_%d.txt",dirstr, tn);
      fileGstat.open (filenamestat);
      for (i=1;i<nuclein+1;i++)
      {
        fileGstat << MGsize[i] << endl;
      }
      fileGstat.close();
    }