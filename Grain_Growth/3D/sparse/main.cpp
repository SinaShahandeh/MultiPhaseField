// 3D simulation with a optimized calculations
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>


// #include <pngwriter.h>

#include "peribc.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int WriteResults();
// ------- main routin --------
int main()
{
  ofstream myfile ;
  // model parameters
  double delt=0.1;
  int timesteps1=11;
  int timesteps=1001;
  
  double L=1;
  double alpha=1.0 ;
  double beta=1.0 ;
  double gamma=1.0;
  double kappa=2.0;
  double epsilon=5.0;
  
  int i,j,k,tn;
  
  // geometry settings
  int p=10; // phase field numbers
  int scale=1;
  int mboxsize=50*scale; // x axis in pixels
  int nboxsize=50*scale; // y axis
  int kboxsize=50*scale;
  double delx=2/scale;      // length unit per pixel
  
  double thresh=0.0001; //threshold value for choosing active nodes
  double *eta;
  eta= new double[mboxsize*nboxsize*kboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*kboxsize*p];
  bool *mbool;
  mbool= new bool[mboxsize*nboxsize*kboxsize*p];
  // number of nucleas at the beginning of simulation
  double nuclein;
  nuclein=mboxsize*nboxsize*kboxsize/200; // ~5 percent of grid points are nuclei
  nuclein=3;
  //setting initial condition  (nuclei of grains)
    int initialpos=int(0.80*kboxsize);
    for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=0;i<mboxsize;i++)
      {
        if (i<initialpos)
        {
          eta[i+j*mboxsize+k*mboxsize*nboxsize]=1;
          eta[i+j*mboxsize+k*mboxsize*nboxsize+1*mboxsize*nboxsize*kboxsize]=0;
        }
        else
        {
          eta[i+j*mboxsize+k*mboxsize*nboxsize]=0;
          eta[i+j*mboxsize+k*mboxsize*nboxsize+1*mboxsize*nboxsize*kboxsize]=1;
        }
      }
    }
  }
/*
int nn,ii,jj,kk;
  double irand,jrand,prand, krand;
  for (nn=0;nn<nuclein;nn++){
    irand=rand();
    jrand=rand();
    krand=rand();
    prand=rand();
    ii=int((nboxsize*irand)/RAND_MAX);
    jj=int((mboxsize*jrand)/RAND_MAX);
    kk=int((kboxsize*krand)/RAND_MAX);
    eta[ii+jj*mboxsize+kk*mboxsize*nboxsize+int(p*prand/RAND_MAX)*mboxsize*nboxsize*kboxsize]=1;
  }
  // particles distribution specification
  /* double diameter=2;
  double particles_fraction=0.00;
  double particlesn=particles_fraction*mboxsize*nboxsize/diameter^2   //
  particles number
  double ppf[nboxsize][mboxsize]
  here goes the function to make particle distribution 
  */
  //dynamics
  double sumterm,sumtermp;
  double detadtM;
  double detadt;
  int pn,psn,pind;
  double delx2=1/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double del2;
  double test;
  //calculating processing time
  clock_t time1;
  time1=clock();
  cout << "Initialization ended." <<endl;
  for (tn=1;tn<timesteps1;tn++)
  {
    //  cout <<test << endl;
    test=0;
    for (k=0;k<kboxsize;k++)
    {
      kn=k*size2;
      knplus1=peribc(k+1,kboxsize)*size2;
      knminus1=peribc(k+1,kboxsize)*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        jnplus1=peribc(j+1,mboxsize)*mboxsize;
        jnminus1=peribc(j-1,mboxsize)*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          inplus1=peribc(i+1,mboxsize);
          inminus1=peribc(i-1,mboxsize);
          // here is the sum of all order parameters^2 for the point i and j
          sumterm=0;
          for (psn=0;psn<p;psn++)
          {
            sumterm=sumterm+eta[i+jn+kn+psn*size3]*eta[i+jn+kn+psn*size3];
          }
          // calculation of nabla square eta
          for (pn=0;pn<p;pn++)
          {
            pnn=pn*size3;
            del2=delx2*((eta[inplus1+jn+kn+pnn]+eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]+eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]+eta[i+jn+knminus1+pnn])-6*eta[i+jn+kn+pnn]);
            sumtermp=eta[i+jn+kn+pnn]*sumterm-pow(eta[i+jn+kn+pnn],3);
            detadtM=-alpha*eta[i+jn+kn+pnn]+beta*pow(eta[i+jn+kn+pnn],3)-kappa*del2;
            detadt=-L*(detadtM+2*gamma*sumtermp);
            if (fabs(detadt)>thresh) // optimization function
            {
              mbool[i+jn+kn+pnn]=true;
              mbool[inplus1+jn+kn+pnn]=true;
              mbool[inminus1+jn+kn+pnn]=true;
              mbool[i+jnplus1+kn+pnn]=true;
              mbool[i+jnminus1+kn+pnn]=true;
              mbool[i+jn+knplus1+pnn]=true;
              mbool[i+jn+knminus1+pnn]=true;
            }
            else
            {
              mbool[i+jn+kn+pnn]=false; 
            }
            test=test+eta[i+jn+kn+pnn];
            
            eta2[i+jn+kn+pnn]=eta[i+jn+kn+pnn]+delt*detadt;
            // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
            if (eta2[i+jn+kn+pnn]>1) eta2[i+jn+kn+pnn]=1;
            if (eta2[i+jn+kn+pnn]<0) eta2[i+jn+kn+pnn]=0;
          }
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    for (pind=0;pind<p;pind++)
    {
      pnn=pind*size3;
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            eta[i+jn+kn+pnn]=eta2[i+jn+kn+pnn];
          }
        }
      }
    }
    cout << tn << "\n";
  }
  cout << "time required for 10 time steps:" << double((clock()-time1))/double(CLOCKS_PER_SEC) << "seconds. \n";
  //optimized loop -----------------------------------------------------------------------------------------------------
  for (tn=1;tn<timesteps;tn++)
  {
    time1=clock();
    for (k=0;k<kboxsize;k++)
    {
      kn=k*size2;
      knplus1=peribc(k+1,kboxsize)*size2;
      knminus1=peribc(k+1,kboxsize)*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        jnplus1=peribc(j+1,mboxsize)*mboxsize;
        jnminus1=peribc(j-1,mboxsize)*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          for (pn=0;pn<p;pn++)
          {
            pnn=pn*size3;
            if (mbool[i+jn+kn+pnn]==true)
            {
              inplus1=peribc(i+1,mboxsize);
              inminus1=peribc(i-1,mboxsize);
              // here is the sum of all order parameters^2 for the point i and j
              sumterm=0;
              for (psn=0;psn<p;psn++)
              {
                sumterm=sumterm+eta[i+jn+kn+psn*size3]*eta[i+jn+kn+psn*size3];
              }
              // calculation of nabla square eta
              del2=delx2*((eta[inplus1+jn+kn+pnn]+eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]+eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]+eta[i+jn+knminus1+pnn])
              -6*eta[i+jn+kn+pnn]);
              sumtermp=eta[i+jn+kn+pnn]*sumterm-pow(eta[i+jn+kn+pnn],3);
              detadtM=-alpha*eta[i+jn+kn+pnn]+beta*pow(eta[i+jn+kn+pnn],3)-kappa*del2;
              detadt=-L*(detadtM+2*gamma*sumtermp);
              if (fabs(detadt)>thresh) // optimization function
              {
                mbool[i+jn+kn+pnn]=true;
                mbool[inplus1+jn+kn+pnn]=true;
                mbool[inminus1+jn+kn+pnn]=true;
                mbool[i+jnplus1+kn+pnn]=true;
                mbool[i+jnminus1+kn+pnn]=true;
                mbool[i+jn+knplus1+pnn]=true;
                mbool[i+jn+knminus1+pnn]=true;
              }
              else
              {
                mbool[i+jn+kn+pnn]=false; 
              }
              eta2[i+jn+kn+pnn]=eta[i+jn+kn+pnn]+delt*detadt;
              // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
              if (eta2[i+jn+kn+pnn]>1) eta2[i+jn+kn+pnn]=1;
              if (eta2[i+jn+kn+pnn]<0) eta2[i+jn+kn+pnn]=0;
            }
          }
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    for (pind=0;pind<p;pind++)
    {
      pnn=pind*size3;
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            eta[i+jn+kn+pnn]=eta2[i+jn+kn+pnn];
          }
        }
      }
    }
    cout << "Time required for time step number [" << tn << " ] is :" << double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds. \n";
    // write array into a file each 100 time steps
    if  (tn % 100 ==0)
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
      ofstream myfile;
      // make a string like "result_5.txt"
      int n;
      char filename[200];
      n=sprintf (filename, "results/result_ %d .txt", tn);
      myfile.open (filename);
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            if (phi[i+jn+kn]<0.99)
            {
              myfile << i <<"   "<< j    <<"   " << k   <<"    " << phi[i+jn+kn] << endl; 
            }
          }
        }
      }
      myfile.close();
      
      // writing
      ofstream myfile2;
      // make a string like "result_5.txt"
      n=sprintf (filename, "results/Fullres_ %d .txt", tn);
      myfile2.open (filename);
      for (k=0;k<1;k++)
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
        myfile2.close();
      }
      
    }
  }
  return 0;
}
