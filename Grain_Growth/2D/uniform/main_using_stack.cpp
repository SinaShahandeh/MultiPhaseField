//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>


// #include <pngwriter.h>

#include "peribc.cpp"

using namespace std;
int peribc();

// ------- main routin --------
int main()
{
  // model parameters
  double delt=0.2;
  int timesteps=1000;
  
  double L=1;
  double alpha=1 ;
  double beta=1 ;
  double gamma=1;
  double kappa=2;
  double epsilon=5;
  
  int i,j,tn;
  
  // geometry settings
  int p=25; // phase field numbers
  int scale=1;
  int mboxsize=100*scale; // x axis in pixels
  int nboxsize=100*scale; // y axis
  double delx=2/scale;      // length unit per pixel
  
  
  double eta[mboxsize][nboxsize][p];
  double eta2[mboxsize][nboxsize][p];
  double phi[mboxsize][nboxsize];
  // number of nucleas at the beginning of simulation
  double nuclein;
  nuclein=mboxsize*nboxsize/200; // ~5 percent of grid points are nuclei
  //setting initial condition  (nuclei of grains)
  int nn,ii,jj;
  double irand,jrand,prand;
  for (nn=0;nn<nuclein;nn++)
  {
    irand=rand();
    jrand=rand();
    prand=rand();
    ii=int((nboxsize*irand)/RAND_MAX);
    jj=int((mboxsize*jrand)/RAND_MAX);
    eta[ii][jj][int(p*prand/RAND_MAX)]=1;
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
  double detadtM[p];
  double detadt[p];
  int pn,psn,pind;
  double delx2=1.0/(delx*delx);
  double del2[p];
  for (tn=1;tn<timesteps;tn++)
  {
    for (i=0;i<mboxsize;i++)
    {
      for (j=0;j<nboxsize;j++)
      {
        // here is the sum of all order parameters^2 for the point i and j
        sumterm=0;
        for (psn=0;psn<p;psn++)
        {
          sumterm=sumterm+eta[i][j][psn]*eta[i][j][psn];
        }
        // calculation of nabla square eta
        for (pn=0;pn<p;pn++)
        {
          
          del2[pn]=delx2*(0.5*(eta[peribc(i+1,nboxsize)][j][pn]-2*eta[i][j][pn]+eta[peribc (i-1,nboxsize)][j][pn])+ 0.25*(eta[peribc(i+2,nboxsize)][j][pn]-2*eta[i][j][pn]+eta[peribc(i-2,nboxsize)] [j][pn]))
          +delx2*(0.5*(eta[i][peribc(j+1,mboxsize)][pn]-2*eta[i][j][pn]+eta[i][peribc(j-1, mboxsize)][pn])
          +0.25*(eta[i][peribc(j+2,mboxsize)][pn]-2*eta[i][j][pn]+eta[i][peribc(j-2, mboxsize)][pn]));
          sumtermp=eta[i][j][pn]*sumterm-pow(eta[i][j][pn],3);
          detadtM[pn]=-alpha*eta[i][j][pn]+beta*pow(eta[i][j][pn],3)-kappa*del2[pn];
          detadt[pn]=-L*(detadtM[pn]+2*gamma*sumtermp);
          eta2[i][j][pn]=eta[i][j][pn]+delt*detadt[pn];
          // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
          if (eta2[i][j][pn]>1) eta2[i][j][pn]=1;
          if (eta2[i][j][pn]<0) eta2[i][j][pn]=0;
        }
      }
    }
    //setting T equal to the new T2 for next time step
    for (i=0;i<mboxsize;i++)
    {
      for (j=0;j<nboxsize;j++)
      {
        for (pind=0;pind<p;pind++)
        {
          eta[i][j][pind]=eta2[i][j][pind];
        }
      }
    }
    // making the phi array
    for (i=0;i<mboxsize;i++)
    {
      for (j=0;j<nboxsize;j++)
      {
        phi[i][j]=0;
        for (pind=0;pind<p;pind++)
        {
          phi[i][j]=phi[i][j]+eta[i][j][pind]*eta[i][j][pind];
        }
      }
    }
    // write array into a file each 100 time steps
    if  (tn % 100 ==0)
    {
      ofstream myfile;
      // make a string like "result_5.txt"
      int n;
      char filename[100];
      n=sprintf (filename, "results/result_ %d .txt", tn);
      myfile.open (filename);
      for (i=0;i<mboxsize;i++)
      {
        for (j=0;j<nboxsize;j++)
        {
          myfile << phi[i][j] << "	"; 
        }
        myfile << "\n";
      }
      myfile.close();
    }
/*
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
*/
    cout << tn << "\n";
  }
  return 0;
}
