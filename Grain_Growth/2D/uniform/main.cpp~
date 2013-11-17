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
  int mboxsize=800*scale; // x axis in pixels
  int nboxsize=800*scale; // y axis
  double delx=2/scale;      // length unit per pixel
  
  double *eta;
  eta= new double[mboxsize*nboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*p];
  double *phi;
  phi= new double[mboxsize*nboxsize];
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
    eta[ii+jj*mboxsize+int(p*prand/RAND_MAX)*mboxsize*nboxsize]=1;
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
  int size=mboxsize*nboxsize;
  int jn, pnn;
  double del2[p];
  for (tn=1;tn<timesteps;tn++)
  {
    for (i=0;i<mboxsize;i++)
    {
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        // here is the sum of all order parameters^2 for the point i and j
        sumterm=0;
        for (psn=0;psn<p;psn++)
        {
          sumterm=sumterm+eta[i+jn+psn*size]*eta[i+jn+psn*size];
        }
        // calculation of nabla square eta
        for (pn=0;pn<p;pn++)
        {
          pnn=pn*size;
          del2[pn]=delx2*(0.5*(eta[peribc(i+1,nboxsize)+jn+pnn]-2*eta[i+jn+pnn]+eta[peribc(i-1,nboxsize)+jn+pnn])
          +0.25*(eta[peribc(i+2,nboxsize)+jn+pnn]-2*eta[i+jn+pnn]+eta[peribc(i-2,nboxsize)+jn+pnn]))
          +delx2*(0.5*(eta[i+peribc(j+1,mboxsize)*mboxsize+pnn]-2*eta[i+jn+pnn]+eta[i+peribc(j-1,mboxsize)*mboxsize+pnn])
          +0.25*(eta[i+peribc(j+2,mboxsize)*mboxsize+pnn]-2*eta[i+jn+pnn]+eta[i+peribc(j-2, mboxsize)*mboxsize+pnn]));
          sumtermp=eta[i+jn+pnn]*sumterm-pow(eta[i+jn+pnn],3);
          detadtM[pn]=-alpha*eta[i+jn+pnn]+beta*pow(eta[i+jn+pnn],3)-kappa*del2[pn];
          detadt[pn]=-L*(detadtM[pn]+2*gamma*sumtermp);
          eta2[i+jn+pnn]=eta[i+jn+pnn]+delt*detadt[pn];
          // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
          if (eta2[i+jn+pnn]>1) eta2[i+jn+pnn]=1;
          if (eta2[i+jn+pnn]<0) eta2[i+jn+pnn]=0;
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    for (i=0;i<mboxsize;i++)
    {
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (pind=0;pind<p;pind++)
        {
          pnn=pind*size;
          eta[i+jn+pnn]=eta2[i+jn+pnn];
        }
      }
    }
    // making the phi array
    for (i=0;i<mboxsize;i++)
    {
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        phi[i+jn]=0;
        for (pind=0;pind<p;pind++)
        {
          pnn=pind*size;
          phi[i+jn]=phi[i+jn]+eta[i+jn+pnn]*eta[i+jn+pnn];
        }
      }
    }
    // write array into a file each 100 time steps
    if  (tn % 100 ==0)
    {
      ofstream myfile;
      // make a string like "result_5.txt"
      int n;
      char filename[200];
      n=sprintf (filename, "results/result_ %d .txt", tn);
      myfile.open (filename);
      for (i=0;i<mboxsize;i++)
      {
        for (j=0;j<nboxsize;j++)
        {
          myfile << phi[i+j*mboxsize] << "	"; 
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
