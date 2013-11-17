// <source lang=cpp>

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
#include "GrainsStat.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int WriteResults();
int GrainsStat();
// ------- main routin --------
int* MGsize;
int main(int a, char** charinput)
{
  char* dirstr;
  dirstr=charinput[1];
  cout <<"Data will be saved in this Directory:---->" << dirstr <<endl;

  bool fullwrite=false;
  // model parameters
  double delt=0.1;
  int timesteps1=99;
  int timesteps=2001;
  int writingtimesteps=100;
  double L=1;
  double m=2.0;
  double gamma=2*1.5*m;
  double kappa=2.0;
  
  int i,j,k,tn;
  
  // geometry settings
  int scale=1;
  int mboxsize=50*scale; // x axis in pixels
  int nboxsize=50*scale; // y axis
  int kboxsize=50*scale;
  double delx=2/scale;      // length unit per pixel
  int p=6; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  double thresh=0.00001; //threshold value for choosing active nodes
  double *eta;
  eta= new double[mboxsize*nboxsize*kboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*kboxsize*p];
  int *inds;
  inds= new int[mboxsize*nboxsize*kboxsize*p];
  bool* mbool;
  mbool= new bool[mboxsize*nboxsize*kboxsize*p];
  double* phi;
  if (fullwrite==true)
  {
    phi= new double[mboxsize*nboxsize*kboxsize];
  }
  int* Maxindfield;
  Maxindfield= new int[mboxsize*nboxsize*kboxsize];
  // number of nucleas at the beginning of simulation
  int nuclein;
  nuclein=int(mboxsize*nboxsize*kboxsize/1000); // ~1 percent of grid points are nuclei
  int nn,ii,jj,kk,pp;
  double irand,jrand,prand, krand;
  
  MGsize=new int[nuclein];
  
  for (nn=1;nn<nuclein+1;nn++){
    irand=rand();
    jrand=rand();
    krand=rand();
    prand=rand();
    ii=int((nboxsize*irand)/RAND_MAX);
    jj=int((mboxsize*jrand)/RAND_MAX);
    kk=int((kboxsize*krand)/RAND_MAX);
    pp=int(p*prand/RAND_MAX)*mboxsize*nboxsize*kboxsize;
    eta[ii+jj*mboxsize+kk*mboxsize*nboxsize+pp]=1;
    inds[ii+jj*mboxsize+kk*mboxsize*nboxsize+pp]=nn;
    inds[peribc(ii+1,mboxsize)+jj*mboxsize+kk*mboxsize*nboxsize+pp]=nn;
    inds[peribc(ii-1,mboxsize)+jj*mboxsize+kk*mboxsize*nboxsize+pp]=nn;
    inds[ii+peribc(jj+1,nboxsize)*mboxsize+kk*mboxsize*nboxsize+pp]=nn;
    inds[ii+peribc(jj-1,nboxsize)*mboxsize+kk*mboxsize*nboxsize+pp]=nn;
    inds[ii+jj*mboxsize+peribc(kk+1,kboxsize)*mboxsize*nboxsize+pp]=nn;
    inds[ii+jj*mboxsize+peribc(kk-1,kboxsize)*mboxsize*nboxsize+pp]=nn;
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
  double sumterm,currenteta;
  double detadtM;
  double detadt;
  int pn,psn,pind;
  double delx2=1/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double del2;
  double sumeta, mineta;
  int currentind, indc, minind;
  //calculating processing time
  clock_t time1;
  time1=clock();
  cout << "Initialization ended." <<endl;
  for (tn=1;tn<timesteps1;tn++)
  {
    #pragma omp parallel for
    for (k=0;k<kboxsize;k++)
    {
      kn=k*size2;
      knplus1=peribc(k+1,kboxsize)*size2;
      knminus1=peribc(k-1,kboxsize)*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        jnplus1=peribc(j+1,nboxsize)*mboxsize;
        jnminus1=peribc(j-1,nboxsize)*mboxsize;
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
            mbool[i+jn+kn+pnn]=false; //firts we make is false and then if detadt>thresh then it becomes true later
            //
            currentind=inds[i+jn+kn+pnn];
            currenteta=eta[i+jn+kn+pnn];
            //searching for neighbors with the same index as currentind
            sumeta=0;
            for (indc=0;indc<p;indc++)
            {
              if(currentind==inds[inplus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inplus1+jn+kn+indc*size3];}
            }
            for (indc=0;indc<p;indc++)
            {
              if(currentind==inds[inminus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inminus1+jn+kn+indc*size3];}
            }
            for (indc=0;indc<p;indc++)
            {
              if(currentind==inds[i+jnplus1+kn+(indc)*size3]){sumeta=sumeta+eta[i+jnplus1+kn+indc*size3];}
            }
            for (indc=0;indc<p;indc++)
            {
              if(currentind==inds[i+jnminus1+kn+indc*size3]){sumeta=sumeta+eta[i+jnminus1+kn+indc*size3];}
            }
            for (indc=0;indc<p;indc++)
            {
              if(currentind==inds[i+jn+knplus1+indc*size3]){sumeta=sumeta+eta[i+jn+knplus1+indc*size3];}
            }
            for (indc=0;indc<p;indc++)
            {
              if(currentind==inds[i+jn+knminus1+indc*size3]){sumeta=sumeta+eta[i+jn+knminus1+indc*size3];}
            }
            del2=delx2*(sumeta-6*currenteta);
            detadtM=m*(currenteta*currenteta*currenteta-currenteta)-kappa*del2;
            detadt=-L*(detadtM+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta));
            if (fabs(detadt)>thresh) // make a padding for the changing order parameter so next time step it will take that into account.
            {
              mbool[i+jn+kn+pnn]=true;
              mineta=eta[inplus1+jn+kn];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[inplus1+jn+kn+indc*size3]){mineta=eta[inplus1+jn+kn+indc*size3]; minind=indc;}
                if(inds[inplus1+jn+kn+indc*size3]==currentind){minind=indc; break;}
              }
              inds[inplus1+jn+kn+minind*size3]=currentind;
              mbool[inplus1+jn+kn+minind*size3]=true;
              mineta=eta[inminus1+jn+kn];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[inminus1+jn+kn+indc*size3]){mineta=eta[inminus1+jn+kn+indc*size3]; minind=indc;}
                if(inds[inminus1+jn+kn+indc*size3]==currentind){minind=indc; break;}
              }
              inds[inminus1+jn+kn+minind*size3]=currentind;
              mbool[inminus1+jn+kn+minind*size3]=true;
              mineta=eta[i+jnplus1+kn];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[i+jnplus1+kn+indc*size3]){mineta=eta[i+jnplus1+kn+indc*size3]; minind=indc;}
                if(inds[i+jnplus1+kn+indc*size3]==currentind){minind=indc; break;}
              }
              inds[i+jnplus1+kn+minind*size3]=currentind;
              mbool[i+jnplus1+kn+minind*size3]=true;
              mineta=eta[i+jnminus1+kn];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[i+jnminus1+kn+indc*size3]){mineta=eta[i+jnminus1+kn+indc*size3]; minind=indc;}
                if(inds[i+jnminus1+kn+indc*size3]==currentind){minind=indc; break;}
              }
              inds[i+jnminus1+kn+minind*size3]=currentind;
              mbool[i+jnminus1+kn+minind*size3]=true;
              mineta=eta[i+jn+knplus1];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[i+jn+knplus1+indc*size3]){mineta=eta[i+jn+knplus1+indc*size3]; minind=indc;}
                if(inds[i+jn+knplus1+indc*size3]==currentind){minind=indc; break;}
              }
              inds[i+jn+knplus1+minind*size3]=currentind;
              mbool[i+jn+knplus1+minind*size3]=true;
              mineta=eta[i+jn+knminus1];  //the first i the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[i+jn+knminus1+indc*size3]){mineta=eta[i+jn+knminus1+indc*size3]; minind=indc;}
                if(inds[i+jn+knminus1+indc*size3]==currentind){minind=indc; break;}
              }
              inds[i+jn+knminus1+minind*size3]=currentind;
              mbool[i+jn+knminus1+minind*size3]=true;
            }
            eta2[i+jn+kn+pnn]=currenteta+delt*detadt;
            // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
            if (eta2[i+jn+kn+pnn]>1) eta2[i+jn+kn+pnn]=1;
            if (eta2[i+jn+kn+pnn]<0) eta2[i+jn+kn+pnn]=0;
          }
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    #pragma omp parallel for
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
    cout << tn << " \n";
  }
  cout << endl << "time required for 'timesteps1' time steps:" << double((clock()-time1))/double(CLOCKS_PER_SEC) << "seconds. \n";
  time1=clock();
  
  //optimized loop -----------------------------------------------------------------------------------------------------
  for (tn=timesteps1;tn<timesteps;tn++)
  {
    #pragma omp parallel for
    for (k=0;k<kboxsize;k++)
    {
      kn=k*size2;
      knplus1=peribc(k+1,kboxsize)*size2;
      knminus1=peribc(k-1,kboxsize)*size2;
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
              mbool[i+jn+kn+pnn]=false; //firts we make is false and then if detadt>thresh then it becomes true later
              inplus1=peribc(i+1,mboxsize);
              inminus1=peribc(i-1,mboxsize);
              // here is the sum of all order parameters^2 for the point i and j
              sumterm=0;
              for (psn=0;psn<p;psn++)
              {
                sumterm=sumterm+eta[i+jn+kn+psn*size3]*eta[i+jn+kn+psn*size3];
              }
              // calculation of nabla square eta
              currentind=inds[i+jn+kn+pnn];
              currenteta=eta[i+jn+kn+pnn];
              //searching for neighbors with the same index as currentind
              sumeta=0;
              for (indc=0;indc<p;indc++)
              {
                if(currentind==inds[inplus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inplus1+jn+kn+indc*size3];}
              }
              for (indc=0;indc<p;indc++)
              {
                if(currentind==inds[inminus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inminus1+jn+kn+indc*size3];}
              }
              for (indc=0;indc<p;indc++)
              {
                if(currentind==inds[i+jnplus1+kn+(indc)*size3]){sumeta=sumeta+eta[i+jnplus1+kn+indc*size3];}
              }
              for (indc=0;indc<p;indc++)
              {
                if(currentind==inds[i+jnminus1+kn+indc*size3]){sumeta=sumeta+eta[i+jnminus1+kn+indc*size3];}
              }
              for (indc=0;indc<p;indc++)
              {
                if(currentind==inds[i+jn+knplus1+indc*size3]){sumeta=sumeta+eta[i+jn+knplus1+indc*size3];}
              }
              for (indc=0;indc<p;indc++)
              {
                if(currentind==inds[i+jn+knminus1+indc*size3]){sumeta=sumeta+eta[i+jn+knminus1+indc*size3];}
              }
              del2=delx2*(sumeta-6*currenteta);
              detadtM=m*(currenteta*currenteta*currenteta-currenteta)-kappa*del2;
              detadt=-L*(detadtM+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta));
              if (fabs(detadt)>thresh) // make a padding for the changing order parameter so next time step it will take that into account.
              {
                mbool[i+jn+kn+pnn]=true;
                mineta=eta[inplus1+jn+kn];  //the first in the table is assumed as smallest 
                minind=0;
                for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                {
                  if(mineta>eta[inplus1+jn+kn+indc*size3]){mineta=eta[inplus1+jn+kn+indc*size3]; minind=indc;}
                  if(inds[inplus1+jn+kn+indc*size3]==currentind){minind=indc; break;}
                }
                inds[inplus1+jn+kn+minind*size3]=currentind;
                mbool[inplus1+jn+kn+minind*size3]=true;
                mineta=eta[inminus1+jn+kn];  //the first in the table is assumed as smallest 
                minind=0;
                for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                {
                  if(mineta>eta[inminus1+jn+kn+indc*size3]){mineta=eta[inminus1+jn+kn+indc*size3]; minind=indc;}
                  if(inds[inminus1+jn+kn+indc*size3]==currentind){minind=indc; break;}
                }
                inds[inminus1+jn+kn+minind*size3]=currentind;
                mbool[inminus1+jn+kn+minind*size3]=true;
                mineta=eta[i+jnplus1+kn];  //the first in the table is assumed as smallest 
                minind=0;
                for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                {
                  if(mineta>eta[i+jnplus1+kn+indc*size3]){mineta=eta[i+jnplus1+kn+indc*size3]; minind=indc;}
                  if(inds[i+jnplus1+kn+indc*size3]==currentind){minind=indc; break;}
                }
                inds[i+jnplus1+kn+minind*size3]=currentind;
                mbool[i+jnplus1+kn+minind*size3]=true;
                mineta=eta[i+jnminus1+kn];  //the first in the table is assumed as smallest 
                minind=0;
                for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                {
                  if(mineta>eta[i+jnminus1+kn+indc*size3]){mineta=eta[i+jnminus1+kn+indc*size3]; minind=indc;}
                  if(inds[i+jnminus1+kn+indc*size3]==currentind){minind=indc; break;}
                }
                inds[i+jnminus1+kn+minind*size3]=currentind;
                mbool[i+jnminus1+kn+minind*size3]=true;
                mineta=eta[i+jn+knplus1];  //the first in the table is assumed as smallest 
                minind=0;
                for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                {
                  if(mineta>eta[i+jn+knplus1+indc*size3]){mineta=eta[i+jn+knplus1+indc*size3]; minind=indc;}
                  if(inds[i+jn+knplus1+indc*size3]==currentind){minind=indc; break;}
                }
                inds[i+jn+knplus1+minind*size3]=currentind;
                mbool[i+jn+knplus1+minind*size3]=true;
                mineta=eta[i+jn+knminus1];  //the first i the table is assumed as smallest 
                minind=0;
                for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                {
                  if(mineta>eta[i+jn+knminus1+indc*size3]){mineta=eta[i+jn+knminus1+indc*size3]; minind=indc;}
                  if(inds[i+jn+knminus1+indc*size3]==currentind){minind=indc; break;}
                }
                inds[i+jn+knminus1+minind*size3]=currentind;
                mbool[i+jn+knminus1+minind*size3]=true;
              }
              eta2[i+jn+kn+pnn]=currenteta+delt*detadt;
              // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
              if (eta2[i+jn+kn+pnn]>1) eta2[i+jn+kn+pnn]=1;
              if (eta2[i+jn+kn+pnn]<0) eta2[i+jn+kn+pnn]=0;
            }
          }
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    #pragma omp parallel for
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
    // write array into a file each 100 time steps
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
  }
return 0;
}


// </source>
