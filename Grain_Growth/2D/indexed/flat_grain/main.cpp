// <source lang=cpp>

// 2D simulation with indexed algorithm
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>


// #include <pngwriter.h>

#include "peribc.cpp"
#include "symbc.cpp"
#include "GrainsStat.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int WriteResults();
double GrainsStat();

double sign(double x){
  if (x>0){return 1;}
  if (x<0){return -1;}
  if (x==0){return 0;}
  return 0;
}

// ------- main routin --------
int* MGsize;
int main(int a, char** charinput)
{
  char* dirstr;
  dirstr=charinput[1];
  cout <<"Data will be saved in this Directory:---->" << dirstr <<endl;
  
  bool fullwrite=true;
  // model parameters
  double delt=0.05;
  int timesteps1=250;
  int timesteps=70000;
  int writingtimesteps=1000;
  int writingstat=200;
  double L=1;
  double m=2.0;
  double gamma=2*1.5*m;
  double kappa=2.0;
  double Pz;
  Pz=double(atoi(charinput[2]))/1000;
  cout << "The friction force is = " << Pz <<"  ----" <<endl;
  double Lf=L, pzi;
  int i,j,tn;
  int maxind;
      double maxeta;
  // geometry settings
  int R=60;
  int scale=2;
  int mboxsize=600*scale; // x axis in pixels
  int nboxsize=600*scale; // y axis
  double delx=2/scale;      // length unit per pixel
  int p=5; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  double thresh=0.00001; //threshold value for choosing active nodes
  double *eta;
  eta= new double[mboxsize*nboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*p];
  int *inds;
  inds= new int[mboxsize*nboxsize*p];
  bool* mbool;
  mbool= new bool[mboxsize*nboxsize];
  int n;
  for (n=0; n<mboxsize*nboxsize*p;n++)
  {
    eta[n]=0;
    inds[n]=0;
  }
  for (n=0; n<mboxsize*nboxsize; n++)
  {
    mbool[n]=true;
  }
  
  double phii;
  double* phi;
  if (fullwrite==true)
  {
    phi= new double[mboxsize*nboxsize];
  }
  int* Maxindfield;
  Maxindfield= new int[mboxsize*nboxsize];
  double meanG;
  // number of nucleas at the beginning of simulation
  int nuclein;
  nuclein=int(mboxsize*nboxsize/100); // ~0.5 percent of grid points are nuclei
  int nn,ii,jj,pp,pos;
  double irand,jrand,prand;
  double sumterm,currenteta;
  double detadtM;
  double detadt;
  int pn,psn,pind;
  double delx2=1/(delx*delx);
  int size2=mboxsize*nboxsize;
  int jn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1;
  double del2;
  double sumeta, mineta, sumeta2;
  int currentind, indc, minind, indscn;
  MGsize=new int[nuclein];
  
  for (nn=2;nn<nuclein+1;nn++){
    irand=rand();
    jrand=rand();
    prand=rand();
    ii=int((nboxsize*irand)/RAND_MAX);
    jj=int((mboxsize*jrand)/RAND_MAX);
    pp=int(p*prand/RAND_MAX)*mboxsize*nboxsize;
    eta[ii+jj*mboxsize+pp]=1;
    inds[ii+jj*mboxsize+pp]=nn;
    inds[symbc(ii+1,mboxsize)+jj*mboxsize+pp]=nn;
    inds[symbc(ii-1,mboxsize)+jj*mboxsize+pp]=nn;
    inds[ii+peribc(jj+1,nboxsize)*mboxsize+pp]=nn;
    inds[ii+peribc(jj-1,nboxsize)*mboxsize+pp]=nn;
  }
  // putting a flat grain in the left side of the domain

  for (j=0;j<nboxsize;j++)
  {
    jn=j*mboxsize;
    for (i=0;i<mboxsize/20;i++)
    {
//         for (pind=0;pind<p;pind++) //make everything zero again
//           {
//             pp=pind*mboxsize*nboxsize;
//             eta[i+j*mboxsize+pp]=0;
//             inds[i+j*mboxsize+pp]=0;
//             inds[symbc(i+1,mboxsize)+j*mboxsize+pp]=0;
//             inds[symbc(i-1,mboxsize)+j*mboxsize+pp]=0;
//             inds[i+peribc(j+1,nboxsize)*mboxsize+pp]=0;
//             inds[i+peribc(j-1,nboxsize)*mboxsize+pp]=0;
//           }
        eta[i+j*mboxsize+0]=1; //put the big grain in layer 0
        inds[i+j*mboxsize+0]=1;
        inds[symbc(i+1,mboxsize)+j*mboxsize+0]=1;
        inds[symbc(i-1,mboxsize)+j*mboxsize+0]=1;
        inds[i+peribc(j+1,nboxsize)*mboxsize+0]=1;
        inds[i+peribc(j-1,nboxsize)*mboxsize+0]=1;
 
    }
  }
 
  // particles distribution specification

  //dynamics
  clock_t time1;
  time1=clock();
  cout << "Initialization ended." <<endl;
  for (tn=1;tn<timesteps1;tn++)
  {
    #pragma omp parallel for
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      jnplus1=peribc(j+1,nboxsize)*mboxsize;
      jnminus1=peribc(j-1,nboxsize)*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
        inplus1=symbc(i+1,mboxsize);
        inminus1=symbc(i-1,mboxsize);
        // here is the sum of all order parameters^2 for the point i and j
        sumterm=0;
        for (psn=0;psn<p;psn++)
        {
          sumterm=sumterm+eta[i+jn+psn*size2]*eta[i+jn+psn*size2];
        }
        // calculation of nabla square eta
        for (pn=0;pn<p;pn++)
        {
          pnn=pn*size2;
          mbool[i+jn]=false; //firts we make is false and then if detadt>thresh then it becomes true later
          //
          currentind=inds[i+jn+pnn];
          currenteta=eta[i+jn+pnn];
          //searching for neighbors with the same index as currentind
          sumeta=0;
          for (indc=0;indc<p;indc++)
          {
            indscn=indc*size2;
            if(currentind==inds[inplus1+jn+indscn]){sumeta=sumeta+eta[inplus1+jn+indc*size2];}
            if(currentind==inds[inminus1+jn+indscn]){sumeta=sumeta+eta[inminus1+jn+indc*size2];}
            if(currentind==inds[i+jnplus1+indscn]){sumeta=sumeta+eta[i+jnplus1+indc*size2];}
            if(currentind==inds[i+jnminus1+indscn]){sumeta=sumeta+eta[i+jnminus1+indc*size2];}
          }
          sumeta2=0;
          for (indc=0;indc<p;indc++)
          {
            indscn=indc*size2;
            if(currentind==inds[inplus1+jnplus1+indscn]){sumeta2=sumeta2+eta[inplus1+jnplus1+indc*size2];}
            if(currentind==inds[inminus1+jnminus1+indscn]){sumeta2=sumeta2+eta[inminus1+jnminus1+indc*size2];}
            if(currentind==inds[inminus1+jnplus1+indscn]){sumeta2=sumeta2+eta[inminus1+jnplus1+indc*size2];}
            if(currentind==inds[inplus1+jnminus1+indscn]){sumeta2=sumeta2+eta[inplus1+jnminus1+indc*size2];}
          }
          del2=delx2*(sumeta+0.25*sumeta2-5*currenteta);
          detadtM=m*(currenteta*currenteta*currenteta-currenteta)-kappa*del2;
          detadt=-L*(detadtM+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta));
          eta2[i+jn+pnn]=currenteta+delt*detadt;
          // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
          if (eta2[i+jn+pnn]>1) eta2[i+jn+pnn]=1;
          if (eta2[i+jn+pnn]<0) eta2[i+jn+pnn]=0;
          if (eta2[i+jn+pnn]>thresh) // make a padding for the changing order parameter so next time step it will take that into account.
          {
            mineta=eta[inplus1+jn];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[inplus1+jn+indc*size2]){mineta=eta[inplus1+jn+indc*size2]; minind=indc;}
              if(inds[inplus1+jn+indc*size2]==currentind){minind=indc; break;}
            }
            inds[inplus1+jn+minind*size2]=currentind;
            mineta=eta[inminus1+jn];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[inminus1+jn+indc*size2]){mineta=eta[inminus1+jn+indc*size2]; minind=indc;}
              if(inds[inminus1+jn+indc*size2]==currentind){minind=indc; break;}
            }
            inds[inminus1+jn+minind*size2]=currentind;
            mineta=eta[i+jnplus1];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[i+jnplus1+indc*size2]){mineta=eta[i+jnplus1+indc*size2]; minind=indc;}
              if(inds[i+jnplus1+indc*size2]==currentind){minind=indc; break;}
            }
            inds[i+jnplus1+minind*size2]=currentind;
            mineta=eta[i+jnminus1];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[i+jnminus1+indc*size2]){mineta=eta[i+jnminus1+indc*size2]; minind=indc;}
              if(inds[i+jnminus1+indc*size2]==currentind){minind=indc; break;}
            }
            inds[i+jnminus1+minind*size2]=currentind;
            // ------------------- second neighbors --------------------------
            // inminus1+jnminus1 -1 -1
            pos=inminus1+jnminus1;
            mineta=eta[pos];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
              if(inds[pos+indc*size2]==currentind){minind=indc; break;}
            }
            inds[pos+minind*size2]=currentind;
            // inplus1+jnminus1 1 -1
            pos=inplus1+jnminus1;
            mineta=eta[pos];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
              if(inds[pos+indc*size2]==currentind){minind=indc; break;}
            }
            inds[pos+minind*size2]=currentind;
            pos=inminus1+jnplus1; // -1 1
            mineta=eta[pos];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
              if(inds[pos+indc*size2]==currentind){minind=indc; break;}
            }
            inds[pos+minind*size2]=currentind;
            pos=inplus1+jnplus1; // 1 1
            mineta=eta[pos];  //the first in the table is assumed as smallest 
            minind=0;
            for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
            {
              if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
              if(inds[pos+indc*size2]==currentind){minind=indc; break;}
            }
            inds[pos+minind*size2]=currentind;
          }
        }
        phii=0;
        for (pind=0;pind<p;pind++)
        {
          pnn=pind*size2;
          phii=phii+eta[i+jn+pnn]*eta[i+jn+pnn];
        }
        if (phii<0.999999) //this point in at interface
        {
          mbool[i+jn]=true;
          mbool[inplus1+jn]=true;
          mbool[inminus1+jn]=true;
          mbool[i+jnplus1]=true;
          mbool[i+jnminus1]=true;
          mbool[inplus1+jnplus1]=true;
          mbool[inplus1+jnminus1]=true;
          mbool[inminus1+jnplus1]=true;
          mbool[inminus1+jnminus1]=true;
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    #pragma omp parallel for
    for (pind=0;pind<p;pind++)
    {
      pnn=pind*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          eta[i+jn+pnn]=eta2[i+jn+pnn];
        }
      }
    }
    
    if (tn%100==0)
      cout << "time required for 'timesteps1' time steps:" << double((clock()-time1))/double(CLOCKS_PER_SEC) << "seconds. \n";
    
  }
  cout <<tn <<endl;
  
  time1=clock();
  
  //optimized loop -----------------------------------------------------------------------------------------------------
  double pastsize, newsize;
  pastsize=100;
  newsize=10;
  tn=0;
  while (newsize<0.98*mboxsize*nboxsize) //newsize-pastsize>0.0000000000001
  {
    tn=tn+1;
    #pragma omp parallel for
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      jnplus1=peribc(j+1,mboxsize)*mboxsize;
      jnminus1=peribc(j-1,mboxsize)*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
        if (mbool[i+jn]==true)//mbool[i+jn+pnn]==true
        {
          inplus1=symbc(i+1,mboxsize);
          inminus1=symbc(i-1,mboxsize);
          mbool[i+jn]=false; //firts we make is false and then if detadt>thresh then it becomes true later
          for (pn=0;pn<p;pn++)
          {
            pnn=pn*size2;
            // here is the sum of all order parameters^2 for the point i and j
            sumterm=0;
            for (psn=0;psn<p;psn++)
            {
              sumterm=sumterm+eta[i+jn+psn*size2]*eta[i+jn+psn*size2];
            }
            currentind=inds[i+jn+pnn];
            currenteta=eta[i+jn+pnn];
            //searching for neighbors with the same index as currentind
            sumeta=0;
            for (indc=0;indc<p;indc++)
            {
              indscn=indc*size2;
              if(currentind==inds[inplus1+jn+indscn]){sumeta=sumeta+eta[inplus1+jn+indc*size2];}
              if(currentind==inds[inminus1+jn+indscn]){sumeta=sumeta+eta[inminus1+jn+indc*size2];}
              if(currentind==inds[i+jnplus1+indscn]){sumeta=sumeta+eta[i+jnplus1+indc*size2];}
              if(currentind==inds[i+jnminus1+indscn]){sumeta=sumeta+eta[i+jnminus1+indc*size2];}
            }
            sumeta2=0;
            for (indc=0;indc<p;indc++)
            {
              indscn=indc*size2;
              if(currentind==inds[inplus1+jnplus1+indscn]){sumeta2=sumeta2+eta[inplus1+jnplus1+indc*size2];}
              if(currentind==inds[inminus1+jnminus1+indscn]){sumeta2=sumeta2+eta[inminus1+jnminus1+indc*size2];}
              if(currentind==inds[inminus1+jnplus1+indscn]){sumeta2=sumeta2+eta[inminus1+jnplus1+indc*size2];}
              if(currentind==inds[inplus1+jnminus1+indscn]){sumeta2=sumeta2+eta[inplus1+jnminus1+indc*size2];}
            }
            del2=delx2*(sumeta+0.25*sumeta2-5*currenteta);
            detadtM=m*(currenteta*currenteta*currenteta-currenteta)+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta)-kappa*del2;
            pzi=3*currenteta*(1-currenteta)*sign(detadtM)*Pz;
            if (fabs(detadtM)<fabs(pzi)){
              Lf=0;
            }
            else{
              Lf=L;
            }
            detadt=-L*(detadtM-pzi);
            eta2[i+jn+pnn]=currenteta+delt*detadt;
            // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
            if (eta2[i+jn+pnn]>1) eta2[i+jn+pnn]=1;
         if (eta2[i+jn+pnn]<0) eta2[i+jn+pnn]=0;
         // creating phi to find interface area from it
         if (eta2[i+jn+pnn]>thresh) // make a padding for the changing order parameter so next time step it will take that into account.
            {
              mineta=eta[inplus1+jn];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[inplus1+jn+indc*size2]){mineta=eta[inplus1+jn+indc*size2]; minind=indc;}
                if(inds[inplus1+jn+indc*size2]==currentind){minind=indc; break;}
              }
              inds[inplus1+jn+minind*size2]=currentind;
              
              mineta=eta[inminus1+jn];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[inminus1+jn+indc*size2]){mineta=eta[inminus1+jn+indc*size2]; minind=indc;}
                if(inds[inminus1+jn+indc*size2]==currentind){minind=indc; break;}
              }
              inds[inminus1+jn+minind*size2]=currentind;
              
              mineta=eta[i+jnplus1];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[i+jnplus1+indc*size2]){mineta=eta[i+jnplus1+indc*size2]; minind=indc;}
                if(inds[i+jnplus1+indc*size2]==currentind){minind=indc; break;}
              }
              inds[i+jnplus1+minind*size2]=currentind;
              
              mineta=eta[i+jnminus1];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[i+jnminus1+indc*size2]){mineta=eta[i+jnminus1+indc*size2]; minind=indc;}
                if(inds[i+jnminus1+indc*size2]==currentind){minind=indc; break;}
              }
              inds[i+jnminus1+minind*size2]=currentind;
              
              // ------------------- second neighbors --------------------------
              // inminus1+jnminus1 -1 -1
              pos=inminus1+jnminus1;
              mineta=eta[pos];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
                if(inds[pos+indc*size2]==currentind){minind=indc; break;}
              }
              inds[pos+minind*size2]=currentind;
              
              // inplus1+jnminus1 1 -1
              pos=inplus1+jnminus1;
              mineta=eta[pos];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
                if(inds[pos+indc*size2]==currentind){minind=indc; break;}
              }
              inds[pos+minind*size2]=currentind;
              
              pos=inminus1+jnplus1; // -1 1
              mineta=eta[pos];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
                if(inds[pos+indc*size2]==currentind){minind=indc; break;}
              }
              inds[pos+minind*size2]=currentind;
              
              pos=inplus1+jnplus1; // 1 1
              mineta=eta[pos];  //the first in the table is assumed as smallest 
              minind=0;
              for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
                if(mineta>eta[pos+indc*size2]){mineta=eta[pos+indc*size2]; minind=indc;}
                if(inds[pos+indc*size2]==currentind){minind=indc; break;}
              }
              inds[pos+minind*size2]=currentind;
            }
          }
        }
        phii=0;
        for (pind=0;pind<p;pind++)
        {
          pnn=pind*size2;
          phii=phii+eta[i+jn+pnn]*eta[i+jn+pnn];
        }
        if (phii<0.9999999) //this point in at interface
        {
          mbool[i+jn]=true;
          mbool[inplus1+jn]=true;
          mbool[inminus1+jn]=true;
          mbool[i+jnplus1]=true;
          mbool[i+jnminus1]=true;
          mbool[inplus1+jnplus1]=true;
          mbool[inplus1+jnminus1]=true;
          mbool[inminus1+jnplus1]=true;
          mbool[inminus1+jnminus1]=true;
        }
      }
    }
    //setting eta equal to the new eta2 for the next time step
    #pragma omp parallel for
    for (pind=0;pind<p;pind++)
    {
      pnn=pind*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          eta[i+jn+pnn]=eta2[i+jn+pnn];
        }
      }
    }
    // write array into a file each 100 time steps
    if  (tn % writingtimesteps ==0)
    {
      cout << "Time required for time step number [" << tn << " ] is :" << double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds. \n";
      time1=clock();
      if (fullwrite==true)
      {
        // making the phi array
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            phi[i+jn]=0;
            for (pind=0;pind<p;pind++)
            {
              pnn=pind*size2;
              phi[i+jn]=phi[i+jn]+eta[i+jn+pnn]*eta[i+jn+pnn];
            }
          }
        }
        char filename[200];
        ofstream myfile2;
        // make a string like "result_5.txt"
        sprintf (filename, "%sFullres_%d.txt",dirstr, tn);
        myfile2.open (filename);
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            myfile2 << phi[i+jn] << "   "; 
          }
          myfile2 << endl;
        }
        myfile2.close();
        //writing order parameters
        /*
        char fileetaname[200];
        for (pind=0;pind<p;pind++)
        {
          pnn=pind*size2;
          ofstream fileeta;
          n=sprintf (fileetaname, "%sEta_%d_%d.txt",dirstr, pind, tn);
          fileeta.open (fileetaname);
          for (j=0;j<nboxsize;j++)
          {
            jn=j*mboxsize;
            for (i=0;i<mboxsize;i++)
            {
              fileeta << eta[i+jn+pnn] << "   "; 
            }
            fileeta << endl;
          }
          fileeta.close();
        }
        //writing indices
        char fileindsname[200];
        for (pind=0;pind<p;pind++)
        {
          pnn=pind*size2;
          ofstream fileinds;
          n=sprintf (fileindsname, "%sinds_%d_%d.txt",dirstr, pind, tn);
          fileinds.open (fileindsname);
          for (j=0;j<nboxsize;j++)
          {
            jn=j*mboxsize;
            for (i=0;i<mboxsize;i++)
            {
              fileinds << inds[i+jn+pnn] << "   "; 
            }
            fileinds << endl;
          }
          fileinds.close();
        }
        */
      }
      // writing indexed matrix
      char filename3[200];
      ofstream myfile3;
      sprintf (filename3, "%sInds_%d.txt",dirstr, tn);
      myfile3.open (filename3);
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          maxeta=eta[i+jn];
          maxind=0;
          for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
          {
            if(maxeta<eta[i+jn+indc*size2]){maxeta=eta[i+jn+indc*size2]; maxind=indc;}
          }
          Maxindfield[i+jn]=inds[i+jn+maxind*size2];
          myfile3 << Maxindfield[i+jn] << "   ";
        }
        myfile3 << endl;
      }
      myfile3.close();
    }
    if (tn % writingstat ==0){
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          maxeta=eta[i+jn];
          maxind=0;
          for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
          {
            if(maxeta<eta[i+jn+indc*size2]){maxeta=eta[i+jn+indc*size2]; maxind=indc;}
          }
          Maxindfield[i+jn]=inds[i+jn+maxind*size2];
        }
      }
     meanG=GrainsStat(Maxindfield, MGsize, nuclein, mboxsize, nboxsize);
      newsize=MGsize[1];
      cout << "Big grain size is : " << newsize*delx*delx <<endl;
      //writing grain statistics data
      char filenamestat[200];
      ofstream fileGstat;
      sprintf (filenamestat, "%sGrainStat_%d.txt",dirstr, tn);
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
