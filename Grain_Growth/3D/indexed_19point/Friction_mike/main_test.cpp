
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
//#include "resume.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int WriteResults();
int GrainsStat();
int resume();

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
  //sprintf(dirstr, "/media/disk/sim_res/test/");
  dirstr=charinput[1];
  cout <<"Data will be saved in this Directory:---->" << dirstr <<endl;

  bool fullwrite=true;
  // model parameters
  double delt=0.1;
  int timesteps1=2;
  int timesteps=2;
  int resumetimestep=43000;
  int writingtimesteps=1000;
  double L=1;
  double m=1.0;
  double gamma=2*1.5*m;
  double kappa=2.0;
  double Pz=0.03;
  double Lf=L;
  int i,j,k,tn;
  int n;

  // geometry settings
  int scale=2;
  int mboxsize=50*scale; // x axis in pixels
  int nboxsize=50*scale; // y axis
  int kboxsize=50*scale;
  double delx=2/scale;      // length unit per pixel
  int p=8; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  double thresh=0.0000001; //threshold value for choosing active nodes
  double *eta;
  eta= new double[mboxsize*nboxsize*kboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*kboxsize*p];
  int *inds;
  inds= new int[mboxsize*nboxsize*kboxsize*p];
  bool* mbool;
  mbool= new bool[mboxsize*nboxsize*kboxsize*p];
  double* phi;

    phi= new double[mboxsize*nboxsize*kboxsize];
  int* Maxindfield;
  Maxindfield= new int[mboxsize*nboxsize*kboxsize];
 for (i=0;i<mboxsize*nboxsize*kboxsize*p;i++)
  {
     eta[i]=0;eta2[i]=0;inds[i]=0;mbool[i]=true;
  }
 for (i=0;i<mboxsize*nboxsize*kboxsize;i++)
  {
     phi[i]=0;Maxindfield[i]=0;
  }
  // number of nucleas at the beginning of simulation
  int nuclein;
  nuclein=int(mboxsize*nboxsize*kboxsize/1000/scale/scale/scale); // ~1 percent of grid points are nuclei
  int nn,ii,jj,kk,pp;
  double irand,jrand,prand, krand;
  
  MGsize=new int[nuclein];
  cout << "step 1:before nucleation" <<endl;
  // resuming from a file
  //  if (resumetimestep>0)  resume( dirstr, inds, eta, resumetimestep, mboxsize, nboxsize, kboxsize);
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
  
  double sumterm,currenteta;
  double detadtM;
  double detadt;
  int pn,psn,pind,pos;
  double delx2=0.166666666666/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double del2,pzi;
  double sumeta, sumeta2, mineta;
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
  //    cout <<" i="<<i <<" j=" <<j <<" k=" << k <<" pn=" <<pn << " tn=" << tn << endl;
cout <<k <<",";    
  kn=k*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          // calculation of nabla square eta

          for (pn=0;pn<p;pn++)
          {
            pnn=pn*size3;
          //  mbool[i+jn+kn+pnn]=false; //firts we make is false and then if detadt>thresh then it becomes true later
            //
       //     currentind=inds[i+jn+kn+pnn];
      //      currenteta=eta[i+jn+kn+pnn];
      
}
}
}
}
}
}
