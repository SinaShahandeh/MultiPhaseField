// distributes particles in randomly
//#include "ParticleDistro.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

//#include "sphere.h"
#include "sphere.cpp"
// #include "peribc.cpp"
using namespace std;
double *ParticleDistro(int mboxsize, int nboxsize, int kboxsize, double Pf, int Pr)
{ 
  int particlesn=int(Pf*double(mboxsize*nboxsize*kboxsize)/(4.0/3.0*3.1415*double(Pr*Pr*Pr))); // number of particles
  // cout <<"Number of particles are" << particlesn << endl;
  
  int nn,ii,jj,kk; //random indexes
  int i,j,k; // loops
  int iind,jind,kind; //periodic indexes
  double irand,jrand,prand, krand;
  double *ppfi;
  double *sphi;
  sphi=sphere(Pr);
  ppfi=new double[mboxsize*nboxsize*kboxsize];
  for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=0;i<mboxsize;i++)
      {
        ppfi[i+j*mboxsize+k*mboxsize*nboxsize]=0;
      }
    }
  }
  
  if (Pf ==-1){
    cout <<"hi, I am here";
    // 0.1* mboxsize is initial position of the interface + particle radius + interface thickness threshold
    ii=int(mboxsize/2);//int(0.1*mboxsize+2*Pr+10);
    jj=int(nboxsize/2);
    kk=int(kboxsize/2);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
  }
  // if Pf is negative then it mean there is only one particle in the middle
  if (Pf ==-4){
    // 0.1* mboxsize is initial position of the interface + particle radius + interface thickness threshold
    ii=int(0.2*mboxsize+Pr+40);
    jj=int(nboxsize/2);
    kk=int(kboxsize/2);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
    // for 4 particles on hex --------------------------------------------------------------------------
        ii=int(0.2*mboxsize+Pr+40);
    jj=0;
    kk=int(kboxsize/2);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
        ii=int(0.2*mboxsize+Pr+40);
    jj=int(nboxsize/4);
    kk=0;
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
        ii=int(0.2*mboxsize+Pr+40);
    jj=int(3*nboxsize/4);
    kk=0;
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
    //------------------------------------------------------------
  }
  
  
    if (Pf ==-6){
    // 0.1* mboxsize is initial position of the interface + particle radius + interface thickness threshold
    ii=int(0.2*mboxsize+Pr+40);
//     jj=int(nboxsize*0.25);
//     kk=int(kboxsize*29/88);
        jj=int(nboxsize*0.63);
    kk=int(kboxsize*22/88);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
    // for 6 particles on hex --------------------------------------------------------------------------
        ii=int(0.2*mboxsize+Pr+40);
//     jj=int(nboxsize*0.25);
//     kk=int(kboxsize*58/88);
        jj=int(nboxsize*0.38);
    kk=int(kboxsize*22/88);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
        ii=int(0.2*mboxsize+Pr+40);
//     jj=int(nboxsize*0.50);
//     kk=int(kboxsize*73/88);
	    jj=int(nboxsize*0.25);
    kk=int(kboxsize*44/88);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
    ii=int(0.2*mboxsize+Pr+40);
//     jj=int(nboxsize*0.75);
//     kk=int(kboxsize*58/88);
        jj=int(nboxsize*0.38);
    kk=int(kboxsize*65/88);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
        ii=int(0.2*mboxsize+Pr+40);
//     jj=int(nboxsize*0.75);
//     kk=int(kboxsize*29/88);
	    jj=int(nboxsize*0.63);
    kk=int(kboxsize*65/88);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
        ii=int(0.2*mboxsize+Pr+40);
//     jj=int(nboxsize*0.50);
//     kk=int(kboxsize*15/88);
	    jj=int(nboxsize*0.75);
    kk=int(kboxsize*44/88);
    for (i=-Pr;i<Pr;i++)
    {
      iind=peribc(ii+i,mboxsize);
      for (j=-Pr;j<Pr;j++)
      {
        jind=peribc(jj+j,nboxsize);
        for (k=-Pr;k<Pr;k++) 
        {
          kind=peribc(kk+k,kboxsize);
          ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr];
        }
      }
    }
    //------------------------------------------------------------
  }
  
  
  
  
  
  
  if (Pf >0){
    srand ( time(NULL) );
    for (nn=0;nn<particlesn;nn++){
      irand=rand();
      jrand=rand();
      krand=rand();
      prand=rand();
      // these are random index of the top corner of a box contaning the particle.
      // The box is generated by sphere.cpp 
      ii=int((mboxsize*irand)/RAND_MAX);
      jj=int((nboxsize*jrand)/RAND_MAX);
      kk=int((kboxsize*krand)/RAND_MAX);
      for (i=-Pr;i<Pr;i++)
      {
        iind=peribc(ii+i,mboxsize);
        for (j=-Pr;j<Pr;j++)
        {
          jind=peribc(jj+j,nboxsize);
          for (k=-Pr;k<Pr;k++)
          {
            kind=peribc(kk+k,kboxsize);
            ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=sphi[(i+Pr)+(j+Pr)*2*Pr+(k+Pr)*2*Pr*2*Pr]+ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize];
            if (ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]>1)
              ppfi[(iind)+(jind)*mboxsize+(kind)*mboxsize*nboxsize]=1;
          }
        }
      }
    }
  }
  return ppfi;
}
