
// calculating energy of system from a saved eta and ind
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>

#include "peribc.cpp"

#include "readstructure.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int readstructure();

// ------- main routin --------
int* MGsize;
int main(int a, char** charinput)
{ 
  char* strinput;
  strinput=charinput[1];
  cout <<"This file is read as input:---->" << strinput <<endl;
  
  int readtimestep;
  readtimestep=atoi(charinput[2]);
  
  ifstream infile;
  infile.open (strinput);
  if (infile.is_open())
    cout << "Successfully opened the input file";
  else
  {
    cout << "Please write the address to the input file. Directory name should be with no space, e.g. /home/data/input.txt  :  " << endl;
    char dirinput[300];
    cin >> dirinput;
    infile.open (dirinput);
    if (infile.is_open()==0) {
      cout << "I could not open the file specified. The program can not continue.";
      exit(1);
    }
  }
  string aline;
  getline( infile, aline);   getline( infile, aline);//header line
  cout <<  "Reading input.txt file:";
  getline( infile, aline); //3
  cout << aline ; //Reading and saving directory:
  getline( infile, aline); //4
  stringstream value(aline);
  char dirstr[200];
  value >>dirstr;
  cout <<"Data will be saved in this Directory:---->" << dirstr <<endl;
  // geometry settings
  int mboxsize=1;
  int nboxsize=1;
  int kboxsize=1;
  getline( infile, aline); getline( infile, aline); //5 ,6
  cout << aline ;
  getline( infile, aline); // 7
  stringstream(aline)>> mboxsize;
  cout <<mboxsize <<endl;
  getline( infile, aline); getline( infile, aline); // 8,9
  cout << aline ;
  getline( infile, aline); // 10
  stringstream(aline)>> nboxsize;
  cout <<nboxsize <<endl;
  getline( infile, aline); getline( infile, aline); // 11,12
  cout << aline ;
  getline( infile, aline); // 13
  stringstream(aline)>> kboxsize;
  cout <<nboxsize <<endl;
  double delx;      // length unit per pixel
  getline( infile, aline); getline( infile, aline); // 14, 15
  cout << aline ;
  getline( infile, aline); // 16
  stringstream(aline)>> delx;
  cout <<delx <<endl;
  double thresh; //threshold value for choosing active nodes
  getline( infile, aline); getline( infile, aline); //17, 18
  cout << aline ;
  getline( infile, aline); // 19
  stringstream(aline)>> thresh;
  cout <<thresh <<endl;
  double delt;
  getline( infile, aline); getline( infile, aline); //20, 21
  cout << aline;
  getline( infile, aline); // 22
  stringstream(aline)>> delt;
  cout <<delt <<endl;
  // number of time steps for first calculation without boolean matrix.
  int timesteps1;
  getline( infile, aline); getline( infile, aline); //23, 24
  cout << aline ;
  getline( infile, aline); // 25
  stringstream(aline)>> timesteps1;
  cout <<timesteps1 <<endl;
  int timesteps;
  getline( infile, aline); getline( infile, aline); //26, 27
  cout << aline ;
  getline( infile, aline); // 28
  stringstream(aline)>> timesteps;
  cout <<timesteps <<endl;
  int writingtimesteps;
  getline( infile, aline); getline( infile, aline); //29, 30
  cout << aline ;
  getline( infile, aline); // 31
  stringstream(aline)>> writingtimesteps;
  cout <<writingtimesteps <<endl;
  int fulletawritesteps;
  getline( infile, aline); getline( infile, aline); //32, 33
  cout << aline ;
  getline( infile, aline); // 34
  stringstream(aline)>> fulletawritesteps;
  cout <<writingtimesteps <<endl;
  double L;
  getline( infile, aline); getline( infile, aline); //35, 36
  cout << aline ;
  getline( infile, aline); // 37
  stringstream(aline)>> L;
  cout <<L <<endl;
  double m;
  getline( infile, aline); getline( infile, aline); //38, 39
  cout << aline ;
  getline( infile, aline); // 40
  stringstream(aline)>> m;
  cout <<m <<endl;
  double gamma=1.5;
  double kappa;
  getline( infile, aline); getline( infile, aline); //41, 42
  cout << aline ;
  getline( infile, aline); // 43
  stringstream(aline)>> kappa;
  cout <<kappa <<endl;
  double Pz;
  getline( infile, aline); getline( infile, aline); //44, 45
  cout << aline ;
  getline( infile, aline); // 46
  stringstream(aline)>> Pz;
  cout <<Pz <<endl;
  int resumetimestep;
  getline( infile, aline); getline( infile, aline); //47, 48
  cout << aline ;
  getline( infile, aline); // 49
  stringstream(aline)>> resumetimestep;
  cout << resumetimestep <<endl;
  
  bool fullwrite=true;
  int p=8; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  
  int i,j,k,tn;
  int n;
  double *eta;
  eta= new double[mboxsize*nboxsize*kboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*kboxsize*p];
  int *inds;
  inds= new int[mboxsize*nboxsize*kboxsize*p];
  
  double* E;
  E= new double[mboxsize*nboxsize*kboxsize];
  
  
  for (i=0; i<mboxsize*nboxsize*kboxsize*p;i++){
    eta[i]=0.0;
    eta2[i]=0.0;
    inds[i]=0;
  }
  for (i=0;i<mboxsize*nboxsize*kboxsize;i++){
    E[i]=0.0;
  }
  // number of nucleas at the beginning of simulation
  int nn,ii,jj,kk,pp;
  double sumterm,currenteta;
  double detadtM;
  double detadt;
  int pn,psnj,psni;
  double delx2=0.166666666666/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double sumeta, sumeta2;
  int currentind, indc;
  double gradient, gradientsum, firstsum;
  // resuming from a file
  if (resumetimestep>0)  readstructure(dirstr, inds, eta, readtimestep, mboxsize, nboxsize, kboxsize,p);
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
        inplus1=peribc(i+1,mboxsize);
        inminus1=peribc(i-1,mboxsize);
        // here is the sum of all order parameters^2 for the point i and j
        sumterm=0;
        for (psni=0;psni<p;psni++)
        {
          for (psnj=psni+1;psnj<p;psnj++)
          {
            sumterm=sumterm+eta[i+jn+kn+psnj*size3]*eta[i+jn+kn+psnj*size3]*eta[i+jn+kn+psni*size3]*eta[i+jn+kn+psni*size3];
          }
        }
        // calculation of gradient of eta
        //searching for neighbors with the same index as currentind
        //  //  //  First nearest neighbors
        gradientsum=0;
        for (psni=0;psni<p;psni++)
        {
          currentind=inds[i+jn+kn+psni*size3];
          currenteta=eta[i+jn+kn+psni*size3];
          sumeta=0;
          for (indc=0;indc<p;indc++) //adds the i+1 etas and substract i-1 etas to calculate gradient
            { 
              if(currentind==inds[inplus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inplus1+jn+kn+indc*size3];}
              if(currentind==inds[inminus1+jn+kn+(indc)*size3]){sumeta=sumeta-eta[inminus1+jn+kn+indc*size3];}
              if(currentind==inds[i+jnplus1+kn+(indc)*size3]){sumeta=sumeta+eta[i+jnplus1+kn+indc*size3];}
              if(currentind==inds[i+jnminus1+kn+indc*size3]){sumeta=sumeta-eta[i+jnminus1+kn+indc*size3];}
              if(currentind==inds[i+jn+knplus1+indc*size3]){sumeta=sumeta+eta[i+jn+knplus1+indc*size3];}
              if(currentind==inds[i+jn+knminus1+indc*size3]){sumeta=sumeta-eta[i+jn+knminus1+indc*size3];}
            }
            gradient=sumeta/2/delx; //gradient for eta with the index currentind in the layer psni
            gradientsum=gradientsum+gradient*gradient;
        }
        firstsum=0;
        for (psni=0;psni<p;psni++)
        {
          firstsum=firstsum+(-eta[i+jn+kn+psni*size3]*eta[i+jn+kn+psni*size3]/2+eta[i+jn+kn+psni*size3]*eta[i+jn+kn+psni*size3]*eta[i+jn+kn+psni*size3]*eta[i+jn+kn+psni*size3]/4);
        }
        E[i+jn+kn]=m*(firstsum+gamma*sumterm+0.250)+kappa/2*gradientsum;
      }
    }
  }
  
  
  
  // total energy of the system
  double TotalE=0;
  for (i=0;i<mboxsize*nboxsize*kboxsize;i++)
  {
    TotalE=TotalE + E[i];
  }
  
  cout <<"Total energy of the structure = "  << TotalE <<endl;
  
  //writing energy density map
  if (fullwrite == true)
  {
    char filename[200];
    ofstream myfile2;
    // make a string like "result_5.txt"
    n=sprintf (filename, "%sEnergy_%d.txt",dirstr, resumetimestep);
    myfile2.open (filename);
    for (k=0;k<kboxsize;k++)
    {
      kn=k*size2;
      for (j=0;j<nboxsize;j++)
      {
        jn=j*mboxsize;
        for (i=0;i<mboxsize;i++)
        {
          myfile2 << E[i+jn+kn] << " "; 
        }
        myfile2 << endl;
      }
    }
    myfile2.close();
    
  }
  
  return 0;
}
