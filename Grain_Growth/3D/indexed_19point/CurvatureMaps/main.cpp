
// <source lang=cpp>

// 3D simulation with a optimized calculations
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <sstream>

// #include <pngwriter.h>

#include "peribc.cpp"
#include "readeta.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int readeta();

// ------- main routin --------
int* MGsize;
int main(int a, char** charinput)
{ 
  char* inputstr;
  inputstr=charinput[1];
  ifstream infile;
  infile.open (inputstr);
  if (infile.is_open())
    cout << "Successfully opened the file";
  else
  {
    cout << "Please write file address to the input file. Directory name should be with no space, e.g. /home/data/input.txt  :  " << endl;
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
  getline( infile, aline); getline( infile, aline); //5, 6
  cout << aline ;
  getline( infile, aline); // 7
  stringstream(aline)>> mboxsize;
  cout <<mboxsize <<endl;
  getline( infile, aline); getline( infile, aline); // 8, 9
  cout << aline ;
  getline( infile, aline); // 10
  stringstream(aline)>> nboxsize;
  cout <<nboxsize <<endl;
  getline( infile, aline); getline( infile, aline); // 11, 12
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
  double gamma=2*1.5*m;
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
  
  int p=8; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  double Lf=L;
  int i,j,k,tn;
  int n;
  double *eta;
  eta= new double[mboxsize*nboxsize*kboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*kboxsize*p];
  int *inds;
  inds= new int[mboxsize*nboxsize*kboxsize*p];
  bool* mbool;
  mbool= new bool[mboxsize*nboxsize*kboxsize*p];
  
  for (i=0;i<mboxsize*nboxsize*kboxsize*p;i++)
  {
    eta[i]=0.0000000000;
    eta2[i]=0.0000000000;
    inds[i]=0;
    mbool[i]=false;
  }
  
  double currenteta, gradeta, grad2eta, K;
  double sumterm, sumtermM, detadt, detadtM, mineta;
  int pn,psn,pind,pos,pcount,minind;
  double delx2=0.166666666666/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double del2;
  double sumeta, sumeta2;
  int currentind, indc;
  int maxind;
  double maxeta;
  //calculating processing time
  clock_t time1;
  time1=clock();
  cout << "Initialization ended. Loading the structure" <<endl;
  
  tn=0;
  int tnrelax=0;
  for (tn=resumetimestep;tn<timesteps;tn=tn+fulletawritesteps)
  {
    readeta( dirstr, inds, eta, tn, mboxsize, nboxsize, kboxsize, p);
    cout << "Reading of order parameter data for time step " << tn << " is done." << endl;
    char filenamestat[200];
    ofstream fileK;
    sprintf (filenamestat, "%sK_%d.txt",dirstr, tn);
    fileK.open (filenamestat);
    
    for (i=0;i<mboxsize*nboxsize*kboxsize*p;i++){
      if (eta[i]>thresh) {mbool[i]=true;}
      if (eta[i]>0.999) {mbool[i]=false;}
    }	  
    //----------------------Relaxing the structure ---------------
    
    for (tnrelax=1;tnrelax<timesteps1;tnrelax++)
    {
      #pragma omp parallel for num_threads(12)
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
		mbool[i+jn+kn+pnn]=false; //firts we make it false and then if detadt>thresh then it becomes true later
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
		//  //  //  First nearest neighbors
		sumeta=0;
		for (indc=0;indc<p;indc++)
		{
		  if(currentind==inds[inplus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inplus1+jn+kn+indc*size3];}
		  if(currentind==inds[inminus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inminus1+jn+kn+indc*size3];}
		  if(currentind==inds[i+jnplus1+kn+(indc)*size3]){sumeta=sumeta+eta[i+jnplus1+kn+indc*size3];}
		  if(currentind==inds[i+jnminus1+kn+indc*size3]){sumeta=sumeta+eta[i+jnminus1+kn+indc*size3];}
		  if(currentind==inds[i+jn+knplus1+indc*size3]){sumeta=sumeta+eta[i+jn+knplus1+indc*size3];}
		  if(currentind==inds[i+jn+knminus1+indc*size3]){sumeta=sumeta+eta[i+jn+knminus1+indc*size3];}
		}
		//  //  //  Second nearest neighbors
		sumeta2=0;
		for (indc=0;indc<p;indc++)
		{
		  // 0 -1 -1
		  if(currentind==inds[i+jnminus1+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnminus1+knminus1+indc*size3];}
		  // -1 0 -1
		  if(currentind==inds[inminus1+jn+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jn+knminus1+indc*size3];}
		  // 1 0 -1
		  if(currentind==inds[inplus1+jn+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jn+knminus1+indc*size3];}
		  // 0 1 -1
		  if(currentind==inds[i+jnplus1+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnplus1+knminus1+indc*size3];}
		  // -1 -1 0
		  if(currentind==inds[inminus1+jnminus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jnminus1+kn+indc*size3];}
		  // 1 -1 0
		  if(currentind==inds[inplus1+jnminus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jnminus1+kn+indc*size3];}
		  // -1 1 0
		  if(currentind==inds[inminus1+jnplus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jnplus1+kn+indc*size3];}
		  // 1 1 0
		  if(currentind==inds[inplus1+jnplus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jnplus1+kn+indc*size3];}
		  // 0 -1 1
		  if(currentind==inds[i+jnminus1+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnminus1+knplus1+indc*size3];}
		  // -1 0 1
		  if(currentind==inds[inminus1+jn+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jn+knplus1+indc*size3];}
		  // 1 0 1
		  if(currentind==inds[inplus1+jn+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jn+knplus1+indc*size3];}
		  // 0 1 1
		  if(currentind==inds[i+jnplus1+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnplus1+knplus1+indc*size3];}
		}
		del2=delx2*(sumeta2+2*sumeta-24*currenteta);
		detadtM=m*(currenteta*currenteta*currenteta-currenteta)+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta)-kappa*del2;
		
		detadt=-L*(detadtM);
		eta2[i+jn+kn+pnn]=currenteta+delt*detadt;
		// wirting data for curvature measurement 
		// make a padding for the changing order parameter so next time step it will take that into account.
		if (eta2[i+jn+kn+pnn]>thresh) 
		{
		  mbool[i+jn+kn+pnn]=true;
		  if (eta2[i+jn+kn+pnn]>0.999) {mbool[i+jn+kn+pnn]=false;}
		  pos=inplus1+jn+kn; // 1 0 0
		  mineta=eta[pos];  //the first in the table is assumed as smallest 
		  minind=0;
		  for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      pos=inminus1+jn+kn; //-1 0 0
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      pos=i+jnplus1+kn; // 0 1 0
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      pos= i+jnminus1+kn; // 0 -1 0
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      pos=i+jn+knplus1; // 0 0 1
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      pos=i+jn+knminus1; // 0 0 -1
	      mineta=eta[pos];  //the first i the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // ----------- second neighbors ------------
	      // i+jnminus1+knminus1  0 -1 -1  
	      pos=i+jnminus1+knminus1;
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inminus1+jn+knminus1 -1 0 -1 
	      pos=inminus1+jn+knminus1;
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inplus1+jn+knminus1 1 0 -1 
	      pos=inplus1+jn+knminus1;
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // i+jnplus1+knminus1 0 1 -1 
	      pos=i+jnplus1+knminus1;
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inminus1+jnminus1+kn -1 -1 0 
	      pos=inminus1+jnminus1+kn;
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inplus1+jnminus1+kn 1 -1 0
	      pos=inplus1+jnminus1+kn;
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inminus1+jnplus1+kn 1 3 2
	      pos=inminus1+jnplus1+kn; // -1 1 0
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inplus1+jnplus1+kn 3 3 2
	      pos=inplus1+jnplus1+kn; // 1 1 0
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // i+jnminus1+knplus1 2 1 3
	      pos=i+jnminus1+knplus1; // 0 -1 1 
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      // inminus1+jn+knplus1 1 2 3
	      pos=inminus1+jn+knplus1; // -1 0 1
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      //  3 2 3
	      pos=inplus1+jn+knplus1; // 1 0 1
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
	      //  2 3 3
	      pos=i+jnplus1+knplus1; // 0 1 1
	      mineta=eta[pos];  //the first in the table is assumed as smallest 
	      minind=0;
	      for (indc=0;indc<p;indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
              {
		if(inds[pos+indc*size3]==currentind){minind=indc; break;}
		if(mineta>eta[pos+indc*size3]){mineta=eta[pos+indc*size3]; minind=indc;}
	      }
	      inds[pos+minind*size3]=currentind;
	      mbool[pos+minind*size3]=true;
		}
	      }
	    }
	  }
	}
      }
      //setting eta equal to the new eta2 for the next time step
      #pragma omp parallel for num_threads(12)
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
      cout <<tnrelax << " " ;
    }
    cout <<endl;
    //-------------------------------------------------------------------------------------
    
    
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
	  for (pn=0;pn<p;pn++)
	  {
	    pnn=pn*size3;
	    currentind=inds[i+jn+kn+pnn];
	    currenteta=eta[i+jn+kn+pnn];
	    if (currenteta>0.2 && currenteta<0.8)
	    {
	      //searching for neighbors with the same index as currentind
	      //  //  //  First nearest neighbors
	      sumeta=0;
	      for (indc=0;indc<p;indc++)
	      {
		if(currentind==inds[inplus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inplus1+jn+kn+indc*size3];}
		if(currentind==inds[inminus1+jn+kn+(indc)*size3]){sumeta=sumeta+eta[inminus1+jn+kn+indc*size3];}
		if(currentind==inds[i+jnplus1+kn+(indc)*size3]){sumeta=sumeta+eta[i+jnplus1+kn+indc*size3];}
		if(currentind==inds[i+jnminus1+kn+indc*size3]){sumeta=sumeta+eta[i+jnminus1+kn+indc*size3];}
		if(currentind==inds[i+jn+knplus1+indc*size3]){sumeta=sumeta+eta[i+jn+knplus1+indc*size3];}
		if(currentind==inds[i+jn+knminus1+indc*size3]){sumeta=sumeta+eta[i+jn+knminus1+indc*size3];}
	      }
	      //  //  //  Second nearest neighbors
	      sumeta2=0;
	      for (indc=0;indc<p;indc++)
	      {
		if(currentind==inds[i+jnminus1+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnminus1+knminus1+indc*size3];}// 0 -1 -1
                if(currentind==inds[inminus1+jn+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jn+knminus1+indc*size3];}                // -1 0 -1
                if(currentind==inds[inplus1+jn+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jn+knminus1+indc*size3];} // 1 0 -1
                if(currentind==inds[i+jnplus1+knminus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnplus1+knminus1+indc*size3];} // 0 1 -1
                if(currentind==inds[inminus1+jnminus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jnminus1+kn+indc*size3];}// -1 -1 0
                if(currentind==inds[inplus1+jnminus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jnminus1+kn+indc*size3];} // 1 -1 0
                if(currentind==inds[inminus1+jnplus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jnplus1+kn+indc*size3];} // -1 1 0
                if(currentind==inds[inplus1+jnplus1+kn+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jnplus1+kn+indc*size3];}   // 1 1 0
                if(currentind==inds[i+jnminus1+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnminus1+knplus1+indc*size3];}   // 0 -1 1
                if(currentind==inds[inminus1+jn+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[inminus1+jn+knplus1+indc*size3];} // -1 0 1
                if(currentind==inds[inplus1+jn+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[inplus1+jn+knplus1+indc*size3];}   // 1 0 1
                if(currentind==inds[i+jnplus1+knplus1+(indc)*size3]){sumeta2=sumeta2+eta[i+jnplus1+knplus1+indc*size3];}     // 0 1 1
	      }
	      del2=delx2*(sumeta2+2.0*sumeta-24.0*currenteta);
	      gradeta=-0.5*sqrt(m/2.0/kappa)*(1-(2.0*currenteta-1)*(2.0*currenteta-1));
	      grad2eta=-(m/2.0/kappa)*(2.0*currenteta-1)*(1-(2.0*currenteta-1)*(2.0*currenteta-1));
	      K=(del2-grad2eta)/gradeta;
	      fileK << i << " " << j << " " << k << " " << currentind << " " << K <<endl;
	    }
	  }
	}
      }
    }
    fileK.close();
    cout << "Writing of curvature data for time step " << tn << " is done." << endl;
  }
  
  return 0;
}
