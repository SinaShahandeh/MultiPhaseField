
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
  cout <<kboxsize <<endl;
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
  
  int p=5; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  double Lf=L;
  int i,j,tn;
  int n;
  double *eta;
  eta= new double[mboxsize*nboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*p];
  int *inds;
  inds= new int[mboxsize*nboxsize*p];
  bool* mbool;
  mbool= new bool[mboxsize*nboxsize*p];
  for (i=0;i<mboxsize*nboxsize*p;i++)
  {
    eta[i]=0.00;
    eta2[i]=0;
  }
  for (i=0;i<mboxsize*nboxsize;i++)
  {
    mbool[i]=true;
  }
  double mineta;
  double sumterm;
  double detadtM;
  double detadt;
  
  double phii;
  double currenteta, gradeta, grad2eta, K;
  int pn,psn,pind,pos,pcount;
  double delx2=2.0000000/3.0000000/(delx*delx);
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn, indscn;
  int inplus1, inminus1, jnplus1,jnminus1;
  double del2;
  double sumeta, sumeta2;
  int currentind, indc;
  int maxind, minind;
  //calculating processing time
  clock_t time1;
  time1=clock();
  cout << "Initialization ended." <<endl;
  tn=0;
  for (tn=resumetimestep;tn<timesteps;tn=tn+fulletawritesteps)
  {
    readeta( dirstr, inds, eta, tn, mboxsize, nboxsize, kboxsize, p);
    cout << "Reading of order parameter data for time step " << tn << " is done." << endl;
   
    
    // ---------------- Relaxing the structure -------------------------------------------------------------------------------------------------
    
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      jnplus1=peribc(j+1,nboxsize)*mboxsize;
      jnminus1=peribc(j-1,nboxsize)*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
	inplus1=peribc(i+1,mboxsize);
	inminus1=peribc(i-1,mboxsize);
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
    int tn2=0;
    //optimized loop -----------------------------------------------------------------------------------------------------
    while (tn2<timesteps1) //newsize-pastsize>0.0000000000001
  {
    tn2++;
    cout << tn2 <<endl;
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
	  inplus1=peribc(i+1,mboxsize);
	  inminus1=peribc(i-1,mboxsize);
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
	    detadt=-L*(detadtM);
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
  }
  
  // --------------------------------------------------------------------------------------------------------------------
  
  
  
  char filenamestat[200];
  ofstream fileK;
  sprintf (filenamestat, "%sK_%d.txt",dirstr, tn);
  fileK.open (filenamestat);
  for (j=0;j<nboxsize;j++)
  {
    jn=j*mboxsize;
    jnplus1=peribc(j+1,nboxsize)*mboxsize;
    jnminus1=peribc(j-1,nboxsize)*mboxsize;
    for (i=0;i<mboxsize;i++)
    {
      inplus1=peribc(i+1,mboxsize);
      inminus1=peribc(i-1,mboxsize);
      // calculation of nabla square eta
      for (pn=0;pn<p;pn++)
      {
	pnn=pn*size2;
	//
	currentind=inds[i+jn+pnn];
	currenteta=eta[i+jn+pnn];
	if (currenteta>0.2 && currenteta<0.8)
	{
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
	  gradeta=-0.5*sqrt(m/2.0/kappa)*(1-(2.0*currenteta-1)*(2.0*currenteta-1));
	  grad2eta=-(m/2.0/kappa)*(2.0*currenteta-1)*(1-(2.0*currenteta-1)*(2.0*currenteta-1));
	  K=(del2-grad2eta)/gradeta;
	  fileK << i << " " << j << " " << currentind << " " << K <<endl;
	}
      }
    }
  }
  fileK.close();
  cout << "Writing of curvature data for time step " << tn << " is done." << endl;
  }
  
  return 0;
}
