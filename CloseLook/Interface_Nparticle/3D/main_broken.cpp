// 3D simulation with a optimized calculations
// for an interface moving in array of particles
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <sstream>
// #include <pngwriter.h>
//#include "peribc.h"
//#include "sphere.h"
//#include "ParticleDistro.h"
//#include "WriteResults.h"

#include "peribc.cpp"
#include "symbc.cpp"
#include "ParticleDistro.cpp"
#include "etavolume.cpp"
#include "calculatephi.cpp"

using namespace std;

// functions
int peribc();
int symbc();
double* shpere();
double* ParticleDistro();
double* calculatephi();
double etavolume();
int WriteResults();

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
  double gamma=2*1.50*m;
  double kappa;
  getline( infile, aline); getline( infile, aline); //41, 42
  cout << aline ;
  getline( infile, aline); // 43
  stringstream(aline)>> kappa;
  cout <<kappa <<endl;
  double epsilon;
  getline( infile, aline); getline( infile, aline); //44, 45
  cout << aline ;
  getline( infile, aline); // 46
  stringstream(aline)>> epsilon;
  cout << epsilon <<endl;
  int Rp;
  getline( infile, aline); getline( infile, aline); //47, 48
  cout << aline ;
  getline( infile, aline); // 49
  stringstream(aline)>> Rp;
  cout << Rp <<endl;
  double fp;
  getline( infile, aline); getline( infile, aline); //47, 48
  cout << aline ;
  getline( infile, aline); // 52
  stringstream(aline)>> fp;
  cout << fp <<endl;
  double DelGinput;
  getline( infile, aline); getline( infile, aline); //47, 48
  cout << aline ;
  getline( infile, aline); // 55
  stringstream(aline)>> DelGinput;
  cout << DelGinput <<endl;
  
  double DelG[2];
  DelG[0]=DelGinput;
  DelG[1]=-DelGinput;
  
  int i,j,k,tn;
  int p=2; // phase field numbers
  double *eta;
  eta= new double[mboxsize*nboxsize*kboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*kboxsize*p];
  bool *mbool;
  mbool= new bool[mboxsize*nboxsize*kboxsize*p];
  double* phi;
  phi= new double[mboxsize*nboxsize*kboxsize];
  float* writedata;
  writedata= new float[mboxsize*nboxsize*kboxsize];
  double vol;
  double sumterm, currenteta;
  double detadtM;
  double detadt;
  int pn,psn,pind;
  double delx2=0.166666666666/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double del2;
  //setting initial condition  (one interface on top of the domain)
  int initialpos=int(0.50*mboxsize);
  
  for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=0;i<mboxsize;i++)
      {
        if (k<initialpos)
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
  // particles distribution specification
  double* ppf;
  //ppf=ParticleDistro(int particlesn, int Pr, int mboxsize, int nboxsize, int kboxsize);
  ppf=ParticleDistro(mboxsize,nboxsize, kboxsize, fp, Rp);
  double Pf=0; //actual particles vlocume fraction
  for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=0;i<mboxsize;i++)
      {
        if (ppf[i+j*mboxsize+k*mboxsize*nboxsize]==1)
        {Pf=Pf+1;}
      }
    }
  }
  Pf=Pf/size3;
  cout << "Actual particles volume fraction = " << Pf << endl;
  //dynamics
  //calculating processing time
  clock_t time1;
  time1=clock();
  cout << "Initialization ended." <<endl;
  //first loop over all the nodes to creates the obtimized matrix mbool
  for (tn=0;tn<timesteps1;tn++)
  {
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
            currenteta=eta[i+jn+kn+pnn];
            del2=delx2*(eta[i+jnplus1+knplus1+pnn]+eta[i+jnminus1+knplus1+pnn]+eta[i+jnminus1+knminus1+pnn]+eta[i+jnplus1+knminus1+pnn]
            +eta[inplus1+jn+knplus1+pnn]+eta[inplus1+jn+knminus1+pnn]+eta[inminus1+jn+knplus1+pnn]+eta[inminus1+jn+knminus1+pnn]
            +eta[inplus1+jnplus1+kn+pnn]+eta[inplus1+jnminus1+kn+pnn]+eta[inminus1+jnplus1+kn+pnn]+eta[inminus1+jnminus1+kn+pnn]
            +2*(eta[inplus1+jn+kn+pnn]+eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]+eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]+eta[i+jn+knminus1+pnn])
            -24*currenteta);
            detadtM=m*(currenteta*currenteta*currenteta-currenteta)+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta)-kappa*del2;
            detadt=-L*(detadtM+epsilon*currenteta*ppf[i+jn+kn]);
            
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
    //cout << tn << "\n";
    
    calculatephi(eta,phi,mboxsize,nboxsize,kboxsize,p);
      int n;
      char filename[200];
      // writing
      ofstream myfile2;
      // make a string like "result_5.txt"
      n=sprintf (filename, "%sFullres_%d.txt", dirstr, tn);
      myfile2.open (filename);
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            myfile2 << phi[i+jn+kn] << " "; 
          }
          myfile2 << endl;
        }
      }
      myfile2.close();      
      cout << "The volume at time step " << tn << " is: " << vol <<endl;
      
    
    
  }
  cout << "time required for 10 time steps:" << double((clock()-time1))/double(CLOCKS_PER_SEC) << "seconds. \n";
  //optimized loop -----------------------------------------------------------------------------------------------------
  tn=0;
  vol=etavolume(eta,mboxsize, nboxsize, kboxsize);
  cout << "Initial volume is:" << vol <<endl;
  ofstream volfile; //file containing volume data logs
  int nvol;
  char volfilename[200];
  nvol=sprintf(volfilename,"%sVollog.log", dirstr);
  volfile.open (volfilename);
  
  while (vol<(0.95-Pf)*size3)
  {
    tn=tn+1;
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
            if (4>3)//(mbool[i+jn+kn+pnn]==true)
            {
              inplus1=peribc(i+1,mboxsize);
              inminus1=peribc(i-1,mboxsize);
              // here is the sum of all order parameters^2 for the point i and j
              sumterm=0;
              for (psn=0;psn<p;psn++)
              {
                sumterm=sumterm+eta[i+jn+kn+psn*size3]*eta[i+jn+kn+psn*size3];
              }
              pnn=pn*size3;
              currenteta=eta[i+jn+kn+pnn];
              del2=delx2*(eta[i+jnplus1+knplus1+pnn]+eta[i+jnminus1+knplus1+pnn]+eta[i+jnminus1+knminus1+pnn]+eta[i+jnplus1+knminus1+pnn]
              +eta[inplus1+jn+knplus1+pnn]+eta[inplus1+jn+knminus1+pnn]+eta[inminus1+jn+knplus1+pnn]+eta[inminus1+jn+knminus1+pnn]
              +eta[inplus1+jnplus1+kn+pnn]+eta[inplus1+jnminus1+kn+pnn]+eta[inminus1+jnplus1+kn+pnn]+eta[inminus1+jnminus1+kn+pnn]
              +2*(eta[inplus1+jn+kn+pnn]+eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]+eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]+eta[i+jn+knminus1+pnn])
              -24*currenteta);
              detadtM=m*(currenteta*currenteta*currenteta-currenteta)+gamma*(currenteta*sumterm-currenteta*currenteta*currenteta)-kappa*del2;
              detadt=-L*(detadtM+epsilon*currenteta*ppf[i+jn+kn]-3*(currenteta*(1-currenteta)*DelG[pn]));
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
    vol=etavolume(eta,mboxsize, nboxsize, kboxsize);
    //cout << "The volume is:" << vol <<endl;
    volfile <<tn*delt <<"    " <<vol*delx*delx*delx <<endl;
    // write array into a file each 100 time steps
    if  (tn % fulletawritesteps ==0)
    {
      calculatephi(eta,phi,mboxsize,nboxsize,kboxsize,p);
      int n;
      char filename[200];
      // writing
      ofstream myfile2;
      // make a string like "result_5.txt"
      n=sprintf (filename, "%sFullres_%d.txt", dirstr, tn);
      myfile2.open (filename);
      for (k=0;k<kboxsize;k++)
      {
        kn=k*size2;
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          for (i=0;i<mboxsize;i++)
          {
            myfile2 << phi[i+jn+kn] << " "; 
          }
          myfile2 << endl;
        }
      }
      myfile2.close();      
      cout << "The volume at time step " << tn << " is: " << vol <<endl;
      
      
      ofstream binphi;
      char binphiname[200];
      sprintf (binphiname, "%sFullres_%d.bin", dirstr, tn);
      binphi.open (binphiname, ios::out | ios::binary);
      for (i=0;i<size3;i++)
        writedata[i]=float(phi[i]);
      binphi.write((char*) writedata,4*size3 );
      binphi.close();
      
    }
    
  }
  
  return 0;
}
