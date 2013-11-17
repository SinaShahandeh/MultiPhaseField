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
#include "energy.cpp"
using namespace std;
// global variables

// functions
int peribc();
int symbc();
double* shpere();
double* ParticleDistro();
double* calculatephi();
double etavolume();
int WriteResults();
double energy();

int main(int a, char** charinput)
{
  
  int var1, var2;
  //  for (var=100;var>10;var=var-10){
    char* inp1;
    inp1=charinput[2];
    var1=double(atoi(inp1));
    
    char* inp2;
    inp2=charinput[3];
    var2=double(atoi(inp2));
    
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
    //sprintf(dirstr,"%s%s/", dirstr, inputstr);
    sprintf(dirstr,"%s%s/", dirstr,inp1 );
    char command1[100];
    sprintf(command1, "mkdir %s", dirstr);
    system (command1);
    
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
    int writingtimestep;
    getline( infile, aline); getline( infile, aline); //29, 30
    cout << aline ;
    getline( infile, aline); // 31
    stringstream(aline)>> writingtimestep;
    cout <<writingtimestep <<endl;
    int fulletawritesteps;
    getline( infile, aline); getline( infile, aline); //32, 33
    cout << aline ;
    getline( infile, aline); // 34
    stringstream(aline)>> fulletawritesteps;
    cout <<fulletawritesteps <<endl;
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
    fp=double(double(var2)/1000);  //--------------------------------------------------------------- << here
    cout << fp <<endl;
    double DelGinput;
    getline( infile, aline); getline( infile, aline); //47, 48
    cout << aline ;
    getline( infile, aline); // 55
    stringstream(aline)>> DelGinput;
    
    DelGinput=double(double(var1)/1000); // ------------------------------------------------------- << here
    cout << DelGinput <<endl;
    
    double asol;
    getline( infile, aline); getline( infile, aline);
    cout << aline ;
    getline( infile, aline); // 40
    stringstream(aline)>> asol;
    //asol=double(var1)/10; 
    cout << asol <<endl;
    
    double bsol;
    getline( infile, aline); getline( infile, aline);
    cout << aline ;
    getline( infile, aline); // 40
    stringstream(aline)>> bsol;
    cout << bsol <<endl;
    
    // geometry settings
    int p=2; // phase field numbers
    int Pr=Rp;               // particle distribution properties
    double particles_fraction=fp;
    
    // model parameters
    double G[2];
    G[0]=-DelGinput;
    G[1]=+DelGinput;
    
    int i,j,k,tn;
    double *eta;
    eta= new double[mboxsize*nboxsize*kboxsize*p];
    double *eta2;
    eta2= new double[mboxsize*nboxsize*kboxsize*p];
    bool *mbool;
    mbool= new bool[mboxsize*nboxsize*kboxsize*p];
    double* phi;
    phi= new double[mboxsize*nboxsize*kboxsize];
    double* ME;
    ME = new double[mboxsize*nboxsize*kboxsize];
    
    double vol;
    double sumterm,sumtermp;
    double detadtM, detadt;
    double E;
    int pn,psn,pind;
    double  delx2=1/(delx*delx);//delx2=0.166666666666/(delx*delx);
    int size3=mboxsize*nboxsize*kboxsize;
    int size2=mboxsize*nboxsize;
    int jn, kn, pnn;
    int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
    double del2;
    bool isintmove=true;
    double grad, Mvel, error, delxgrad, detadtpast, veloc, Ps;
    int itr;
    delxgrad=1.0/2.0/delx;
    //setting initial condition  (one interface on top of the domain)
    int initialpos=int(0.20*mboxsize);
    
    for (k=0;k<kboxsize;k++)
    {
      for (j=0;j<nboxsize;j++)
      {
	for (i=0;i<mboxsize;i++)
	{
	  if (i>initialpos)
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
    ppf=ParticleDistro(mboxsize,nboxsize, kboxsize, particles_fraction, Pr);
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
    cout << "Actual particles volume fraction = " << Pf/size3 << endl;
    //dynamics
    //calculating processing time
    clock_t time1;
    time1=clock();
    cout << "Initialization ended." <<endl;
    //first loop over all the nodes to creates the obtimized matrix mbool
    for (tn=0;tn<timesteps1;tn++)
    {
      #pragma omp parallel for
      for (k=0;k<kboxsize;k++)
      {
	kn=k*size2;
	knplus1=symbc(k+1,kboxsize)*size2;
	knminus1=symbc(k-1,kboxsize)*size2;
	for (j=0;j<nboxsize;j++)
	{
	  jn=j*mboxsize;
	  jnplus1=symbc(j+1,mboxsize)*mboxsize;
	  jnminus1=symbc(j-1,mboxsize)*mboxsize;
	  for (i=0;i<mboxsize;i++)
	  {
	    inplus1=symbc(i+1,mboxsize);
	    inminus1=symbc(i-1,mboxsize);
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
	      del2=delx2*((eta[inplus1+jn+kn+pnn]+eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]+eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]+eta[i+jn+knminus1+pnn])-6*eta[i+jn+kn+pnn]);
	      sumtermp=eta[i+jn+kn+pnn]*sumterm-eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn];
	      detadtM=-m*eta[i+jn+kn+pnn]+m*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]-kappa*del2;	      
	      detadt=-L*(detadtM+gamma*sumtermp+epsilon*eta[i+jn+kn+pnn]*ppf[i+jn+kn]);
	      //cout << ppf[i+jn+kn] << endl;
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
    }
    cout << "time required for 10 time steps:" << double((clock()-time1))/double(CLOCKS_PER_SEC) << "seconds. \n";
    //optimized loop -----------------------------------------------------------------------------------------------------
    tn=0;
    vol=etavolume(eta,mboxsize, nboxsize, kboxsize);
    cout << "Initial volume is:" << vol*delx*delx*delx <<endl;
    //open volume log file
    ofstream volfile; //file containing volume data logs
    char volfilename[200];
    sprintf (volfilename, "%sVollog.log",dirstr);
    volfile.open (volfilename);
    volfile << Pf/size3 << " " << Pf/size3 << endl;
    // open energy log file
    ofstream engfile; //file containing volume data logs
    char engfilename[200];
    sprintf (engfilename, "%sEnglog.log",dirstr);
    engfile.open (engfilename);
    double pastvol=size3;
    double newvol=vol;
    do
    {
      tn=tn+1;
      time1=clock();
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
		inplus1=symbc(i+1,mboxsize);
		inminus1=symbc(i-1,mboxsize);
		// here is the sum of all order parameters^2 for the point i and j
		
		//               calculation of nabla square eta
		Mvel=0;
		error=10;
		itr=0;
		while (error>0.000001 and itr<50){
		  sumterm=0;
		for (psn=0;psn<p;psn++)
		{
		  sumterm=sumterm+eta2[i+jn+kn+psn*size3]*eta2[i+jn+kn+psn*size3];
		}
		del2=delx2*((eta2[inplus1+jn+kn+pnn]+eta2[inminus1+jn+kn+pnn]+eta2[i+jnplus1+kn+pnn]+eta2[i+jnminus1+kn+pnn]+eta2[i+jn+knplus1+pnn]+eta2[i+jn+knminus1+pnn])
		-6*eta2[i+jn+kn+pnn]);
		
		//               del2=delx2*(eta[i+jnplus1+knplus1+pnn]+eta[i+jnminus1+knplus1+pnn]+eta[i+jnminus1+knminus1+pnn]+eta[i+jnplus1+knminus1+pnn]
		//               +eta[inplus1+jn+knplus1+pnn]+eta[inplus1+jn+knminus1+pnn]+eta[inminus1+jn+knplus1+pnn]+eta[inminus1+jn+knminus1+pnn]
		//               +eta[inplus1+jnplus1+kn+pnn]+eta[inplus1+jnminus1+kn+pnn]+eta[inminus1+jnplus1+kn+pnn]+eta[inminus1+jnminus1+kn+pnn]
		//               +2*(eta[inplus1+jn+kn+pnn]+eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]+eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]+eta[i+jn+knminus1+pnn])
		//               -24*eta[i+jn+kn+pnn]);
		
		sumtermp=eta2[i+jn+kn+pnn]*sumterm-eta2[i+jn+kn+pnn]*eta2[i+jn+kn+pnn]*eta2[i+jn+kn+pnn];
		detadtM=m*(-eta2[i+jn+kn+pnn]+eta2[i+jn+kn+pnn]*eta2[i+jn+kn+pnn]*eta2[i+jn+kn+pnn])-kappa*del2;
		
		grad=delxgrad*sqrt((eta2[inplus1+jn+kn+pnn]-eta2[inminus1+jn+kn+pnn])*(eta2[inplus1+jn+kn+pnn]-eta2[inminus1+jn+kn+pnn])
		+(eta2[i+jnplus1+kn+pnn]-eta2[i+jnminus1+kn+pnn])*(eta2[i+jnplus1+kn+pnn]-eta2[i+jnminus1+kn+pnn])
		+(eta2[i+jn+knplus1+pnn]-eta2[i+jn+knminus1+pnn])*(eta2[i+jn+knplus1+pnn]-eta2[i+jn+knminus1+pnn]));
		
		
		  detadtpast=Mvel;
		  veloc=Mvel/grad;
		  if (fabs(grad)<0.00005){
		    veloc=0;
		  }
		  Ps=asol*veloc/(1+bsol*veloc*veloc);
		  detadt=-L*(detadtM+gamma*sumtermp+epsilon*eta2[i+jn+kn+pnn]*ppf[i+jn+kn]+3*eta2[i+jn+kn+pnn]*(1-eta2[i+jn+kn+pnn])*(-G[pn]+Ps));
		  Mvel=detadt; // detadt for next time step
		  error=fabs(detadtpast-Mvel);
		  itr++;
		  eta2[i+jn+kn+pnn]=eta[i+jn+kn+pnn]+delt*detadt;
		}
		
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
      volfile <<tn*delt <<" " <<vol*delx*delx*delx << endl;
      // energy(double* eta,int mboxsize, int nboxsize, int kboxsize, int p, double delx, double m, double kappa
      //      E=energy(eta, mboxsize, nboxsize, kboxsize, p, delx, m, kappa, ME);
      //      engfile <<tn*delt <<" "<<E*delx*delx<< endl;
      // write array into a file each 100 time steps
      if  (tn % writingtimestep ==0)
      {
	pastvol=newvol;
	newvol=etavolume(eta,mboxsize, nboxsize, kboxsize);
	if ((pastvol-newvol)<0.0001*size3)
	  isintmove=false;
	cout << "The volume at " << tn << " is: " << vol*delx*delx*delx <<endl;
	phi=calculatephi(eta,phi, mboxsize,nboxsize,kboxsize,p);
	int n;
	char filename[200];
	// writing
	ofstream myfile2;
	// make a string like "result_5.txt"
	sprintf (filename, "%sFullres_%d.txt",dirstr, tn);
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
	//       char zipcommand[200];
	//       sprintf (zipcommand, "zip %sFullres_%d.zip %sFullres_%d.txt",dirstr, tn, dirstr, tn);
	//       system (zipcommand);
	//       sprintf (zipcommand, "rm %sFullres_%d.txt",dirstr, tn);
	//       system (zipcommand);
	
	
	/*
	 *     // write energy matrix
	 *     char engfilestr[200];
	 *     // writing
	 *     ofstream engfile;
	 *     // make a string like "result_5.txt"
	 *     sprintf (engfilestr, "%sEngDen_%d.txt",dirstr, tn);
	 *     engfile.open (engfilestr);
	 *     for (k=0;k<kboxsize;k++)
	 *     {
	 *       kn=k*size2;
	 *       for (j=0;j<nboxsize;j++)
	 *       {
	 *         jn=j*mboxsize;
	 *         for (i=0;i<mboxsize;i++)
	 *         {
	 *           engfile << ME[i+jn+kn] << " "; 
      }
      engfile << endl;
      }
      }
      engfile.close();*/
	
	
      }
      
    }
    while (vol>(0.15)*size3 && isintmove==true);
    cout << "Final volume is:" << vol*delx*delx*delx << endl;
    volfile.close();
    engfile.close();
    //  }
    
    return 0;
}
