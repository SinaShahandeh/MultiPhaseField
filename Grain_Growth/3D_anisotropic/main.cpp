
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
#include "GrainsStat.cpp"
#include "resume.cpp"
#include "readeta.cpp"
//#include "WriteResults.h"

using namespace std;
int peribc();
int WriteResults();
int GrainsStat();
int resume();
int readeta();

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
  char command[300];
  sprintf(command, "mkdir %s",dirstr);
  system(command);
  
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
  cout <<fulletawritesteps <<endl;
  int Ktimestep;
  getline( infile, aline); getline( infile, aline); //32, 33
  cout << aline ;
  getline( infile, aline); // 34
  stringstream(aline)>> Ktimestep;
  cout <<Ktimestep <<endl;
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
  
  bool fullwrite=true;
  int p=8; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
  int nuclein;
  nuclein=int(mboxsize*nboxsize*kboxsize/1000/8); // ~1 percent of grid points are nuclei
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
  double* phi;
  if (fullwrite==true)
  {
    phi= new double[mboxsize*nboxsize*kboxsize];
  }
  int* Maxindfield;
  Maxindfield= new int[mboxsize*nboxsize*kboxsize];
  float* K;
  K=new float[mboxsize*nboxsize*kboxsize*p];
  
  for (i=0; i<mboxsize*nboxsize*kboxsize*p;i++){
    eta[i]=0.0;
    eta2[i]=0.0;
    inds[i]=0;
    mbool[i]=true;
    K[i]=-100.00;
  }
  for (i=0;i<mboxsize*nboxsize*kboxsize;i++){
    phi[i]=0.0;
    Maxindfield[i]=0;
  }
  // number of nucleas at the beginning of simulation
  int nn,ii,jj,kk,pp;
  double irand,jrand,prand, krand;
  MGsize=new int[nuclein];
  srand ( time(NULL) );
  // resuming from a file
  bool alreadyread=false;
if ((resumetimestep>0)&& (resumetimestep % fulletawritesteps==0) ){
  readeta( dirstr, inds, eta, resumetimestep, mboxsize, nboxsize, kboxsize, p);
  alreadyread=true;
}
  if ( (resumetimestep>0) && (resumetimestep % writingtimesteps==0) && (alreadyread==false))  resume( dirstr, inds, eta, resumetimestep, mboxsize, 
nboxsize, 
kboxsize);
 
  if (resumetimestep==0){
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
  }
  //dynamics
  double sumterm,currenteta;
  double gradeta, grad2eta;
  double detadtM;
  double detadt;
  int pn,psn,pind,pos;
  double delx2=0.166666666666/(delx*delx);
  int size3=mboxsize*nboxsize*kboxsize;
  int size2=mboxsize*nboxsize;
  int jn, kn, pnn;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double del2,pzi,Pz_applied;
  double sumeta, sumeta2, mineta;
  int currentind, indc, minind;
  int maxind;
  double maxeta;
  int indmax1, indmax2, pnl;
  int pnlmax1, pnlmax2;
  double maxeta1, maxeta2;
  double etalist[p];
  int indlist[p];
  double Lgb;


  bool isKwrite = false;
  char fileKname[200];
  ofstream fileK;
  
  //calculating processing time
  clock_t time1;
  time1=clock();
  cout << "Initialization ended." <<endl;
  
  Pz_applied=0.0;
  //optimized loop -----------------------------------------------------------------------------------------------------
  for (tn=resumetimestep+1;tn<timesteps;tn++)
  {
    // to make sure we don't calculate K if not necessary
    if (isKwrite==true){
      isKwrite= false;
    }
    if (tn % Ktimestep ==0){
      isKwrite = true;
    }
    if (tn>resumetimestep+timesteps1)
      Pz_applied=Pz;
    #pragma omp parallel for num_threads(8)
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
	    // list of all order parameters at this point
		  for (pnl = 0; pnl < p; pnl++)
		  {
			  etalist[pnl] = eta[i + jn + kn + pnl*size3];
			  indlist[pnl] = inds[i + jn + kn + pnl*size3];
		  }

		  for (pn=0;pn<p;pn++)
	    {
	      pnn=pn*size3;
	      if (mbool[i+jn+kn+pnn]==true)
	      {
		mbool[i+jn+kn+pnn]=false; //first we make it false and then if detadt>thresh then it becomes true later
		inplus1=peribc(i+1,mboxsize);
		inminus1=peribc(i-1,mboxsize);
		// here is the sum of all order parameters^2 for the point i and j and k
		sumterm=0;
		for (psn=0;psn<p;psn++)
		{
			sumterm = sumterm + eta[i + jn + kn + psn*size3] * eta[i + jn + kn + psn*size3];
		}
		// chosing the two tlargest order parameters
		maxeta1 = etalist[0];
		pnlmax1 = 0;
		for (pnl = 0; pnl < p; pnl++){
			if (etalist[pnl]>maxeta1) {
				maxeta1 = eta[pnl];
				pnlmax1 = pnl;
			}
		}
		etalist[pnlmax1] = -1; //set the max to lowest value so we can find the second largest
		maxeta2 = etalist[0]; //find the second largest
		pnlmax2 = 0;
		for (pnl = 0; pnl < p; pnl++){
			if (etalist[pnl]>maxeta2) {
				maxeta2 = etalist[pnl];
				pnlmax2 = pnl;
			}
		}
		indmax1 = indlist[pnlmax1]; //indecies of two largest order parameters
		indmax2 = indlist[pnlmax2];
		// conditional statemets to find the mobility of the interface at this point.
		Lgb = L;
		if (tn > resumetimestep + timesteps1){
			if (indmax1 <= nuclein / 2 & indmax2 <= nuclein / 2){ Lgb = L; }
			if (indmax1 > nuclein / 2 & indmax2 > nuclein / 2){ Lgb = L; }
			if (indmax1 > nuclein / 2 & indmax2 <= nuclein / 2){ Lgb = L / 100; }
			if (indmax1 <= nuclein / 2 & indmax2 > nuclein / 2){ Lgb = L / 100; }
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
		pzi=3*currenteta*(1-currenteta)*sign(detadtM)*Pz_applied;
		if (fabs(detadtM)<fabs(pzi)){
		  Lf=0;
		}
		else{
		  Lf=Lgb;
		}
		detadt=-Lf*(detadtM-pzi);
		eta2[i+jn+kn+pnn]=currenteta+delt*detadt;
		if (eta2[i+jn+kn+pnn]>1) eta2[i+jn+kn+pnn]=1;
		if (eta2[i+jn+kn+pnn]<0) eta2[i+jn+kn+pnn]=0;
		// wirting data for curvature measurement 
		if (isKwrite==true && currenteta>0.2 && currenteta<0.8){
		  gradeta=-0.5*sqrt(m/2.0/kappa)*(1-(2.0*currenteta-1)*(2.0*currenteta-1));
		  grad2eta=-(m/2.0/kappa)*(2.0*currenteta-1)*(1-(2.0*currenteta-1)*(2.0*currenteta-1));
		  K[i+jn+kn+pnn]=(float)((del2-grad2eta)/gradeta);
		}
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
      #pragma omp parallel for num_threads(8)
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
		myfile2 << phi[i+jn+kn] << " "; 
	      }
	      myfile2 << endl;
	    }
	  }
	  myfile2.close();
	}
      }
      if  (tn % writingtimesteps ==0)
      {
	// writing indexed matrix
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
	    myfile3 << Maxindfield[i+jn+kn] << " ";
	    
	    }
	    myfile3 << endl;
	  }
	}
	myfile3.close();
      }
      if  (tn % 50 ==0)
      {
	int mboolsum=0;
	for (i=0;i<size3*p;i++)
	  if (mbool[i]==true) mboolsum++;
	  cout << "Fraction of system with mbool=true is " << double(mboolsum)/double(size3)/double(p) <<endl;
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
	    }
	  }
	}
	int meanG;
	meanG=GrainsStat(Maxindfield, MGsize, nuclein, mboxsize, nboxsize, kboxsize);
	cout << "average grain size is : " << meanG <<endl;
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
      cout <<tn << endl;
      // writing all order parameters and indecis
      if (tn % fulletawritesteps ==0)
      {
	char fileetaname[200];
	for (pind=0;pind<p;pind++)
	{
	  pnn=pind*size3;
	  ofstream fileeta;
	  sprintf (fileetaname, "%sEta_%d_%d.txt",dirstr, pind, tn);
	  fileeta.open (fileetaname);
	  for (k=0;k<kboxsize;k++)
	  {
	    kn=k*size2;
	    for (j=0;j<nboxsize;j++)
	    {
	      jn=j*mboxsize;
	      for (i=0;i<mboxsize;i++)
	      {
		fileeta << eta[i+jn+kn+pnn] << " "; 
	      }
	      fileeta << endl;
	    }
	  }
	  fileeta.close();
	}
	// now for inds
	char fileindsname[200];
	for (pind=0;pind<p;pind++)
	{
	  pnn=pind*size3;
	  ofstream fileinds;
	  sprintf (fileindsname, "%sInds_%d_%d.txt",dirstr, pind, tn);
	  fileinds.open (fileindsname);
	  for (k=0;k<kboxsize;k++)
	  {
	    kn=k*size2;
	    for (j=0;j<nboxsize;j++)
	    {
	      jn=j*mboxsize;
	      for (i=0;i<mboxsize;i++)
	      {
		fileinds << inds[i+jn+kn+pnn] << " "; 
	      }
	      fileinds << endl;
	    }
	  }
	  fileinds.close();
	}
	
      }
      // open curvature file
      if (tn % Ktimestep ==0){
	sprintf (fileKname, "%sK_%d.txt",dirstr, tn);
	fileK.open (fileKname);
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
		if ((K[i+jn+kn+pnn]==-100.00)){
		  
		}
		else{
		  fileK << i << " " << j << " " << k << " " << inds[i+jn+kn+pnn] << " " << K[i+jn+kn+pnn] <<endl;
		}
		K[i+jn+kn+pnn]=-100.00;
	      }
	    }
	  }
	}
	fileK.close();
      }
      
  }
  return 0;
}


// </source>
