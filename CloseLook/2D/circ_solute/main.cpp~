// A circular shape grain moving with curvature and pinning force also acting on it.
//parameters :  a.out [dome radius] [friction force]
// sparce structure with boolian variable.
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>

// #include <pngwriter.h>
#include "symbc.cpp"
//#include "WriteResults.h"
#include "etavolume.cpp"
#include "createcircle.cpp"
#include "createdome.cpp"
using namespace std;
int symbc();
int WriteResults();
double* createcircle();
double etavolume();

double sign(double x){
  if (x>0){return 1;}
  if (x<0){return -1;}
  if (x=0){return 0;}
  return 0;
}

// ------- main routin --------
int main(int a, char** charinput)
{
  
  char* dirstr;
  //sprintf(dirstr, "/media/disk/sim_res/test/");
  dirstr=charinput[1];
  cout <<"Data will be saved in this Directory:---->" << dirstr <<endl;
  
  ofstream myfile ;
  // model parameters
  double delt=0.01;
  int timesteps1=500;
  int writingtimestep=10000;
  
  double L=1;
  double alpha=4.0;
  double beta=alpha;
  double gamma=alpha*1.5;
  double kappa=8.0;
  
  double Lf;
  double Pf;
  // solte drag parameters
  double asol=2.00;
  double bsol;
  
  char* inp;
  inp=charinput[2];
  bsol=double(atoi(inp));
  cout <<"value of b is:---->" << bsol <<endl;
  
  
  int i,j,tn;
  // geometry settings
  // geometry settings

  double delx=1.0;      // length unit per pixel
  int r=int(atoi(charinput[3])/delx); //dome radius in pixel
  cout << "r=" << r << endl;
  cout << "delx=" << delx << endl;
  int p=2; // phase field numbers
  int mboxsize=20+2*r; // x axis in pixels
  int nboxsize=20+2*r; // y axis
  double endvol=0.005;
  double thresh=0.000000001; //threshold value for choosing active nodes
  double* eta;
  eta= new double[mboxsize*nboxsize*p];
  double *eta2;
  eta2= new double[mboxsize*nboxsize*p];
  bool *mbool;
  mbool= new bool[mboxsize*nboxsize*p];
  
  double* Mdetadt;
  Mdetadt=new double[mboxsize*nboxsize*p];
  
  for (i=0;i<mboxsize*nboxsize*p;i++)
  {
    eta[i]=0;
    eta2[i]=0;
    Mdetadt[i]=0;
  }
  int nn,ii,jj;
  int inplus1, inminus1, jnplus1,jnminus1, knplus1, knminus1;
  double irand,jrand,prand;
  
  eta=createcircle(eta, mboxsize, nboxsize,r);
  //  eta=createdome(eta, mboxsize, nboxsize,r);
  // particles distribution specification
  
  //dynamics
  double sumterm,sumtermp;
  double detadtM;
  double detadt;
  int pn,psn,pind;
  double delx2=2.00000000/3.000000/(delx*delx);
  double veloc;
  double grad;
  int size=mboxsize*nboxsize;
  int jn, pnn;
  double del2;
  double err;
  int ittr;
  //calculating processing time
  clock_t time1;
  time1=clock();
  
  ofstream volfile;
  
  char filename[200];
  sprintf (filename,"%s/%d/vollog.log", dirstr,atoi(charinput[2]));
  volfile.open (filename);
  double vol=mboxsize*nboxsize, initvol;
  double pastvol=vol;
  initvol=etavolume(eta,mboxsize, nboxsize)*delx*delx;
  cout << "initial volume = " << initvol << endl;
  
  tn=0;
  while (vol>endvol*initvol)
  {
    tn=tn+1;
    vol=etavolume(eta,mboxsize, nboxsize)*delx*delx;
    volfile << tn*delt << " " << vol << endl;
    
    if  ((tn-1) % writingtimestep ==0)
    {
      double *phi;
      phi= new double[mboxsize*nboxsize];
      // making the phi array
      for (i=0;i<mboxsize;i++)
      {
        for (j=0;j<nboxsize;j++)
        {
          jn=j*mboxsize;
          phi[i+jn]=0;
          for (pind=0;pind<p;pind++)
          {
            pnn=pind*size;
            phi[i+jn]=phi[i+jn]+eta[i+jn+pnn]*eta[i+jn+pnn];
          }
        }
      }
      // writing
      ofstream myfile;
      char filename[200];
      sprintf (filename, "%s/%d/Fullres_%d.txt",dirstr, atoi(charinput[2]), tn);
      myfile.open (filename);
      for (i=0;i<mboxsize;i++)
      {
        for (j=0;j<nboxsize;j++)
        {
          myfile << eta[i+j*mboxsize+0*size] << " "; 
        }
        myfile << "\n";
      }
      myfile.close();
      cout << "Time required for time step number [ " << tn << " ] is :" << double(clock()-time1)/double(CLOCKS_PER_SEC) << " seconds. r= " << sqrt(vol/3.1415) <<"\n";
      time1=clock();
      if ((pastvol-vol)<5){
        cout << "Interface is not moving!!" << endl;
        //         break;
      }
      pastvol=vol;
      /*
       *      ofstream myfile3;
       *      char filename3[200];
       *      sprintf (filename3, "%s/Mdetadt_%d.txt",dirstr, tn);
       *      myfile3.open (filename3);
       *      for (i=0;i<mboxsize;i++)
       *      {
       *        for (j=0;j<nboxsize;j++)
       *        {
       *          myfile3 << Mdetadt[i+j*mboxsize+0*size] << " "; 
       }
       myfile3 << "\n";
       }
       myfile3.close();
       */
       }
       
       
       #pragma omp parallel for
       for (j=0;j<mboxsize;j++)
       {
         jn=j*mboxsize;
         jnplus1=symbc(j+1,mboxsize)*mboxsize;
         jnminus1=symbc(j-1,mboxsize)*mboxsize;
         for (i=0;i<nboxsize;i++)
         {
           inplus1=symbc(i+1,nboxsize);
           inminus1=symbc(i-1,nboxsize);
           // calculation of nabla square eta
           for (pn=0;pn<p;pn++)
           {
             pnn=pn*size;
             sumterm=0;// here is the sum of all order parameters^2 for the point i and j
             for (psn=0;psn<p;psn++)
             {
               sumterm=sumterm+eta[i+jn+psn*size]*eta[i+jn+psn*size];
             }
             err=10; ittr=0;
             // del2=delx2*((eta[inplus1+jn+pnn]+eta[inminus1+jn+pnn]+eta[i+jnplus1+pnn]+eta[i+jnminus1+pnn])-4*eta[i+jn+pnn]);
             del2=delx2*((eta[inplus1+jn+pnn]+eta[inminus1+jn+pnn]+eta[i+jnplus1+pnn]+eta[i+jnminus1+pnn])
             +0.25*(eta[inplus1+jnplus1+pnn]+eta[inminus1+jnminus1+pnn]+eta[inminus1+jnplus1+pnn]+eta[inplus1+jnminus1+pnn])
             -5*eta[i+jn+pnn]);
             //             gradx=delxgrad*((eta[inplus1+jn+pnn]-eta[inminus1+jn+pnn]));
             //             grady=delxgrad*((eta[i+jnplus1+pnn]-eta[i+jnminus1+pnn]));
             //             grad=sqrt(gradx*gradx+grady*grady);
	    //find gradient of eta for the velocity
             grad=0.5*(1-(2*eta[i+jn+pnn]-1)*(2*eta[i+jn+pnn]-1))*sqrt(alpha/2/kappa);
             sumtermp=eta[i+jn+pnn]*sumterm-(eta[i+jn+pnn]*eta[i+jn+pnn]*eta[i+jn+pnn]);
             detadtM=-alpha*eta[i+jn+pnn]+beta*(eta[i+jn+pnn]*eta[i+jn+pnn]*eta[i+jn+pnn])+2*gamma*sumtermp-kappa*del2;
            //iteration over detadt to get the accurate velocity.
             while (err>0.000001 && ittr<50)
             {
               veloc=Mdetadt[i+jn+pnn]/grad;
               if (fabs(grad)<0.00005){
                 veloc=0;
               }
               Pf=asol*veloc/(1+bsol*veloc*veloc);
               detadt=-L*(detadtM+3*eta[i+jn+pnn]*(1-eta[i+jn+pnn])*Pf); //finding new detadt with the presence of drag, new detadt changed velocity so we need to iterate.
               err=fabs(Mdetadt[i+jn+pnn]-detadt);
               ittr++;
               Mdetadt[i+jn+pnn]=(detadt+Mdetadt[i+jn+pnn])/2.00; // detadt for next time step
             }
             eta2[i+jn+pnn]=eta[i+jn+pnn]+delt*detadt;
             // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
             if (eta2[i+jn+pnn]>1) eta2[i+jn+pnn]=1;
             if (eta2[i+jn+pnn]<0) eta2[i+jn+pnn]=0;
           }
         }
       }
       //setting eta equal to the new eta2 for the next time step
       for (pind=0;pind<p;pind++)
       {
         pnn=pind*size;
         //       maxdetadt=0;
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
    volfile.close();
    cout << "Calculation for r=" << r << " is finished." << endl;
    return 0;
  }
  
  
