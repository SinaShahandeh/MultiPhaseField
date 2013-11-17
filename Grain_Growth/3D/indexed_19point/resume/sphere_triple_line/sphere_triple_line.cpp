// This file makes an initial condition of 2 grains with another spherical in the middle.
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <sstream>
using namespace std;

int main(int a, char** charinput)
{
  char* dirstr;
  dirstr=charinput[1];
  
  int i, j, k, jn, kn, pp;
  int mboxsize=200;
  int nboxsize=200;
  int kboxsize=200;
  int r=80;
  int size2=mboxsize*nboxsize;
  int* maxind;
  maxind=new int[mboxsize*nboxsize*kboxsize];
  // grain 1
  for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=0;i<mboxsize/2;i++)
      {
	maxind[i+j*mboxsize+k*size2]=2;
      }
    }
  }
  // grain 2
  for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=mboxsize/2+1;i<mboxsize;i++)
      {
	maxind[i+j*mboxsize+k*size2]=3;
      }
    }
  }
  // sphere has index of 1
  for (k=0;k<kboxsize;k++)
  {
    for (j=0;j<nboxsize;j++)
    {
      for (i=0;i<mboxsize;i++)
      {
	if (((i-mboxsize/2)*(i-mboxsize/2)+(j-nboxsize/2)*(j-mboxsize/2)+(k-kboxsize/2)*(k-kboxsize/2)) < r*r)
	{
	  maxind[i+j*mboxsize+k*size2]=1;
	}
      }
    }
  }
  
  char fileindsname[200];
  ofstream fileinds;
  sprintf (fileindsname, "%sInds_%d.txt",dirstr, 1);
  fileinds.open (fileindsname);
  for (k=0;k<kboxsize;k++)
  {
    kn=k*mboxsize*nboxsize;
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
	fileinds << maxind[i+jn+kn] << " "; 
      }
      fileinds << endl;
    }
  }
  fileinds.close();
  
}
