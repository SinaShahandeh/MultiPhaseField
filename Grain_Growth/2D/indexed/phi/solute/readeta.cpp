
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int resume(char* dirstr, int* inds, double* eta, int resumetimestep, int mboxsize, int nboxsize)
{
  int i,j,jn;
  int size2=mboxsize*nboxsize;
  int ind=0;
  ifstream infile;
  char filename[300];
  sprintf (filename, "%sInds_%d.txt",dirstr, resumetimestep);
  infile.open (filename);
  
  if (!infile.is_open()) 
  {
    cout << "Could not open file!";
    // exit();
  }
  
  string aline;
  j=0;
  while (getline( infile, aline))
  {
    stringstream ss( aline );
    jn=j*mboxsize;
    for (i=0;i<mboxsize;i++)
    {
      ss>> ind;
      inds[i+jn]=ind;
      eta[i+jn]=1.00;
    }
    j++;
  }
  //count number of grains: in an easy way number of grains can not be larger than largest number in inds
  int maxind;
  maxind=inds[0];
  for (i=0;i<mboxsize*nboxsize;i++)
  {
    if (maxind>inds[i]) maxind=inds[i];
  }
  
  return maxind;
}
