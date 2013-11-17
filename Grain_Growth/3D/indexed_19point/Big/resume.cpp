
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int resume(char* dirstr, int* inds, double* eta, int resumetimestep, int mboxsize, int nboxsize, int kboxsize)
{
  int i,j,k,jn;
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
  k=0;
  j=0;
  while (getline( infile, aline))
  {
    stringstream ss( aline );
    jn=j*mboxsize;
    for (i=0;i<mboxsize;i++)
    {
      ss>> ind;
      inds[i+jn+k*size2]=ind;
      eta[i+jn+k*size2]=1.00;
      if (ind==0)
        eta[i+jn+k*size2]=0.0;
    }
    j++;
    if (j%nboxsize==0)
    {
      j=0;
      k++;
    }
  }
  return 0;
}
