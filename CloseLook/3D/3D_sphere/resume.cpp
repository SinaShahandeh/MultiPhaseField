
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int resume(char* dirstr, int* inds, double* eta, int resumetimestep, int mboxsize, int nboxsize, int kboxsize)
{
  int i,j,k,jn;
  int size2=mboxsize*nboxsize;
  
  ifstream infile;
  char filename[300];
  sprintf (filename, "%sInds_%d.txt",dirstr, resumetimestep);
  infile.open (filename);
  
  if (!infile.is_open()) 
  {
    cout << "Could not open file!";
  }

  string aline;
  while (getline( infile, aline))
  {
    stringstream ss( aline );
    string field;
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
        getline( ss, field, "   " );
        stringstream fs( field );
        int ind = 0; 
        fs >> ind;
        inds[i+jn]=ind;
        eta[i+jn+k*size2]=1.00;
      }
    }
    k=k+1;
  }
  return 0;
}