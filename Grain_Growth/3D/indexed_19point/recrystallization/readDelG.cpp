
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int readDelG(char* dirstr, double* DelG, int mboxsize, int nboxsize, int kboxsize)
{
  int i,j,k,jn;
  int size2=mboxsize*nboxsize;
  double Gval=0.0;
  ifstream infile;
  char filename[300];
  sprintf (filename, "%sDelG.txt",dirstr);
  infile.open (filename);
  
  if (!infile.is_open()) 
  {
    cout << "Could not open DelG.txt file! "  << endl;
    exit(1);
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
      ss>> Gval;
      DelG[i+jn+k*size2]=Gval*2;
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
