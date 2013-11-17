
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int readstructure(char* dirstr, int* inds, double* eta, int resumetimestep, int mboxsize, int nboxsize, int kboxsize, int p)
{
  int i,j,k,jn,pind,pn;
  int size2=mboxsize*nboxsize;
  int ind=0;
  double etai;
  ifstream infile;
  char filename[300];
  
  for (pind=0;pind<p;pind++)
  {
    pn=pind*mboxsize*nboxsize*kboxsize;
    sprintf (filename, "%sInds_%d_%d.txt",dirstr, pind,resumetimestep);
    infile.open (filename);
    
    if (!infile.is_open()) 
    {
      cout << "Could not open Inds file " << filename << endl;
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
        ss>> ind;
        inds[i+jn+k*size2+pn]=ind;
      }
      j++;
      if (j%nboxsize==0)
      {
        j=0;
        k++;
      }
    }
    infile.close();
  }

  for (pind=0;pind<p;pind++)
  {
    pn=pind*mboxsize*nboxsize*kboxsize;
    sprintf (filename, "%sEta_%d_%d.txt",dirstr, pind, resumetimestep);
    infile.open (filename);
    
    if (!infile.is_open()) 
    {
      cout << "Could not open Eta file " << filename << endl;
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
        ss>> etai;
        eta[i+jn+k*size2+pn]=etai;
      }
      j++;
      if (j%nboxsize==0)
      {
        j=0;
        k++;
      }
    }
    infile.close();
  }

  return 0;
}
