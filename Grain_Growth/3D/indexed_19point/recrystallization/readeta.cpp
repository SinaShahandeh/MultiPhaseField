
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int readeta(char* dirstr, int* inds, double* eta, int resumetimestep, int mboxsize, int nboxsize, int kboxsize, int p)
{
  int i,j,k,jn,pn,pnn;
  int size2=mboxsize*nboxsize;
  int size3=mboxsize*nboxsize*kboxsize;
  int indijk=0;
  double etaijk=0.0;
  ifstream infile;
  char filename[300];
  char filenamezip[300];
  char command[300];
  bool istxtexist=false;
  //reading inds
  for (pn=0;pn<p;pn++)
  {
    pnn=pn*size3;
    sprintf (filename, "%sInds_%d_%d.txt",dirstr, pn, resumetimestep);
    infile.open (filename);
    istxtexist=true;
    if (!infile.is_open())
    {
     istxtexist=false;
      cout << "Could not open txt file to resume for time step " <<resumetimestep << " and p= " << pn << ". I try to unzip." << endl;
      sprintf (filenamezip, "%sInds_%d_%d.zip",dirstr, pn, resumetimestep);
      sprintf (command, "unzip %s -d %s", filenamezip, dirstr);
      int n=system(command);
      infile.open (filename);
      if (!infile.is_open())
      {
        cout << "Even could not open zip file! time step= " <<resumetimestep << " and p= " << pn << endl;
        exit(1);
      }
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
        ss>> indijk;
        inds[i+jn+k*size2+pnn]=indijk;
      }
      j++;
      if (j%nboxsize==0)
      {
        j=0;
        k++;
      }
    }
    infile.close();
    if (istxtexist==false){
	sprintf (command, "rm %s", filename);
	cout << command << endl;
        int n=system(command);
    }
  }
  // reading etas
  for (pn=0;pn<p;pn++)
  {
     pnn=pn*size3;
    sprintf (filename, "%sEta_%d_%d.txt",dirstr, pn, resumetimestep);
    infile.open (filename);
    istxtexist=true;
    if (!infile.is_open())
    {
      istxtexist=false;
      cout << "Could not open txt file to resume for time step " <<resumetimestep << " and p= " << pn << ". I try to unzip." << endl;
      sprintf (filenamezip, "%sEta_%d_%d.zip",dirstr, pn, resumetimestep);
      sprintf (command, "unzip %s -d %s", filenamezip, dirstr);
      int n=system(command) ;
      infile.open (filename);
      if (!infile.is_open())
      {
        cout << "Even could not open zip file! time step= " <<resumetimestep << " and p= " << pn << endl;
        exit(1);
      }
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
        ss>> etaijk;
        eta[i+jn+k*size2+pnn]=etaijk;
      }
      j++;
      if (j%nboxsize==0)
      {
        j=0;
        k++;
      }
    }
    infile.close();
    if (istxtexist==false){
	sprintf (command, "rm %s", filename);
	cout << command << endl;
        int n=system(command);
    }

  }
  return 0;
}
