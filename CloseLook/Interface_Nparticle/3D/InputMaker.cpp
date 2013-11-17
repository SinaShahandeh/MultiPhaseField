#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <sstream>

int InputMaker(int a, char** charinput){
  double var;
  for (var=0.001;var<0.11;var=var+0.01;){
    
  char outstr[100];
  sprintf(outstr,"%s%d%s", "input_", int(var), ".txt")
  char* inputstr;
  inputstr=new char[50];
  sprintf(inputstr, "%s", "input_template.txt");
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
  cout << fp <<endl;
  double DelGinput;
  getline( infile, aline); getline( infile, aline); //47, 48
  cout << aline ;
  getline( infile, aline); // 55
  stringstream(aline)>> DelGinput;
  cout << DelGinput <<endl;

  
  }
}