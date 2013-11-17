#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <sstream>
using namespace std;
int main(int a, char** charinput){
  //char* inputstr;
  //inputstr=charinput[1];
  cout << "What is the name of the simulation series input file?" << endl;
  char inputstr[200];
  cin >> inputstr;
    // inputstr="Np300_m1_k2_r5"
  int fvsteps=1;
  int DelGsteps=4; 
  int fv[1]={0}; //divided by 1000
  int DelG[4]={200, 150, 100, 50 };  //divided by 1000
  
  for (int ifv=0;ifv<fvsteps;ifv++){
    for (int iG=0;iG<DelGsteps;iG++){
      
      char filename[200];
      ofstream pbsfile;
      sprintf (filename, "%s_fv%d_DelG%d.pbs",inputstr, fv[ifv],DelG[iG]);
      pbsfile.open (filename);
      //writing
      pbsfile << "#PBS -l walltime=1:00:00:00" <<endl;
      pbsfile << "#PBS -l mem=5GB " << endl;
      pbsfile << "#PBS -l nodes=1:ppn=12 " << endl;
      pbsfile << " " << endl;
      pbsfile << "#PBS -m bea " << endl;
      pbsfile << "#PBS -j oe " << endl;
      pbsfile << "#PBS -N Np300_m2_k4_r5_fv" << fv[ifv] << "_DelG" << DelG[iG] << endl;
      pbsfile << " " << endl;
      pbsfile << "cd $PBS_O_WORKDIR " << endl;
      pbsfile << "CURR_DIR=`pwd` " << endl;
      pbsfile << "USERID=`whoami` " << endl;
      pbsfile << "JOB_OUTPUT=RESULTS " << endl;
      pbsfile << "echo 'Current working directory is \"'$CURR_DIR'\"' " << endl;
      pbsfile << "echo \"Running on `hostname`\" " << endl;
      pbsfile << "echo \"Starting run at: `date`\" " << endl;
      
      pbsfile << "fv=" << fv[ifv] << endl;
      pbsfile << "DelG=" << DelG[iG] << endl;
      pbsfile << "jobname=\"" << inputstr << "_fv${fv}_DelG${DelG}\" " << endl;
      pbsfile << "echo $jobname " << endl;
      pbsfile << "echo \"---------------------\" " << endl;
      pbsfile << "dirscratch=\"/global/scratch/ssina/Np/${jobname}/\" " << endl;
      pbsfile << "mkdir $dirscratch " << endl;
      pbsfile << "cd $dirscratch" << endl;
      pbsfile << " echo `pwd`" << endl;
      pbsfile << "cp $CURR_DIR/Nparticle.out $dirscratch " << endl;
      pbsfile << "cp $CURR_DIR/main.cpp  $dirscratch" << endl;
      pbsfile << "cp $CURR_DIR/" << inputstr << " $dirscratch" << endl;
      pbsfile << "mv "<< inputstr << " $jobname" << endl;
      pbsfile << "./Nparticle.out \"$jobname\" $DelG $fv " << endl;
      pbsfile.close();
      
      cout << "File for simulation " << filename << " was writtten." <<endl;
      
      char runcommand[200];
      sprintf(runcommand,"qsub %s", filename);
      cout << runcommand << endl;
      system( runcommand);
    }
  }
}


