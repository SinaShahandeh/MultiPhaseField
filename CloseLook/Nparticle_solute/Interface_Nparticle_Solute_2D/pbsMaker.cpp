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
  int var1[3]={1, 2, 3};
  int var2[3]={10, 5, 0};
  
  for (int ivar1=0;ivar1<3;ivar1++){
    for (int ivar2=0;ivar2<3;ivar2++){
      
      char filename[200];
      ofstream pbsfile;
      sprintf (filename, "%s_a%d_fv%d.pbs",inputstr, var1[ivar1],var2[ivar2]);
      pbsfile.open (filename);
      //writing
      pbsfile << "#PBS -l walltime=0:10:00:00" <<endl;
      pbsfile << "#PBS -l mem=5GB " << endl;
      pbsfile << "#PBS -l nodes=1:ppn=1 " << endl;
      pbsfile << " " << endl;
      pbsfile << "#PBS -m bea " << endl;
      pbsfile << "#PBS -j oe " << endl;
      pbsfile << "#PBS -N " << inputstr << var1[ivar1] << "_DelG" << var2[ivar2] << endl;
      pbsfile << " " << endl;
      pbsfile << "cd $PBS_O_WORKDIR " << endl;
      pbsfile << "CURR_DIR=`pwd` " << endl;
      pbsfile << "USERID=`whoami` " << endl;
      pbsfile << "JOB_OUTPUT=RESULTS " << endl;
      pbsfile << "echo 'Current working directory is \"'$CURR_DIR'\"' " << endl;
      pbsfile << "echo \"Running on `hostname`\" " << endl;
      pbsfile << "echo \"Starting run at: `date`\" " << endl;
      
      pbsfile << "var1=" << var1[ivar1] << endl;
      pbsfile << "var2=" << var2[ivar2] << endl;
      pbsfile << "jobname=\"" << inputstr << "_fv${var1}_DelG${var2}\" " << endl;
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
      pbsfile << "./Nparticle.out \"$jobname\" $var1 $var2 " << endl;
      pbsfile.close();
      
      cout << "File for simulation " << filename << " was writtten." <<endl;
      
      char runcommand[200];
      sprintf(runcommand,"qsub %s", filename);
      cout << runcommand << endl;
      system( runcommand);
    }
  }
}


