// find 10 largest eta from the neighbors and returns indecies of those important etas
int sorting(double* neighboreta,int* neighborinds, int* maxinds, int p)
{
  int i,ni, MaxIndexInNeighbors,mxc;
  double maxeta;
  int indsnum=0; //number of distinct order parameters
  bool alreadythere;
  for (ni=0;ni<10;ni++)
  {
    maxeta=neighboreta[0];
    MaxIndexInNeighbors=0;
    for (i=1;i<9*p;i++)
    {
      if (neighboreta[i]>maxeta)
      {
        maxeta=neighboreta[i];
        MaxIndexInNeighbors=i;
      }
    }
    //so we set vaule of maximum eta in the neighbor list to -1 which is smaller than all of the rest so next time second maximum wil be selected
    neighboreta[MaxIndexInNeighbors]=-1; 
    // add the max value to the maxinds array if it has not been added before
    alreadythere=false;
    for(mxc=0;mxc<indsnum;mxc++)
    {
      if (neighborinds[MaxIndexInNeighbors]==maxinds[mxc])
        alreadythere=true;
    }
    if (alreadythere==false)
    {
      maxinds[ni]=neighborinds[MaxIndexInNeighbors]; //output of the function
      indsnum=indsnum+1;
    }
  }
  //make remaining indecies in the maxinds to be 0.
  for (ni=indsnum;ni<10;ni++)
      maxinds[ni]=0;
  
  return indsnum;
}
