// box class keeps the internal array of order parameters
using namespace std;

class box
{
  double *eta;
  double *eta2;
  double *phii;
  double *sumtermi;

  bool *mbool;
  int boxmargin; //box is larger in each direction than calculatable area with the margin value
  public:
    int mboxsize;
    int nboxsize;
    int upleft[2];
    int downright[2];

    double *getsquare ();
    void *doatimestep();
    void *setboxcordinate();
    void *makeanucli;
    bool isempty;
    //constructor
    box(){
    }
    //deconstructor
    ~box(){
      delete *eta;
      delete *eta2;
    }
}

void box::nucleate(int uplefti, int upleftj, int size){
  mboxsize=size; //width and height of the box
  nboxsize=size;
  boxmargin=11; 
  upleft[0]=uplefti;
  upleft[1]=upleftj;
  downright[0]=uplefti+size;
  downright[1]=upleftj+size;
  eta= new double[mboxsize*nboxsize];
  eta2= new double[mboxsize*nboxsize];
  mbool= new bool[mboxsize*nboxsize];
  eta[(mboxsize*nboxize+1)/2]=1.0;
  mbool(mboxsize*nboxize+1)/2]=true;
}

void *box::getsquare(){
  int i, j, jn;
  phii= new double[mboxsize*nboxsize];
  for (i=0;i<mboxsize;i++)
  {
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      phii[i+jn]=phii[i+jn]+eta[i+jn]*eta[i+jn];
      if (phi[i+jn]>0.0001 & phi[i+jn]<0.9999) // optimization function
      {mbool[i+jn]=true;
      // @@@@@ Also add neighbors to be true
      }
      else
      {mbool[i+jn]=false;}
    }
  }
  return phii;
}

void box::deletephi(){
  delete phi;
}

bool box::isempty(){
  if (mboxsize<boxmargin){if (nboxsize<boxmargin){return true;}}
  return false;
}
void box::adapt(){
  int i,j,jn;
  int loweri, upperi, lowerj, upperj;
  int boxmargin0; //current margin araund the order parameter.
  loweri=mboxsize;
  lowerj=mboxsize;
  upperi=0;
  upperj=0;
  //finding what is the upper and lower bounding indecies for mbool which is our calculation domain.
  for (i=0;i<mboxsize;i++)
  {
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      if (mbool[i+jn]=true){
        if (i<loweri){loweri=i}
        if (i>upperi){upperi=i}
        if (j<lowerj){lowerj=j}
        if (j>upperj){upperj=j}
      }
    }
  }
  mboxsize=((loweri-upperi)+2*boxmargin); //new size of the box
  nboxsize=((lowerj-upperj)+2*boxmargin);
  double* neweta;
  neweta= new double[mboxsize*nboxsize];
  int leftshift=upperi-boxmargin;
  int upshift=upperj-boxmargin;
  for (i=1;i<mboxsize;i++)
  {
    for (j=0;j<nboxsize;j++)
    {
      if ((leftshift+i)>=0 && (upshift+j)>=0){
        jn=(j+upshift)*mboxsize;
        neweta[i+jn]=eta[(leftshift+i)*(jn)];
      }
    }
  }
}
// Solves evolution equation and finds order parameter at next time step. 
void box::doatimestep{
  int i, j , in;
  double del2, sumtermp,detadtM, detadt;
  double delx2=1.0/(delx*delx);
  for (i=1;i<mboxsize;i++)
  {
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      // the sum of all order parameters^2 for the point i and j is calculated in the main function and was sent to the class before using the "sumterm".
      // calculation of nabla square eta
      del2[pn]=delx2*((eta[peribc(i+1,nboxsize)+jn+pnn]+eta[peribc(i-1,nboxsize)+jn+pnn]+eta[i+peribc(j+1,mboxsize)*mboxsize+pnn]+eta[i+peribc(j-1,mboxsize)*mboxsize+pnn])
      +0.25*(eta[peribc(i+1,nboxsize)+peribc(j+1,mboxsize)*mboxsize+pnn]+eta[peribc(i+1,nboxsize)+peribc(j-1,mboxsize)*mboxsize+pnn]+eta[peribc(i-1,nboxsize)+peribc(j+1,mboxsize)*mboxsize+pnn]+eta[peribc(i-1,nboxsize)+peribc(j-1,mboxsize)*mboxsize+pnn])
      -5*eta[i+jn+pnn]);
      sumtermp=eta[i+jn]*sumtermi-pow(eta[i+jn],3);
      detadtM=-alpha*eta[i+jn]+beta*pow(eta[i+jn],3)-kappa*del2[pn];
      detadt=-L*(detadtM[pn]+2*gamma*sumtermp);
      eta2[i+jn]=eta[i+jn]+delt*detadt[pn];
      // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
      if (eta2[i+jn]>1) eta2[i+jn]=1;
      if (eta2[i+jn]<0) eta2[i+jn]=0;
    }
  }
  //setting eta equal to the new eta2 for the next time step
  for (i=0;i<mboxsize;i++)
  {
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      for (pind=0;pind<p;pind++)
      {
        pnn=pind*size;
        eta[i+jn]=eta2[i+jn];
      }
    }
  }
}