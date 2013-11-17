// this function calculates energy of the system into a matrix and then sum it up
using namespace std;
double energy(double* eta,int mboxsize, int nboxsize, int p, double delx, double m, double kappa, double* ME)
{
  int i, j, jn, pn, pnn, jnplus1, jnminus1, inplus1, inminus1;
  int size2=mboxsize*nboxsize;
  double grad2, sumij,f0;
  double gg;
  double E;
  for (i=0;i<size2;i++)
    ME[i]=0.0;
  
  for (j=0;j<nboxsize;j++)
  {
    jn=j*mboxsize;
    jnplus1=symbc(j+1,nboxsize)*mboxsize;
    jnminus1=symbc(j-1,nboxsize)*mboxsize;
    for (i=0;i<mboxsize;i++)
    {
      inplus1=symbc(i+1,mboxsize);
      inminus1=symbc(i-1,mboxsize);
      sumij=0.0;
      f0=0.0;
      grad2=0.0;
      for (pn=0;pn<p;pn++)
      {
	pnn=pn*size2;
	for (int pnj=pn+1;pnj<p;pnj++)
	{
	  sumij=sumij+eta[i+jn+pnj*size2]*eta[i+jn+pnj*size2]*eta[i+jn+pnn]*eta[i+jn+pnn];
	}
	f0=f0+(-eta[i+jn+pnn]*eta[i+jn+pnn]*0.50+ eta[i+jn+pnn]*eta[i+jn+pnn]*eta[i+jn+pnn]*eta[i+jn+pnn]*0.250);
	gg=(eta[inplus1+jn+pnn]-eta[inminus1+jn+pnn])*(eta[inplus1+jn+pnn]-eta[inminus1+jn+pnn])+(eta[i+jnplus1+pnn]-eta[i+jnminus1+pnn])*(eta[i+jnplus1+pnn]-eta[i+jnminus1+pnn]);
	grad2=grad2+0.250/delx/delx*gg;
      }
      f0=f0+1.50*sumij;
      ME[i+jn]=m*(0.250+f0)+ kappa/2*grad2;
    }
  }
  
  // total energy of the system
  E=0.0;
  for (i=0;i<size2;i++)
    E=E+ME[i];
  
  return E;
}
