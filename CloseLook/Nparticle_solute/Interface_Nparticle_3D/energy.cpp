// this function calculates energy of the system into a matrix and then sum it up
using namespace std;
double energy(double* eta,int mboxsize, int nboxsize, int kboxsize, int p, double delx, double m, double kappa, double* ME)
{
  int i, j, jn, k, kn, pn, pnn, knplus1,knminus1, jnplus1, jnminus1, inplus1, inminus1;
  int size2=mboxsize*nboxsize;
  int size3=mboxsize*nboxsize*kboxsize;
  double grad2, sumij,f0;
  double gg,gg1,gg2;
  double E;
  for (k=0;k<kboxsize;k++)
  {
    kn=k*size2;
    knplus1=peribc(k+1,kboxsize)*size2;
    knminus1=peribc(k-1,kboxsize)*size2;
    for (j=0;j<nboxsize;j++)
    {
      jn=j*mboxsize;
      jnplus1=peribc(j+1,mboxsize)*mboxsize;
      jnminus1=peribc(j-1,mboxsize)*mboxsize;
      for (i=0;i<mboxsize;i++)
      {
        inplus1=symbc(i+1,mboxsize);
        inminus1=symbc(i-1,mboxsize);
        sumij=0.0;
        f0=0.0;
        grad2=0.0;
	pnn=size3;
	sumij=eta[i+jn+kn+0*size3]*eta[i+jn+kn+0*size3]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn];
	gg1=1/2/delx*(eta[inplus1+jn+kn+pnn]-eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]-eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]-eta[i+jn+knminus1+pnn]);
	gg2=1/2/delx*(eta[inplus1+jn+kn+0]-eta[inminus1+jn+kn+0]+eta[i+jnplus1+kn+0]-eta[i+jnminus1+kn+0]+eta[i+jn+knplus1+0]-eta[i+jn+knminus1+0]);
	grad2=gg1*gg1+gg2*gg2;
	f0=(-eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*0.50+ eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*0.250)+
	(-eta[i+jn+kn+0]*eta[i+jn+kn+0]*0.50+ eta[i+jn+kn+0]*eta[i+jn+kn+0]*eta[i+jn+kn+0]*eta[i+jn+kn+0]*0.250);
	f0=f0+1.50*sumij;
	ME[i+jn+kn]=m*(0.250+f0)+ kappa/2*grad2;
	
//         for (pn=0;pn<p;pn++)
//         {
//           pnn=pn*size3;
//           for (int pnj=pn+1;pnj<p;pnj++)
//           {
//             sumij=sumij+eta[i+jn+kn+pnj*size3]*eta[i+jn+kn+pnj*size3]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn];
//           }
//           f0=f0+(-eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*0.50+ eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*eta[i+jn+kn+pnn]*0.250);
//           gg=eta[inplus1+jn+kn+pnn]-eta[inminus1+jn+kn+pnn]+eta[i+jnplus1+kn+pnn]-eta[i+jnminus1+kn+pnn]+eta[i+jn+knplus1+pnn]-eta[i+jn+knminus1+pnn];
//           grad2=grad2+0.250/delx/delx*gg*gg;
//         }
//         f0=f0+1.50*sumij;
//         ME[i+jn+kn]=m*(0.250+f0)+ kappa/2*grad2;
      }
    }
  }
  // total energy of the system
  E=0.0;
  for (i=0;i<size3;i++)
    E=E+ME[i];
  
  return E;
}
