using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace CMPE_2D_GG
{
    public partial class Form1 : Form
    {


        //getting an index of the domain and returning outindex considering periodic boundary condition
        static public int peribc(int inindex, int size)
        {
            int outindex;
            outindex = inindex;
            if (inindex < 0)
            {
                outindex = inindex + size;
            }
            if (inindex >= size)
            {
                outindex = inindex - size;
            }
            return (outindex);
        }

        static public double sign(double x)
        {
            if (x > 0) { return 1; }
            if (x < 0) { return -1; }
            if (x == 0) { return 0; }
            return 0;
        }



        string strSave;
        public Form1()
        {
            InitializeComponent();
        }

        private void btnSave_Click(object sender, EventArgs e)
        {
            folderBrowserDialog1.ShowDialog();
            strSave = folderBrowserDialog1.SelectedPath;
            lblSave.Text = strSave;
        }

        int inplus1, inminus1, jnplus1, jnminus1;
        private void btnStart_Click(object sender, EventArgs e)
        {

            bool fullwrite = true;
            // model parameters
            double delt = Convert.ToDouble(txtDelt.Text);
            int timesteps1 = Convert.ToInt32(txtt0.Text);
            int timesteps = Convert.ToInt32(txttimesteps.Text);
            int writingtimesteps = 10;
            int fulletatimesteps = Convert.ToInt32(txtsaveres.Text);
            int writegrainstat = Convert.ToInt32(txtSaveGsize.Text);
            int resumetimestep = 0;
            double L = Convert.ToDouble(txtL.Text);
            double m = Convert.ToDouble(txtm.Text); ;
            double kappa = Convert.ToDouble(txtkappa.Text);
            double gamma = 2 * 1.5 * m;
            double Pz;
            Pz = Convert.ToDouble(txtPz.Text);
            double Lf = L, pzi;
            int i, j, tn;
            // geometry settings
            int mboxsize = Convert.ToInt32(txtmboxsize.Text); // x axis in pixels
            int nboxsize = Convert.ToInt32(txtnboxsize.Text); // y axis
            double delx = Convert.ToDouble(txtDelx.Text);      // length unit per pixel
            int p = 5; // depth of indexed matrix. 6 should be more than enough for 3D system since quadrouple points have 4 order parameters in common.
            double thresh = 0.00000001; //threshold value for choosing active nodes
            double[] eta = new double[mboxsize * nboxsize * p];
            double[] eta2 = new double[mboxsize * nboxsize * p];
            int[] inds = new int[mboxsize * nboxsize * p];
            bool[] mbool = new bool[mboxsize * nboxsize];
            for (i = 0; i < mboxsize * nboxsize * p; i++)
            {
                eta[i] = 0.00;
                eta2[i] = 0;
            }
            for (i = 0; i < mboxsize * nboxsize; i++)
            {
                mbool[i] = true;
            }
            double phii;
            double[] phi;
            phi = new double[mboxsize * nboxsize];
            int[] Maxindfield;
            Maxindfield = new int[mboxsize * nboxsize];
            double meanG;

            //visualiation objects declarations
            var img = new Bitmap(mboxsize, nboxsize);

            // number of nucleas at the beginning of simulation
            int nuclein;
            nuclein = (int)(mboxsize * nboxsize / 200); // ~1 percent of grid points are nuclei
            int nn, ii, jj, pp, pos;

            int[] MGsize = new int[nuclein];
            Random rnd = new Random();
            if (resumetimestep == 0)
            {
                for (nn = 1; nn < nuclein + 1; nn++)
                {
                    ii = Convert.ToInt32((nboxsize * rnd.NextDouble()));
                    jj = (int)((mboxsize * rnd.NextDouble()));
                    pp = (int)(p * rnd.NextDouble()) * mboxsize * nboxsize;
                    eta[ii + jj * mboxsize + pp] = 1;
                    inds[ii + jj * mboxsize + pp] = nn;
                    inds[peribc(ii + 1, mboxsize) + jj * mboxsize + pp] = nn;
                    inds[peribc(ii - 1, mboxsize) + jj * mboxsize + pp] = nn;
                    inds[ii + peribc(jj + 1, nboxsize) * mboxsize + pp] = nn;
                    inds[ii + peribc(jj - 1, nboxsize) * mboxsize + pp] = nn;
                }
            }
            //if (resumetimestep>0){
            // readeta(dirstr, inds, eta, resumetimestep, mboxsize, nboxsize, 1,  p);
            //}
            // particles distribution specification
            /* double diameter=2;
             *  double particles_fraction=0.00;
             *  double particlesn=particles_fraction*mboxsize*nboxsize/diameter^2   //
             *  particles number
             *  double ppf[nboxsize][mboxsize]
             *  here goes the function to make particle distribution 
             */
            //dynamics
            double sumterm, currenteta;
            double detadtM;
            double detadt;
            int pn, psn, pind;
            double delx2 = 2.0000000 / 3.0000000 / (delx * delx);
            int size2 = mboxsize * nboxsize;
            int jn, pnn;
            double del2;
            double sumeta, mineta, sumeta2;
            int currentind, indc, minind, indscn;
            int maxind;
            double maxeta;
            for (tn = resumetimestep; tn < timesteps1; tn++)
            {
                //     #pragma omp parallel for
                for (j = 0; j < nboxsize; j++)
                {

                    jn = j * mboxsize;
                    jnplus1 = Form1.peribc(j + 1, nboxsize) * mboxsize;
                    jnminus1 = peribc(j - 1, nboxsize) * mboxsize;
                    for (i = 0; i < mboxsize; i++)
                    {
                        inplus1 = peribc(i + 1, mboxsize);
                        inminus1 = peribc(i - 1, mboxsize);
                        // here is the sum of all order parameters^2 for the point i and j
                        sumterm = 0;
                        for (psn = 0; psn < p; psn++)
                        {
                            sumterm = sumterm + eta[i + jn + psn * size2] * eta[i + jn + psn * size2];
                        }
                        // calculation of nabla square eta
                        for (pn = 0; pn < p; pn++)
                        {
                            pnn = pn * size2;
                            mbool[i + jn] = false; //firts we make is false and then if detadt>thresh then it becomes true later
                            //
                            currentind = inds[i + jn + pnn];
                            currenteta = eta[i + jn + pnn];
                            //searching for neighbors with the same index as currentind
                            sumeta = 0;
                            for (indc = 0; indc < p; indc++)
                            {
                                indscn = indc * size2;
                                if (currentind == inds[inplus1 + jn + indscn]) { sumeta = sumeta + eta[inplus1 + jn + indc * size2]; }
                                if (currentind == inds[inminus1 + jn + indscn]) { sumeta = sumeta + eta[inminus1 + jn + indc * size2]; }
                                if (currentind == inds[i + jnplus1 + indscn]) { sumeta = sumeta + eta[i + jnplus1 + indc * size2]; }
                                if (currentind == inds[i + jnminus1 + indscn]) { sumeta = sumeta + eta[i + jnminus1 + indc * size2]; }
                            }
                            sumeta2 = 0;
                            for (indc = 0; indc < p; indc++)
                            {
                                indscn = indc * size2;
                                if (currentind == inds[inplus1 + jnplus1 + indscn]) { sumeta2 = sumeta2 + eta[inplus1 + jnplus1 + indc * size2]; }
                                if (currentind == inds[inminus1 + jnminus1 + indscn]) { sumeta2 = sumeta2 + eta[inminus1 + jnminus1 + indc * size2]; }
                                if (currentind == inds[inminus1 + jnplus1 + indscn]) { sumeta2 = sumeta2 + eta[inminus1 + jnplus1 + indc * size2]; }
                                if (currentind == inds[inplus1 + jnminus1 + indscn]) { sumeta2 = sumeta2 + eta[inplus1 + jnminus1 + indc * size2]; }
                            }
                            del2 = delx2 * (sumeta + 0.25 * sumeta2 - 5 * currenteta);
                            detadtM = m * (currenteta * currenteta * currenteta - currenteta) - kappa * del2;
                            detadt = -L * (detadtM + gamma * (currenteta * sumterm - currenteta * currenteta * currenteta));
                            eta2[i + jn + pnn] = currenteta + delt * detadt;

                            if (eta2[i + jn + pnn] > thresh) // make a padding for the changing order parameter so next time step it will take that into account.
                            {
                                mineta = eta[inplus1 + jn];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[inplus1 + jn + indc * size2]) { mineta = eta[inplus1 + jn + indc * size2]; minind = indc; }
                                    if (inds[inplus1 + jn + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[inplus1 + jn + minind * size2] = currentind;
                                mineta = eta[inminus1 + jn];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[inminus1 + jn + indc * size2]) { mineta = eta[inminus1 + jn + indc * size2]; minind = indc; }
                                    if (inds[inminus1 + jn + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[inminus1 + jn + minind * size2] = currentind;
                                mineta = eta[i + jnplus1];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[i + jnplus1 + indc * size2]) { mineta = eta[i + jnplus1 + indc * size2]; minind = indc; }
                                    if (inds[i + jnplus1 + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[i + jnplus1 + minind * size2] = currentind;
                                mineta = eta[i + jnminus1];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[i + jnminus1 + indc * size2]) { mineta = eta[i + jnminus1 + indc * size2]; minind = indc; }
                                    if (inds[i + jnminus1 + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[i + jnminus1 + minind * size2] = currentind;
                                // ------------------- second neighbors --------------------------
                                // inminus1+jnminus1 -1 -1
                                pos = inminus1 + jnminus1;
                                mineta = eta[pos];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                    if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[pos + minind * size2] = currentind;
                                // inplus1+jnminus1 1 -1
                                pos = inplus1 + jnminus1;
                                mineta = eta[pos];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                    if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[pos + minind * size2] = currentind;
                                pos = inminus1 + jnplus1; // -1 1
                                mineta = eta[pos];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                    if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[pos + minind * size2] = currentind;
                                pos = inplus1 + jnplus1; // 1 1
                                mineta = eta[pos];  //the first in the table is assumed as smallest 
                                minind = 0;
                                for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                {
                                    if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                    if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                }
                                inds[pos + minind * size2] = currentind;
                            }
                        }
                        phii = 0;
                        for (pind = 0; pind < p; pind++)
                        {
                            pnn = pind * size2;
                            phii = phii + eta[i + jn + pnn] * eta[i + jn + pnn];
                        }
                        if (phii < 0.999999) //this point in at interface
                        {
                            mbool[i + jn] = true;
                            mbool[inplus1 + jn] = true;
                            mbool[inminus1 + jn] = true;
                            mbool[i + jnplus1] = true;
                            mbool[i + jnminus1] = true;
                            mbool[inplus1 + jnplus1] = true;
                            mbool[inplus1 + jnminus1] = true;
                            mbool[inminus1 + jnplus1] = true;
                            mbool[inminus1 + jnminus1] = true;
                        }
                    }
                }
                //setting eta equal to the new eta2 for the next time step
                // #pragma omp parallel for
                for (pind = 0; pind < p; pind++)
                {
                    pnn = pind * size2;
                    for (j = 0; j < nboxsize; j++)
                    {
                        jn = j * mboxsize;
                        for (i = 0; i < mboxsize; i++)
                        {
                            eta[i + jn + pnn] = eta2[i + jn + pnn];
                        }
                    }
                }

            }

            //optimized loop -----------------------------------------------------------------------------------------------------
            double pastsize, newsize;
            pastsize = 100;
            newsize = pastsize + 2;
            while (tn < timesteps) //newsize-pastsize>0.0000000000001
            {
                lbltn.Text = Convert.ToString(tn);
                lbltn.Refresh();
                //#pragma omp parallel for
                for (j = 0; j < nboxsize; j++)
                {
                    jn = j * mboxsize;
                    jnplus1 = peribc(j + 1, mboxsize) * mboxsize;
                    jnminus1 = peribc(j - 1, mboxsize) * mboxsize;
                    for (i = 0; i < mboxsize; i++)
                    {
                        if (mbool[i + jn] == true)//mbool[i+jn+pnn]==true
                        {
                            inplus1 = peribc(i + 1, mboxsize);
                            inminus1 = peribc(i - 1, mboxsize);
                            mbool[i + jn] = false; //firts we make is false and then if detadt>thresh then it becomes true later
                            for (pn = 0; pn < p; pn++)
                            {
                                pnn = pn * size2;
                                // here is the sum of all order parameters^2 for the point i and j
                                sumterm = 0;
                                for (psn = 0; psn < p; psn++)
                                {
                                    sumterm = sumterm + eta[i + jn + psn * size2] * eta[i + jn + psn * size2];
                                }
                                currentind = inds[i + jn + pnn];
                                currenteta = eta[i + jn + pnn];
                                //searching for neighbors with the same index as currentind
                                sumeta = 0;
                                for (indc = 0; indc < p; indc++)
                                {
                                    indscn = indc * size2;
                                    if (currentind == inds[inplus1 + jn + indscn]) { sumeta = sumeta + eta[inplus1 + jn + indc * size2]; }
                                    if (currentind == inds[inminus1 + jn + indscn]) { sumeta = sumeta + eta[inminus1 + jn + indc * size2]; }
                                    if (currentind == inds[i + jnplus1 + indscn]) { sumeta = sumeta + eta[i + jnplus1 + indc * size2]; }
                                    if (currentind == inds[i + jnminus1 + indscn]) { sumeta = sumeta + eta[i + jnminus1 + indc * size2]; }
                                }
                                sumeta2 = 0;
                                for (indc = 0; indc < p; indc++)
                                {
                                    indscn = indc * size2;
                                    if (currentind == inds[inplus1 + jnplus1 + indscn]) { sumeta2 = sumeta2 + eta[inplus1 + jnplus1 + indc * size2]; }
                                    if (currentind == inds[inminus1 + jnminus1 + indscn]) { sumeta2 = sumeta2 + eta[inminus1 + jnminus1 + indc * size2]; }
                                    if (currentind == inds[inminus1 + jnplus1 + indscn]) { sumeta2 = sumeta2 + eta[inminus1 + jnplus1 + indc * size2]; }
                                    if (currentind == inds[inplus1 + jnminus1 + indscn]) { sumeta2 = sumeta2 + eta[inplus1 + jnminus1 + indc * size2]; }
                                }
                                del2 = delx2 * (sumeta + 0.25 * sumeta2 - 5 * currenteta);
                                detadtM = m * (currenteta * currenteta * currenteta - currenteta) + gamma * (currenteta * sumterm - currenteta * currenteta * currenteta) - kappa * del2;
                                pzi = 3 * currenteta * (1 - currenteta) * sign(detadtM) * Pz;
                                if (Math.Abs(detadtM) < Math.Abs(pzi))
                                {
                                    Lf = 0;
                                }
                                else
                                {
                                    Lf = L;
                                }
                                detadt = -Lf * (detadtM - pzi);
                                eta2[i + jn + pnn] = currenteta + delt * detadt;
                                // to make sure eta is not outside the equilibrium values. This increases stability of calculation by controlling bounds of the eta whithin equilibrium values
                                if (eta2[i + jn + pnn] > 1) eta2[i + jn + pnn] = 1;
                                if (eta2[i + jn + pnn] < 0) eta2[i + jn + pnn] = 0;
                                // creating phi to find interface area from it
                                if (eta2[i + jn + pnn] > thresh) // make a padding for the changing order parameter so next time step it will take that into account.
                                {
                                    mineta = eta[inplus1 + jn];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[inplus1 + jn + indc * size2]) { mineta = eta[inplus1 + jn + indc * size2]; minind = indc; }
                                        if (inds[inplus1 + jn + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[inplus1 + jn + minind * size2] = currentind;

                                    mineta = eta[inminus1 + jn];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[inminus1 + jn + indc * size2]) { mineta = eta[inminus1 + jn + indc * size2]; minind = indc; }
                                        if (inds[inminus1 + jn + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[inminus1 + jn + minind * size2] = currentind;

                                    mineta = eta[i + jnplus1];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[i + jnplus1 + indc * size2]) { mineta = eta[i + jnplus1 + indc * size2]; minind = indc; }
                                        if (inds[i + jnplus1 + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[i + jnplus1 + minind * size2] = currentind;

                                    mineta = eta[i + jnminus1];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[i + jnminus1 + indc * size2]) { mineta = eta[i + jnminus1 + indc * size2]; minind = indc; }
                                        if (inds[i + jnminus1 + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[i + jnminus1 + minind * size2] = currentind;

                                    // ------------------- second neighbors --------------------------
                                    // inminus1+jnminus1 -1 -1
                                    pos = inminus1 + jnminus1;
                                    mineta = eta[pos];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                        if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[pos + minind * size2] = currentind;

                                    // inplus1+jnminus1 1 -1
                                    pos = inplus1 + jnminus1;
                                    mineta = eta[pos];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                        if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[pos + minind * size2] = currentind;

                                    pos = inminus1 + jnplus1; // -1 1
                                    mineta = eta[pos];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                        if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[pos + minind * size2] = currentind;

                                    pos = inplus1 + jnplus1; // 1 1
                                    mineta = eta[pos];  //the first in the table is assumed as smallest 
                                    minind = 0;
                                    for (indc = 0; indc < p; indc++) //loop around the neighbor and kick the minimum out and replace it with new growing eta
                                    {
                                        if (mineta > eta[pos + indc * size2]) { mineta = eta[pos + indc * size2]; minind = indc; }
                                        if (inds[pos + indc * size2] == currentind) { minind = indc; break; }
                                    }
                                    inds[pos + minind * size2] = currentind;
                                }
                            }
                        }
                        phii = 0;
                        for (pind = 0; pind < p; pind++)
                        {
                            pnn = pind * size2;
                            phii = phii + eta[i + jn + pnn] * eta[i + jn + pnn];
                        }
                        if (phii < 0.9999999) //this point in at interface
                        {
                            mbool[i + jn] = true;
                            mbool[inplus1 + jn] = true;
                            mbool[inminus1 + jn] = true;
                            mbool[i + jnplus1] = true;
                            mbool[i + jnminus1] = true;
                            mbool[inplus1 + jnplus1] = true;
                            mbool[inplus1 + jnminus1] = true;
                            mbool[inminus1 + jnplus1] = true;
                            mbool[inminus1 + jnminus1] = true;
                        }
                    }
                }
                //setting eta equal to the new eta2 for the next time step
                // #pragma omp parallel for
                for (pind = 0; pind < p; pind++)
                {
                    pnn = pind * size2;
                    for (j = 0; j < nboxsize; j++)
                    {
                        jn = j * mboxsize;
                        for (i = 0; i < mboxsize; i++)
                        {
                            eta[i + jn + pnn] = eta2[i + jn + pnn];
                        }
                    }
                }
                tn = tn + 1;

                // writing results
                if (tn % writegrainstat == 0)
                {
                    int rgb;
                    if (fullwrite == true)
                    {
                        // making the phi array
                        for (j = 0; j < nboxsize; j++)
                        {
                            jn = j * mboxsize;
                            for (i = 0; i < mboxsize; i++)
                            {
                                phi[i + jn] = 0;
                                for (pind = 0; pind < p; pind++)
                                {
                                    pnn = pind * size2;
                                    phi[i + jn] = phi[i + jn] + eta[i + jn + pnn] * eta[i + jn + pnn];
                                }
                                rgb = (int)(phi[i + jn] * 255);
                                img.SetPixel(i, j, Color.FromArgb(255, rgb, rgb, rgb));
                            }
                        }
                    }
                    pictureBox1.Image = img;
                    pictureBox1.Refresh();
                    lbltn.Text = Convert.ToString(tn);
                    lbltn.Refresh();
                    img.Save(strSave + "/" + Convert.ToString(tn) + ".png", System.Drawing.Imaging.ImageFormat.Png);
                }

            }


        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void label11_Click(object sender, EventArgs e)
        {

        }
    }

}
