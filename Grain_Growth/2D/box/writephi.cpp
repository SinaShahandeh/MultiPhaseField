void writephi(*phi,tn)
{
    // write array into a file each 100 time steps
      ofstream myfile;
      // make a string like "result_5.txt"
      int n;
      char filename[200];
      n=sprintf (filename, "results/result_ %d .txt", tn);
      myfile.open (filename);
      for (i=0;i<mboxsize;i++)
      {
        for (j=0;j<nboxsize;j++)
        {
          myfile << phi[i+j*mboxsize] << "      "; 
        }
        myfile << "\n";
      }
      myfile.close();
}