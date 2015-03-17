void chan_tree(Int_t runnum=14131025)
{
  TFile * outfile = new TFile("chtr.root","RECREATE");
  char filename[64];
  sprintf(filename,"chanfiles/run%d_4.dat",runnum);
  TTree * tr = new TTree();
  tr->ReadFile(filename,"bx/I:chBBC/I:chZDC/I:chVPD/I:channel/D:count/D");
  tr->Write("tr");
};
