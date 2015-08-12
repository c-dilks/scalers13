// fill-by-fill or pattern by pattering averaging for S_LL distributions

void AverageSLL(const char * filename="fit_result.zdcx.vpdx.root")
{
  const Int_t NBINS = 100;
  const Double_t BOUND = 6e-3;

  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("restr");

  Int_t i;
  Int_t fi;
  Int_t aa;
  Int_t pattern;
  Double_t cons;
  Double_t cons_err;
  Double_t epsi;
  Double_t epsi_err;
  Double_t asym;
  Double_t asym_err;
  tr->SetBranchAddress("i",&i);
  tr->SetBranchAddress("fi",&fi);
  tr->SetBranchAddress("aa",&aa);
  tr->SetBranchAddress("pattern",&pattern);
  tr->SetBranchAddress("cons",&cons);
  tr->SetBranchAddress("cons_err",&cons_err);
  tr->SetBranchAddress("epsi",&epsi);
  tr->SetBranchAddress("epsi_err",&epsi_err);
  tr->SetBranchAddress("asym",&asym);
  tr->SetBranchAddress("asym_err",&asym_err);

  Int_t pat[8] = {13,14,23,24,31,32,41,42};

  //+++++++++++++++++++++++++++++++++++++++++++

  // one entry = scaler asym for one run:
  TH2D * dist_by_run[3][8]; // [aa] [pattern]
  TString dist_by_run_n[3][8];
  TString dist_by_run_t[3][8];

  // one entry = scaler asym averaged over a fill:
  TH2D * dist_by_fill[3][8];
  TString dist_by_fill_n[3][8];
  TString dist_by_fill_t[3][8];
  
  Int_t aa,pp;
  for(aa=0; aa<3; aa++)
  {
    for(pp=0; pp<8; pp++)
    {
      dist_by_run_n[aa][pp] = Form("dist_by_run_a%d_pat%d",aa+1,pat[pp]);
      dist_by_fill_n[aa][pp] = Form("dist_by_fill_a%d_pat%d",aa+1,pat[pp]);
      dist_by_run_t[aa][pp] = Form("S_{%d} by run for pattern %d",aa+1,pat[pp]);
      dist_by_fill_t[aa][pp] = Form("S_{%d} by fill for pattern %d",aa+1,pat[pp]);



