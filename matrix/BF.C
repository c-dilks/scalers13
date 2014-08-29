// bunch fitting algorithm
// -- see equ/bunch_fitter for details
//
// reads hadded matx root file "rootfiles/all.root"
//
// tbit
//  0 - bbc
//  1 - zdc
//  2 - vpd
// cbit
//  0 - east
//  1 - west
//  2 - coin

void BF(const char * var="mul",
        Int_t numer_tbit=1, Int_t numer_cbit=2, Int_t denom_tbit=2, Int_t denom_cbit=2)
{
  // open hadded matx tree and check its ordering by running DrawMatrix.C
  TFile * infile = new TFile("rootfiles/all.root","READ");
  TTree * matx = (TTree*) infile->Get("matx");
  TFile * countsfile = new TFile("../counts.root","READ");
  TTree * counts = (TTree*) countsfile->Get("sca");
  printf("executing DrawMatrix.C to check ordering of matx tree\n");
  gROOT->ProcessLine(".x DrawMatrix.C");


  // trigger bit character strings (tbit)
  char tbit_s[3][4];
  sprintf(tbit_s[0],"bbc");
  sprintf(tbit_s[1],"zdc");
  sprintf(tbit_s[2],"vpd");


  // combination bit character strings (cbit)
  char cbit_s[3][4];
  sprintf(cbit_s[0],"e");
  sprintf(cbit_s[1],"w");
  sprintf(cbit_s[2],"x");

  // set matx branch addresses
  Int_t i,fi,tbit,cbit,bx;
  Double_t t,val;
  matx->SetBranchAddress("i",&i);
  matx->SetBranchAddress("fi",&fi);
  matx->SetBranchAddress("t",&t);
  matx->SetBranchAddress("tbit",&tbit);
  matx->SetBranchAddress("cbit",&cbit);
  matx->SetBranchAddress("bx",&bx);
  if(!strcmp(var,"raw")) matx->SetBranchAddress("raw",&val);
  else if(!strcmp(var,"acc")) matx->SetBranchAddress("acc",&val);
  else if(!strcmp(var,"mul")) matx->SetBranchAddress("mul",&val);
  else if(!strcmp(var,"fac")) matx->SetBranchAddress("fac",&val);
  else 
  {
    fprintf(stderr,"ERROR: val not valid\n");
    return;
  };


  // define ratio array "rat" and corresponding uncertainties "unc"
  // -- from documentation:  r^i = rat;  sigma_{r^i} = unc
  Int_t imax_tmp = matx->GetMaximum("i");
  const Int_t imax = imax_tmp;
  Double_t rat[imax][120]; // rat = num / den
  Double_t den[imax][120]; // denominator
  Double_t num[imax][120]; // numerator
  Double_t unc[imax][120]; // variances of rat
  for(Int_t ii=0; ii<imax; ii++)
  {
    for(Int_t bb=0; bb<120; bb++)
    {
      num[ii][bb]=0;
      den[ii][bb]=0;
    };
  };
  for(Int_t k=0; k<matx->GetEntries(); k++)
  {
    matx->GetEntry(k);
    if(tbit==numer_tbit && cbit==numer_cbit) num[i-1][bx-1] = val;
    if(tbit==denom_tbit && cbit==denom_cbit) den[i-1][bx-1] = val;
  };
  for(Int_t ii=0; ii<imax; ii++)
  {
    for(Int_t bb=0; bb<120; bb++)
    {
      if(num[ii][bb]>0 && den[ii][bb]>0) 
      {
        rat[ii][bb]=num[ii][bb]/den[ii][bb];
        unc[ii][bb]=rat[ii][bb]*sqrt(1/num[ii][bb]+1/den[ii][bb]);
      }
      else rat[ii][bb]=0;
    };
  };

  
  // rat vs. bXing
  TH1D * rat_v_bx[imax];
  TObjArray * rat_v_bx_arr = new TObjArray();
  char rat_v_bx_n[imax][16];
  char rat_v_bx_t[imax][128];
  for(Int_t ii=0; ii<imax; ii++)
  {
    sprintf(rat_v_bx_n[ii],"rat_v_bx_%d",ii);
    sprintf(rat_v_bx_t[ii],"%s%s/%s%s vs. bXing for i=%d",
      tbit_s[numer_tbit],cbit_s[numer_cbit],tbit_s[denom_tbit],cbit_s[denom_cbit],ii);
    rat_v_bx[ii] = new TH1D(rat_v_bx_n[ii],rat_v_bx_t[ii],120,0,120);
    for(Int_t bb=0; bb<120; bb++)
    {
      rat_v_bx[ii]->SetBinContent(bb-1,rat[ii][bb]);
      rat_v_bx[ii]->SetBinError(bb-1,unc[ii][bb]);
    };
    rat_v_bx_arr->AddLast(rat_v_bx[ii]);
  };


  // set counts branch addresses and fill blue/yell/kicked arrays
  Int_t i_c,bx_c,blue,yell;
  Bool_t kicked;
  Int_t blue_arr[imax][120]; // spin arrays
  Int_t yell_arr[imax][120];
  Bool_t kicked_arr[imax][120];
  counts->SetBranchAddress("i",&i_c);
  counts->SetBranchAddress("bx",&bx_c);
  counts->SetBranchAddress("blue",&blue);
  counts->SetBranchAddress("yell",&yell);
  counts->SetBranchAddress("kicked",&kicked);
  for(Int_t k=0; k<counts->GetEntries(); k++)
  {
    counts->GetEntry(k);
    blue_arr[i_c-1][bx_c]=blue;
    yell_arr[i_c-1][bx_c]=yell;
    kicked[i_c-1][bx_c]=kicked;
  };

  
  // fill spin array H 
  // NOTE: H is automatically set to zero if kicked==true!!!!!
  Int_t H[4][imax][120]; // [asym] [run] [bx]
  for(Int_t ii=0; ii<imax; ii++)
  {
    for(Int_t bb=0; bb<120; bb++)
    {
      if(kicked[ii][bb]==0)
      {
        H[1][ii][bb] = yell_arr[ii][bb];
        H[2][ii][bb] = blue_arr[ii][bb];
        H[3][ii][bb] = yell_arr[ii][bb] * blue_arr[ii][bb];
      }
      else for(Int_t aa=1; aa<4; aa++) H[aa][ii][bb]=0;
    };
  };


  // compute sigma functions for each run, summing over bXings
  Double_t sigma_id[4][imax]; // [asym] [run]
  Double_t sigma_r[4][imax];
  Double_t sigma_h[4][imax];
  Double_t sigma_hr[4][imax];
  for(Int_t aa=1; aa<4; aa++)
  {
    for(Int_t ii=0; ii<imax; ii++)
    {
      sigma_id[aa][ii]=0;
      sigma_r[aa][ii]=0;
      sigma_h[aa][ii]=0;
      sigma_hr[aa][ii]=0;
      for(Int_t bb=0; bb<120; bb++)
      {
        if(H[aa][ii][bb]!=0)
        {
          sigma_id[aa][ii] += 1 / pow(unc[ii][bb],2);
          sigma_r[aa][ii] += rat[ii][bb] / pow(unc[ii][bb],2);
          sigma_h[aa][ii] += H[aa][ii][bb] / pow(unc[ii][bb],2);
          sigma_hr[aa][ii] += H[aa][ii][bb] * rat[ii][bb] / pow(unc[ii][bb],2);
        };
      };
    };
  };


  // sigma function vs. run
  TH1D * sigma_id_dist[4]; // [asym]
  TH1D * sigma_r_dist[4];
  TH1D * sigma_h_dist[4];
  TH1D * sigma_hr_dist[4];
  char sigma_id_dist_n[4][16];
  char sigma_r_dist_n[4][16];
  char sigma_h_dist_n[4][16];
  char sigma_hr_dist_n[4][16];
  char sigma_id_dist_t[4][128];
  char sigma_r_dist_t[4][128];
  char sigma_h_dist_t[4][128];
  char sigma_hr_dist_t[4][128];
  for(Int_t aa=1; aa<4; aa++)
  {
    sprintf(sigma_id_dist_n[aa],"sigma_id_a%d",aa);
    sprintf(sigma_r_dist_n[aa],"sigma_r_a%d",aa);
    sprintf(sigma_h_dist_n[aa],"sigma_h_a%d",aa);
    sprintf(sigma_hr_dist_n[aa],"sigma_hr_a%d",aa);
    sprintf(sigma_id_dist_t[aa],"#Sigma(1) vs. run // asym=%d",aa);
    sprintf(sigma_r_dist_t[aa],"#Sigma(r) vs. run // asym=%d",aa);
    sprintf(sigma_h_dist_t[aa],"#Sigma(H) vs. run // asym==%d",aa);
    sprintf(sigma_hr_dist_t[aa],"#Sigma(Hr) vs. run // asym==%d",aa);
    sigma_id_dist[aa] = new TH1D(sigma_id_dist_n[aa],sigma_id_dist_t[aa],imax,1,imax+1);
    sigma_r_dist[aa] = new TH1D(sigma_r_dist_n[aa],sigma_r_dist_t[aa],imax,1,imax+1);
    sigma_h_dist[aa] = new TH1D(sigma_h_dist_n[aa],sigma_h_dist_t[aa],imax,1,imax+1);
    sigma_hr_dist[aa] = new TH1D(sigma_hr_dist_n[aa],sigma_hr_dist_t[aa],imax,1,imax+1);
    for(Int_t ii=0; ii<imax; ii++)
    {
      sigma_id_dist[aa]->SetBinContent(ii+1,sigma_id[aa][ii]);
      sigma_r_dist[aa]->SetBinContent(ii+1,sigma_r[aa][ii]);
      sigma_h_dist[aa]->SetBinContent(ii+1,sigma_h[aa][ii]);
      sigma_hr_dist[aa]->SetBinContent(ii+1,sigma_hr[aa][ii]);
    };
  };









  // write output
  TFile * outfile = new TFile("fit_result.root","RECREATE");
  rat_v_bx_arr->Write("rat_v_bx_arr",TObject::kSingleKey);
  for(Int_t aa=1; aa<4; aa++) 
  {
    sigma_id_dist[aa]->Write();
    sigma_r_dist[aa]->Write();
    sigma_h_dist[aa]->Write();
    sigma_hr_dist[aa]->Write();
  };

  



};
