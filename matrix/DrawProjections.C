// draws distribution of asymmetry determined from bunch fit
// 
// GENERALLY YOU WANT TO USE THIS VERSION
//
// requires the following files, all of which are produced from RunPatterns:
//  - fit_result.[numer].[denom].root
//  - pats/fit_result.[numer].[denom].pat*.root
//  - these files are then to be stored in the subdirectory 
//    output/omit_[Nomit]_bXings_after_aborts/
//
//  - if Nomit = -1 (default value), then it just reads the files from ./

void DrawProjections(const char * numer="vpdx",
                     const char * denom="zdcx",
                     Int_t Nomit=-1)
{
  const Float_t WIDTH = 4e-3;
  gStyle->SetOptFit(1);
  char filename[256];
  if(Nomit==-1)
    sprintf(filename,"fit_result.%s.%s.root",numer,denom);
  else 
    sprintf(filename,"output/omit_%d_bXings_after_aborts/fit_result.%s.%s.root",
      Nomit,numer,denom);
  TFile * infile = new TFile(filename,"READ");
  TH1D * hist_all = (TH1D*) infile->Get("/asymmetry/asym_a3_v_run");

  TFile * patfile[8];
  char patfile_n[8][256];
  Int_t pat_arr[8] = {13,14,23,24,31,32,41,42};
  Color_t colours[8] = {kOrange,kRed,kMagenta,kBlue,kCyan+1,kGreen+1,kYellow+2,kViolet-6}; 
  TH1D * hist_pat[8];
  for(Int_t pp=0; pp<8; pp++)
  {
    if(Nomit==-1)
      sprintf(patfile_n[pp],"pats/fit_result.%s.%s.pat%d.root",numer,denom,pat_arr[pp]);
    else
      sprintf(patfile_n[pp],"output/omit_%d_bXings_after_aborts/pats/fit_result.%s.%s.pat%d.root",
        Nomit,numer,denom,pat_arr[pp]);
    patfile[pp] = new TFile(patfile_n[pp],"READ");
    hist_pat[pp] = (TH1D*) patfile[pp]->Get("/asymmetry/asym_a3_v_run");
  };

  char omit_str[64];
  if(Nomit==-1) strcpy(omit_str,"");
  else strcpy(omit_str," :: omitting %d bx after aborts",Nomit);

  char asym_dist_all_n[256];
  char asym_dist_all_t[256];
  sprintf(asym_dist_all_n,"asym_dist_all");
  //sprintf(asym_dist_all_t,"Run 13 S_{LL} distribution for %s/%s%s",numer,denom,omit_str);
  sprintf(asym_dist_all_t,"Run 13 S_{LL} distribution for rate-safe %s/%s%s",numer,denom,omit_str);
  TH1D * asym_dist_all = new TH1D(asym_dist_all_n,asym_dist_all_t,100,-1*WIDTH,WIDTH);
  asym_dist_all->SetLineWidth(3);

  char asym_dist_pat_n[8][256];
  char asym_dist_pat_t[8][256];
  TH1D * asym_dist_pat[8];
  for(Int_t pp=0; pp<8; pp++)
  {
    sprintf(asym_dist_pat_n[pp],"asym_dist_pat%d",pp);
    //sprintf(asym_dist_pat_t[pp],"Run 13 S_{LL} distribution for %s/%s%s",numer,denom,omit_str);
    sprintf(asym_dist_pat_t[pp],"Run 13 S_{LL} distribution for rate-safe %s/%s%s",numer,denom,omit_str);
    asym_dist_pat[pp] = new TH1D(asym_dist_pat_n[pp],asym_dist_pat_t[pp],100,-1*WIDTH,WIDTH);
    asym_dist_pat[pp]->SetLineColor(colours[pp]);
  };

  Double_t bc;

  for(Int_t b=1; b<=hist_all->GetNbinsX(); b++)
  {
    bc = hist_all->GetBinContent(b);
    asym_dist_all->Fill(bc);
  };

  
  for(Int_t pp=0; pp<8; pp++)
  {
    for(Int_t b=1; b<=hist_pat[pp]->GetNbinsX(); b++)
    {
      bc = hist_pat[pp]->GetBinContent(b);
      if(bc!=0) asym_dist_pat[pp]->Fill(bc);
    };
  };

  // run 13 fits (2 gaussians)
  //TF1 * gaus1 = new TF1("gaus1","gaus",-1*WIDTH,0);
  //TF1 * gaus2 = new TF1("gaus2","gaus",0,WIDTH);
  TF1 * twogaus = 
    new TF1("twogaus","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)");
  twogaus->SetParameter(0,80); // normalisation
  twogaus->SetParameter(3,80); 
  twogaus->SetParameter(1,0.002); // mean
  twogaus->SetParameter(4,-0.002); 
  twogaus->SetParameter(2,0.0005); // sigma
  twogaus->SetParameter(5,0.0005);

  twogaus->SetParNames("N_{L}","#mu_{L}","#sigma_{L}","N_{R}","#mu_{R}","#sigma_{R}");

  //asym_dist_all->Fit(gaus1,"","",-1*WIDTH,0);
  //asym_dist_all->Fit(gaus2,"","",0,WIDTH);
  asym_dist_all->Fit(twogaus,"","",-1*WIDTH,WIDTH);

  c1->Close();


  TCanvas * cc = new TCanvas("cc","cc",1000,1000);
  cc->SetGrid(1,0);
  gStyle->SetOptStat(1100);
  //gStyle->SetStatFontSize(0.1);
  Float_t size = 0.05;
  asym_dist_all->GetXaxis()->SetLabelSize(size);
  asym_dist_all->GetYaxis()->SetLabelSize(size);
  asym_dist_all->Draw();
  //gaus1->Draw("same");
  //gaus2->Draw("same");
  twogaus->Draw("same");
  for(Int_t pp=0; pp<8; pp++)
  {
    asym_dist_pat[pp]->GetXaxis()->SetLabelSize(size);
    asym_dist_pat[pp]->GetYaxis()->SetLabelSize(size);
    asym_dist_pat[pp]->Draw("same");
  };
  asym_dist_all->Draw("same");

  /*
  Double_t sigma1 = gaus1->GetParameter(2);
  Double_t sigma2 = gaus2->GetParameter(2);
  Double_t sigma = (sigma1>sigma2)?sigma1:sigma2;
  Double_t mean = asym_dist_all->GetMean();
  Double_t sys = sigma + fabs(mean);
  printf("sigma1 = %f\nsigma2 = %f\n",sigma1,sigma2);
  printf("sigma = %f\nmean = %f\nsystematic = %f\n",sigma,mean,sys);
  */ 


  char printname[64];
  if(Nomit==-1) strcpy(printname,"projection.png");
  else sprintf(printname,"projection.omit%d.png",Nomit);
  cc->Print(printname,"png");
    
};
