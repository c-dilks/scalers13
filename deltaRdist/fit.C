void fit(Int_t num_post_aborts=0) {
  TString fn = Form("data_%d.root",num_post_aborts);
  TFile * f = new TFile(fn.Data(),"READ");
  TH1D * d = (TH1D*) f->Get("deltadist_rsc_zdc_vpd_3");

  Float_t paramset[6];

  TString formu;
  Int_t ngaus = 0;

  if(num_post_aborts <= 1) {
    ngaus = 2;
    paramset[0] = 80; // normalisation
    paramset[1] = 0.002; // mean
    paramset[2] = 0.0005; // sigma
    paramset[3] = 80; // norm
    paramset[4] = -0.002; // mean
    paramset[5] = 0.0005; // sigma
  } else if(num_post_aborts == 2 ||
            num_post_aborts == 1000) {
    ngaus = 2;
    paramset[0] = 80; // normalisation
    paramset[1] = 0.001; // mean
    paramset[2] = 0.0002; // sigma
    paramset[3] = 80; // norm
    paramset[4] = -0.001; // mean
    paramset[5] = 0.0002; // sigma
  } else if(num_post_aborts == 3 ||
            num_post_aborts == 4 ||
            num_post_aborts == 5 ||
            num_post_aborts == 6 ||
            num_post_aborts == 10 ||
            num_post_aborts == 11 ||
            num_post_aborts == 12 ||
            num_post_aborts == 13 ||
            num_post_aborts == 14 ||
            num_post_aborts == 18 ||
            num_post_aborts == 19 ||
            num_post_aborts == 20 ||
            num_post_aborts == 1000) {
    ngaus = 1;
    paramset[0] = 80; // normalisation
    paramset[1] = 0; // mean
    paramset[2] = 0.0002; // sigma
    paramset[3] = 0; // ignored
    paramset[4] = 0; // ignored
    paramset[5] = 0; // ignored
  } else if(num_post_aborts == 7 ||
            num_post_aborts == 8 ||
            num_post_aborts == 9 ||
            num_post_aborts == 12 ||
            num_post_aborts == 15 ||
            num_post_aborts == 16 ||
            num_post_aborts == 17 ||
            num_post_aborts == 1000) {
    ngaus = 2;
    paramset[0] = 80; // normalisation
    paramset[1] = 0.0005; // mean
    paramset[2] = 0.0002; // sigma
    paramset[3] = 80; // norm
    paramset[4] = -0.0005; // mean
    paramset[5] = 0.0002; // sigma
  } 

  Bool_t bad = false;
  if(num_post_aborts==4 || 
     num_post_aborts==10 ||
     num_post_aborts==12 ||
     num_post_aborts==18 ||
     num_post_aborts==20 ||
     num_post_aborts==10000) {
    printf("\n\nFIT FAILING\n\n");
    bad = true;
  };

  TF1 * ftn; 

  switch(ngaus) {
    case 1:
      formu = "[0]*exp(-0.5*((x-[1])/[2])^2)";
      ftn = new TF1("fit",formu.Data());
      ftn->SetParNames("N","#mu","#sigma");
      break;
    case 2:
      formu = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)"; 
      ftn = new TF1("fit",formu.Data());
      ftn->SetParNames("N_{L}","#mu_{L}","#sigma_{L}","N_{R}","#mu_{R}","#sigma_{R}");
      break;
    default:
      fprintf(stderr,"error: ngaus=0\n");
      formu = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)"; 
      ftn = new TF1("fit",formu.Data());
      ftn->SetParNames("N_{L}","#mu_{L}","#sigma_{L}","N_{R}","#mu_{R}","#sigma_{R}");
      break;
  };

  /*
  if(num_post_aborts==4) {
    ftn->SetParLimits(2,0,1e-4)
    ftn->SetParLimits(5,0,1e-4)
  };
  */

  for(int p=0; p<(3*ngaus); p++) ftn->SetParameter(p,paramset[p]);


  Float_t xmax = d->GetXaxis()->GetXmax();
  Float_t xmin = d->GetXaxis()->GetXmin();

  d->Fit(ftn,"","",xmin,xmax);
  if(c1) c1->Close();

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  TCanvas * c = new TCanvas("c","c",1000,500);
  c->SetGrid(1,1);
  d->Draw();
  TString pngname = num_post_aborts<10 ?
    Form("output_0%d.png",num_post_aborts) : Form("output_%d.png",num_post_aborts);
  c->Print(pngname.Data(),"png");

  gSystem->RedirectOutput("fitdata.dat","a");
  switch(ngaus) {
    case 1:
      printf("%d %d %f %f\n",
        num_post_aborts,
        bad?-1:0,
        ftn->GetParameter(1),
        ftn->GetParameter(2)
      );
      break;
    case 2:
      printf("%d %d %f %f\n",
        num_post_aborts,
        bad?-1:1,
        ftn->GetParameter(1),
        ftn->GetParameter(2)
      );
      printf("%d %d %f %f\n",
        num_post_aborts,
        bad?-1:2,
        ftn->GetParameter(4),
        ftn->GetParameter(5)
      );
      break;
  };
  gSystem->RedirectOutput(0);
};
