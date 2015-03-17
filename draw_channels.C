void draw_channels(const char * filename="chtr.root")
{
  gStyle->SetOptStat(0);
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("tr");
  
  TFile * outfile = new TFile("test_tree.root","RECREATE");
  TTree * outtr = new TTree("outtr","outtr");
  Double_t br_asym,br_integral;
  Int_t br_iterator,br_chBBC,br_chZDC,br_chVPD;
  outtr->Branch("asym",&br_asym,"asym/D");
  outtr->Branch("integral",&br_integral,"integral/D");
  outtr->Branch("iterator",&br_iterator,"iterator/I");
  outtr->Branch("chBBC",&br_chBBC,"chBBC/I");
  outtr->Branch("chZDC",&br_chZDC,"chZDC/I");
  outtr->Branch("chVPD",&br_chVPD,"chVPD/I");

  const Int_t nchBBC = 5; //5
  const Int_t nchZDC = 8; //8
  const Int_t nchVPD = 4; //4

  Int_t chBBCseq[nchBBC] = {0,1,2,3,7};
  Int_t chZDCseq[nchZDC] = {0,1,2,3,4,5,6,7};
  Int_t chVPDseq[nchVPD] = {0,1,2,3};

  //Int_t chBBCseq[nchBBC] = {0};
  //Int_t chZDCseq[nchZDC] = {0};
  //Int_t chVPDseq[nchVPD] = {0};

  TH1D * nbxbx[nchBBC][nchZDC][nchVPD];
  TH1D * total = new TH1D("total","total",120,0,120);
  TH1D * totalcl = new TH1D("totalcl","totalcl",120,0,120);
  TGraph * asym_total = new TGraph();
  TGraph * asym_nbxbx = new TGraph();
  char asym_str[256];
  char asym_total_t[256];
  char asym_nbxbx_t[256];
  strcpy(asym_str,"[ #int_{0}^{60}N_{bx}d(bx) - #int_{61}^{119}N_{bx}d(bx) ] / #int_{1}^{119}N_{bx}d(bx)");
  sprintf(asym_total_t,"%s for #Sigma",asym_str);
  sprintf(asym_nbxbx_t,"%s for summand",asym_str);
  asym_total->SetMarkerColor(kRed);
  asym_total->SetMarkerStyle(kFullCircle);
  asym_total->SetTitle(asym_total_t);
  asym_nbxbx->SetMarkerColor(kBlue);
  asym_nbxbx->SetMarkerStyle(kFullCircle);
  asym_nbxbx->SetTitle(asym_nbxbx_t);
  Int_t asym_total_i = 0;
  Int_t asym_nbxbx_i = 0;
  char nbxbx_n[nchBBC][nchZDC][nchVPD][32];
  char nbxbx_t[nchBBC][nchZDC][nchVPD][128];
  char chcut[nchBBC][nchZDC][nchVPD][128];
  Double_t maxval[nchBBC][nchZDC][nchVPD];
  for(Int_t bbc=0; bbc<nchBBC; bbc++)
  {
    for(Int_t zdc=0; zdc<nchZDC; zdc++)
    {
      for(Int_t vpd=0; vpd<nchVPD; vpd++)
      {
        maxval[bbc][zdc][vpd]=0;

        if(nchBBC>1 && nchZDC>1 && nchVPD>1)
        {
        sprintf(nbxbx_n[bbc][zdc][vpd],"nbxbx_b%d_z%d_v%d",
          chBBCseq[bbc],chZDCseq[zdc],chVPDseq[vpd]);
        sprintf(nbxbx_t[bbc][zdc][vpd],"N_{bx} vs. bXing no. -- chBBC==%d chZDC==%d chVPD==%d",
          chBBCseq[bbc],chZDCseq[zdc],chVPDseq[vpd]);
          sprintf(chcut[bbc][zdc][vpd],"count*(chBBC==%d&&chZDC==%d&&chVPD==%d)",
            chBBCseq[bbc],chZDCseq[zdc],chVPDseq[vpd]);
        }
        else if(nchBBC==1 && nchZDC>1 && nchVPD>1)
        {
        sprintf(nbxbx_n[bbc][zdc][vpd],"nbxbx_b%d_z%d_v%d",
          chBBCseq[bbc],chZDCseq[zdc],chVPDseq[vpd]);
        sprintf(nbxbx_t[bbc][zdc][vpd],"N_{bx} vs. bXing no. -- chBBC==any chZDC==%d chVPD==%d",
          chZDCseq[zdc],chVPDseq[vpd]);
          sprintf(chcut[bbc][zdc][vpd],"count*(chZDC==%d&&chVPD==%d)",
            chZDCseq[zdc],chVPDseq[vpd]);
        }
        else if(nchBBC>1 && nchZDC==1 && nchVPD>1)
        {
        sprintf(nbxbx_n[bbc][zdc][vpd],"nbxbx_b%d_z%d_v%d",
          chBBCseq[bbc],chZDCseq[zdc],chVPDseq[vpd]);
        sprintf(nbxbx_t[bbc][zdc][vpd],"N_{bx} vs. bXing no. -- chBBC==%d chZDC==any chVPD==%d",
          chBBCseq[bbc],chVPDseq[vpd]);
          sprintf(chcut[bbc][zdc][vpd],"count*(chBBC==%d&&chVPD==%d)",
            chBBCseq[bbc],chVPDseq[vpd]);
        }
        else if(nchBBC>1 && nchZDC>1 && nchVPD==1)
        {
        sprintf(nbxbx_n[bbc][zdc][vpd],"nbxbx_b%d_z%d_v%d",
          chBBCseq[bbc],chZDCseq[zdc],chVPDseq[vpd]);
        sprintf(nbxbx_t[bbc][zdc][vpd],"N_{bx} vs. bXing no. -- chBBC==%d chZDC==%d chVPD==any",
          chBBCseq[bbc],chZDCseq[zdc]);
          sprintf(chcut[bbc][zdc][vpd],"count*(chBBC==%d&&chZDC==%d)",
            chBBCseq[bbc],chZDCseq[zdc]);
        };

        nbxbx[bbc][zdc][vpd] = new TH1D(nbxbx_n[bbc][zdc][vpd],nbxbx_t[bbc][zdc][vpd],120,0,120);
        tr->Project(nbxbx_n[bbc][zdc][vpd],"bx",chcut[bbc][zdc][vpd]);
        maxval[bbc][zdc][vpd] = nbxbx[bbc][zdc][vpd]->Integral();
        printf("%s integral=%f\n",chcut[bbc][zdc][vpd],maxval[bbc][zdc][vpd]);
      };
    };
  };

  TCanvas * nbxbx_canv = new TCanvas("nbxbx_canv","nbxbx_canv",500,500);
  TCanvas * bits_canv = new TCanvas("bits_canv","bits_canv",850,520);
  TCanvas * asym_canv = new TCanvas("asym_canv","asym_canv",650,700);
  bits_canv->Divide(4,2);
  asym_canv->Divide(1,2);
  asym_canv->GetPad(1)->SetGrid(1,1);
  asym_canv->GetPad(2)->SetGrid(1,1);

  TH1D * bbcChDist = new TH1D("bbcChDist","bbc scaler bits",8,0,8);
  TH1D * zdcChDist = new TH1D("zdcChDist","zdc scaler bits",8,0,8);
  TH1D * vpdChDist = new TH1D("vpdChDist","vpd scaler bits",8,0,8);
  TH1D * bbcLogDist = new TH1D("bbcLogDist","bbc logic (0=zero 1=logical 2=illogical)",3,0,3);
  TH1D * zdcLogDist = new TH1D("zdcLogDist","zdc logic (0=zero 1=logical 2=illogical)",3,0,3);
  TH1D * vpdLogDist = new TH1D("vpdLogDist","vpd logic (0=zero 1=logical 2=illogical)",3,0,3);
  TH1D * logic_dist = new TH1D("logic_dist","isIllogical bit boolean",2,0,2);


  Int_t iterator=0;
  Int_t iterator_max = nchBBC*nchZDC*nchVPD;
  Double_t curr_max=0;
  Int_t bbc_set,zdc_set,vpd_set;
  Double_t minb,maxb;
  char iterator_title[32];
  Double_t numer,denom;
  while(iterator < iterator_max)
  {
    for(Int_t bbc=0; bbc<nchBBC; bbc++)
    {
      for(Int_t zdc=0; zdc<nchZDC; zdc++)
      {
        for(Int_t vpd=0; vpd<nchVPD; vpd++)
        {
          if(maxval[bbc][zdc][vpd] > curr_max)
          {
            curr_max = maxval[bbc][zdc][vpd];
            bbc_set = bbc;
            zdc_set = zdc;
            vpd_set = vpd;
          };
        };
      };
    };
    maxval[bbc_set][zdc_set][vpd_set] = 0;
    curr_max = 0;

    if(iterator>=0)
    {
      bbcChDist->Fill(chBBCseq[bbc_set]);
      zdcChDist->Fill(chZDCseq[zdc_set]);
      vpdChDist->Fill(chVPDseq[vpd_set]);
      bbcLogDist->Fill(Logical(0,chBBCseq[bbc_set]));
      zdcLogDist->Fill(Logical(1,chZDCseq[zdc_set]));
      vpdLogDist->Fill(Logical(2,chVPDseq[vpd_set]));
      if(Logical(0,chBBCseq[bbc_set])==2 ||
         Logical(1,chZDCseq[zdc_set])==2 ||
         Logical(2,chVPDseq[vpd_set])==2) logic_dist->Fill(1);
      else logic_dist->Fill(0);

      total->Add(nbxbx[bbc_set][zdc_set][vpd_set]);
      totalcl = (TH1D*) total->Clone();
      totalcl->Scale(1.0/totalcl->GetMaximum());
      totalcl->SetTitle(nbxbx_t[bbc_set][zdc_set][vpd_set]);

      numer = nbxbx[bbc_set][zdc_set][vpd_set]->Integral(1,31) +
              nbxbx[bbc_set][zdc_set][vpd_set]->Integral(41,61);
      denom = nbxbx[bbc_set][zdc_set][vpd_set]->Integral(61,112);
      br_asym = (numer-denom)/(numer+denom);
      asym_nbxbx->SetPoint(asym_nbxbx_i,iterator,br_asym);
      asym_nbxbx_i++;

      numer = totalcl->Integral(1,31) +
              totalcl->Integral(41,61);
      denom = totalcl->Integral(61,112);
      asym_total->SetPoint(asym_total_i,iterator,(numer-denom)/(numer+denom));
      asym_total_i++;
      //asym_total->GetHistogram()->SetMinimum(0.996);
      //asym_total->GetHistogram()->SetMaximum(1.001);

      sprintf(iterator_title,"iterator=%d",iterator);
      totalcl->GetXaxis()->SetTitle(iterator_title);

      br_iterator = iterator;
      br_integral = nbxbx[bbc_set][zdc_set][vpd_set]->Integral(1,120);
      br_chBBC = chBBCseq[bbc_set];
      br_chZDC = chZDCseq[zdc_set];
      br_chVPD = chVPDseq[vpd_set];
      outtr->Fill();


      minb = 1;
      maxb = 0;
      for(Int_t b=1; b<=120; b++)
      {
        if(!(b>=27&&b<=42) && !(b>100))
        {
          minb = (totalcl->GetBinContent(b)<minb) ? totalcl->GetBinContent(b):minb;
          maxb = (totalcl->GetBinContent(b)>maxb) ? totalcl->GetBinContent(b):maxb;
        };
      };
      totalcl->GetYaxis()->SetRangeUser(minb-0.00001,maxb);
      bbcChDist->GetYaxis()->SetRangeUser(0,35);
      zdcChDist->GetYaxis()->SetRangeUser(0,23);
      vpdChDist->GetYaxis()->SetRangeUser(0,45);
      bbcLogDist->GetYaxis()->SetRangeUser(0,98);
      zdcLogDist->GetYaxis()->SetRangeUser(0,82);
      vpdLogDist->GetYaxis()->SetRangeUser(0,121);
      logic_dist->GetYaxis()->SetRangeUser(0,97);

      nbxbx_canv->cd(); totalcl->Draw();
      bits_canv->cd(1); bbcChDist->Draw();
      bits_canv->cd(2); zdcChDist->Draw();
      bits_canv->cd(3); vpdChDist->Draw();
      bits_canv->cd(5); bbcLogDist->Draw();
      bits_canv->cd(6); zdcLogDist->Draw();
      bits_canv->cd(7); vpdLogDist->Draw();
      bits_canv->cd(4); logic_dist->Draw();
      asym_canv->cd(1); asym_nbxbx->Draw("AP");
      asym_canv->cd(2); asym_total->Draw("AP");
      nbxbx_canv->Update();
      bits_canv->Update();
      asym_canv->Update();
      //gSystem->Sleep(50);

      printf("%d\t%d %d %d\n",iterator,chBBCseq[bbc_set],chZDCseq[zdc_set],chVPDseq[vpd_set]);
    };
    iterator++;
  };

  //asym_total->Fit("pol0","","",0,iterator);
  asym_nbxbx->Fit("pol0","","",0,110);
  asym_canv->Update();

  //return;

  TCanvas * print_canv = new TCanvas("print_canv","print_canv",500,500);
  char printname[64];
  for(Int_t bbc=0; bbc<nchBBC; bbc++)
  {
    for(Int_t zdc=0; zdc<nchZDC; zdc++)
    {
      for(Int_t vpd=0; vpd<nchVPD; vpd++)
      {
        if(bbc==0 && zdc==0 && vpd==0) strcpy(printname,"out.pdf(");
        else strcpy(printname,"out.pdf");
        nbxbx[bbc][zdc][vpd]->Draw();
        print_canv->Print(printname,"pdf");
      };
    };
  };
  print_canv->Clear();
  print_canv->Print("out.pdf)","pdf");
  print_canv->Close();
  outtr->Write();


  char asym_vs_chBBC_t[256];
  char asym_vs_chZDC_t[256];
  char asym_vs_chVPD_t[256];
  sprintf(asym_vs_chBBC_t,"%s vs. %s",asym_str,"chBBC");
  sprintf(asym_vs_chZDC_t,"%s vs. %s",asym_str,"chZDC");
  sprintf(asym_vs_chVPD_t,"%s vs. %s",asym_str,"chVPD");
  TH2D * asym_vs_chBBC = new TH2D("asym_vs_chBBC",asym_vs_chBBC_t,8,0,8,100,-0.5,0.5);
  TH2D * asym_vs_chZDC = new TH2D("asym_vs_chZDC",asym_vs_chZDC_t,8,0,8,100,-0.5,0.5);
  TH2D * asym_vs_chVPD = new TH2D("asym_vs_chVPD",asym_vs_chVPD_t,8,0,8,100,-0.5,0.5);
  outtr->Project("asym_vs_chBBC","asym:chBBC","");
  outtr->Project("asym_vs_chZDC","asym:chZDC","");
  outtr->Project("asym_vs_chVPD","asym:chVPD","");
  TProfile * asym_vs_chBBC_pfx = asym_vs_chBBC->ProfileX();
  TProfile * asym_vs_chZDC_pfx = asym_vs_chZDC->ProfileX();
  TProfile * asym_vs_chVPD_pfx = asym_vs_chVPD->ProfileX();
  asym_vs_chBBC_pfx->SetLineWidth(2);
  asym_vs_chZDC_pfx->SetLineWidth(2);
  asym_vs_chVPD_pfx->SetLineWidth(2);
  TCanvas * asym_vs_chBBC_canv = new TCanvas("asym_vs_chBBC_canv","asym_vs_chBBC_canv",600,600);
    asym_vs_chBBC_canv->SetGrid(1,1);
    asym_vs_chBBC->Draw("colz");
    asym_vs_chBBC_pfx->Draw("same");
  TCanvas * asym_vs_chZDC_canv = new TCanvas("asym_vs_chZDC_canv","asym_vs_chZDC_canv",600,600);
    asym_vs_chZDC_canv->SetGrid(1,1);
    asym_vs_chZDC->Draw("colz");
    asym_vs_chZDC_pfx->Draw("same");
  TCanvas * asym_vs_chVPD_canv = new TCanvas("asym_vs_chVPD_canv","asym_vs_chVPD_canv",600,600);
    asym_vs_chVPD_canv->SetGrid(1,1);
    asym_vs_chVPD->Draw("colz");
    asym_vs_chVPD_pfx->Draw("same");
  asym_vs_chBBC_canv->Write();
  asym_vs_chZDC_canv->Write();
  asym_vs_chVPD_canv->Write();
  asym_vs_chBBC_canv->Close();
  asym_vs_chZDC_canv->Close();
  asym_vs_chVPD_canv->Close();
  nbxbx_canv->Write();
  bits_canv->Write();
  asym_canv->Write();
}

Int_t Logical(Int_t det, Int_t sbc)
{
  if(sbc==0) return 0;
  else if(det<2)
  {
    if(sbc==1 || sbc==2 || sbc==7) return 1;
    else return 2;
  }
  else return 1;
};
