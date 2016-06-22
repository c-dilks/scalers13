// investigates the "pattern of threes" seen in the ratio of scaler counts
// vs. bX

void Threes(const char * filename="fit_result.zdcx.vpdx.root")
{
  // open ratio vs. bx plots
  TFile * infile = new TFile(filename,"READ");
  TObjArray * oa = (TObjArray*) infile->Get("rat_v_bx_arr");
  Int_t ENT_tmp = oa->GetEntries();
  const Int_t ENT = ENT_tmp;
  TH1D * rat[ENT];
  Int_t i,j,k,phase,bx,bn;
  Double_t bc;
  for(i=0; i<ENT; i++) rat[i] = (TH1D*)(oa->At(i));

  // pattern no. --> integer between 0-7
  Int_t key[100];
  Int_t rev[8];
  for(j=0; j<100; j++) key[j]=0;
  key[13]=0; rev[0]=13;  // 1 
  key[14]=1; rev[1]=14;  // 2
  key[23]=2; rev[2]=23;  // 3
  key[24]=3; rev[3]=24;  // 4
  key[31]=4; rev[4]=31;  //  5
  key[32]=5; rev[5]=32;  //  7
  key[41]=6; rev[6]=41;  //  6
  key[42]=7; rev[7]=42;  //  8

  
  // obtain spin patterns
  TFile * sumsfile = new TFile("../sums.root","READ");
  TTree * sums = (TTree*) sumsfile->Get("sum");
  Int_t pattern_tr,index;
  Int_t pattern[ENT];
  sums->SetBranchAddress("i",&index);
  sums->SetBranchAddress("pattern",&pattern_tr);
  for(i=0; i<sums->GetEntries(); i++)
  {
    sums->GetEntry(i);
    if(index>0 && index<=ENT) pattern[index-1] = pattern_tr;
    else 
    {
      fprintf(stderr,"sums tree index doesn't agree with obj array\n");
      return;
    };
  };


  /* bin contents - average over 6 bXs*/
  Double_t con[20][6][6][ENT]; // [(bx+phase)/6] [(bx+phase)%6] [phase] [run]

  /* average over 6 bXs */
  Double_t ave[20][6][ENT]; // [(bx+phase)/6] [phase] [run]


  // initialise "ave"
  for(i=0; i<ENT; i++)
  {
    for(j=0; j<20; j++)
    {
      for(phase=0; phase<6; phase++)
      {
        ave[j][phase][i] = 0;
      };
    };
  };
     

  // fill "con" and "ave" arrays
  for(i=0; i<ENT; i++)
  {
    for(phase=0; phase<6; phase++)
    {
      for(bn=1; bn<=120; bn++)
      {
        bx=bn-1;
        bc = rat[i]->GetBinContent(bn);
        con[((bx+phase)%120)/6][((bx+phase)%120)%6][phase][i] = bc;
        ave[((bx+phase)%120)/6][phase][i] += bc;
      };
      for(j=0; j<20; j++) 
      {
        ave[j][phase][i] /= 6;
        for(k=1; k<6; k++) con[j][k][phase][i] /= ave[j][phase][i];
      };
    };
  };


  // initialise analysis distributions
  Int_t NBINS = 100;
  //Double_t LBOUND = -1.5e-3;
  //Double_t UBOUND = 1.5e-3;
  Double_t LBOUND = 1-6e-3;
  Double_t UBOUND = 1+6e-3;
  TH1D * dist[8][6]; // [spin pattern] [phase]
  TString dist_n[8][6];
  TString dist_t[8][6];
  TH1D * distL[8][6];
  TString distL_n[8][6];
  TH1D * distR[8][6];
  TString distR_n[8][6];
  TH2D * dist_vs_run[8][6];
  TString dist_vs_run_n[8][6];
  Int_t s;
  for(s=0; s<8; s++)
  {
    for(phase=0; phase<6; phase++)
    {
      dist_n[s][phase] = Form("dist_pat%d_phase%d",rev[s],phase);
      distL_n[s][phase] = Form("distL_pat%d_phase%d",rev[s],phase);
      distR_n[s][phase] = Form("distR_pat%d_phase%d",rev[s],phase);
      dist_vs_run_n[s][phase] = Form("dist_vs_run_pat%d_phase%d",rev[s],phase);
      dist_t[s][phase] = 
        Form("ratio / <ratio over 6 bXs> :: pattern=%d :: phase=%d",rev[s],phase);
      dist[s][phase] = new TH1D(dist_n[s][phase].Data(),dist_t[s][phase].Data(),
        NBINS,LBOUND,UBOUND);
      distL[s][phase] = new TH1D(distL_n[s][phase].Data(),dist_t[s][phase].Data(),
        NBINS,LBOUND,UBOUND);
      distR[s][phase] = new TH1D(distR_n[s][phase].Data(),dist_t[s][phase].Data(),
        NBINS,LBOUND,UBOUND);
      dist_vs_run[s][phase] = new TH2D(dist_vs_run_n[s][phase].Data(),dist_vs_run_n[s][phase].Data(),
        ENT/8+1,0,ENT/8+1,NBINS/2,LBOUND,UBOUND);

      distL[s][phase]->SetLineColor(kRed);
      distR[s][phase]->SetLineColor(kBlue);
      dist[s][phase]->SetLineColor(kBlack);
      dist[s][phase]->SetLineWidth(2);
    };
  };


  // fill analysis distributions
  Int_t i_for_pat[8];
  Int_t bbx;
  for(Int_t qq=0; qq<8; qq++) i_for_pat[qq]=0;
  for(i=0; i<ENT; i++)
  {
    for(j=0; j<20; j++)
    {
      for(k=0; k<6; k++)
      {
        for(phase=0; phase<6; phase++)
        {
          if(pattern[i]!=0)
          {
            bbx = 6*j+k;
            if((bbx>=6 && bbx<=24) || (bbx>=46 && bbx<=104))
            {
              if(k<3) distL[key[pattern[i]]][phase]->Fill(con[j][k][phase][i]);
              else    distR[key[pattern[i]]][phase]->Fill(con[j][k][phase][i]);
              dist[key[pattern[i]]][phase]->Fill(con[j][k][phase][i]);
              dist_vs_run[key[pattern[i]]][phase]->Fill(i_for_pat[key[pattern[i]]],con[j][k][phase][i]);
            };
          };
        };
      };
    };
    i_for_pat[key[pattern[i]]]++;
  };


  // write analysis distributions
  TFile * outfile = new TFile("three.root","RECREATE");
  TObjArray * dist_arr[6]; // [phase];
  TString dist_arr_n[6];
  TCanvas * canv_dist[6]; // [phase]
  TString canv_dist_n[6];
  TCanvas * canv_dist_vs_run[6]; // [phase]
  TString canv_dist_vs_run_n[6];
  for(phase=0; phase<6; phase++)
  {
    dist_arr_n[phase] = Form("dist_arr_phase%d",phase);
    dist_arr[phase] = new TObjArray();
    for(s=0; s<8; s++) dist_arr[phase]->Add(dist[s][phase]);
    dist_arr[phase]->Write(dist_arr_n[phase].Data(),TObject::kSingleKey);
  };
  Int_t pads[8] = {1,2,3,4,5,7,6,8};
  for(phase=0; phase<6; phase++)
  {
    canv_dist_n[phase] = Form("canv_dist_phase%d",phase);
    canv_dist[phase] = new TCanvas(canv_dist_n[phase].Data(),canv_dist_n[phase].Data(),1200,400);
    canv_dist[phase]->Divide(4,2);
    for(s=0; s<8; s++)
    {
      canv_dist[phase]->cd(pads[s]);
      dist[s][phase]->Draw();
      distL[s][phase]->Draw("same");
      distR[s][phase]->Draw("same");
    };
    canv_dist[phase]->Write(canv_dist_n[phase]);
  };
  for(phase=0; phase<6; phase++)
  {
    canv_dist_vs_run_n[phase] = Form("canv_dist_vs_run_phase%d",phase);
    canv_dist_vs_run[phase] = new TCanvas(canv_dist_vs_run_n[phase].Data(),canv_dist_vs_run_n[phase].Data(),1000,500);
    canv_dist_vs_run[phase]->Divide(4,2);
    for(s=0; s<8; s++)
    {
      canv_dist_vs_run[phase]->cd(pads[s]);
      dist_vs_run[s][phase]->Draw("colz");
    };
    canv_dist_vs_run[phase]->Write(canv_dist_vs_run_n[phase]);
  };

  // print pngs
  for(phase=0; phase<6; phase++)
  {
    canv_dist_n[phase]+=".png";
    canv_dist[phase]->Print(canv_dist_n[phase],"png");
  };
};
