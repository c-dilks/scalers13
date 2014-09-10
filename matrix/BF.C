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
        Int_t numer_tbit=1, Int_t numer_cbit=2, Int_t denom_tbit=2, Int_t denom_cbit=2,
        Bool_t evaluateChiSquare=false, Int_t specificPattern=0)
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  // open hadded matx tree and check its ordering by running DrawMatrix.C
  TFile * infile = new TFile("rootfiles/all.root","READ");
  TTree * matx = (TTree*) infile->Get("matx");
  TFile * countsfile = new TFile("../counts.root","READ");
  TTree * counts = (TTree*) countsfile->Get("sca");
  TFile * sumsfile = new TFile("../sums.root","READ");
  TTree * sums_tr = (TTree*) sumsfile->Get("sum");
  //printf("executing DrawMatrix.C to check ordering of matx tree\n");
  //gROOT->ProcessLine(".x DrawMatrix.C");


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

  char detstr[16];
  printf("\n--------------------------\n");
  printf("ratio r^i = %s%s / %s%s\n",
    tbit_s[numer_tbit],cbit_s[numer_cbit],
    tbit_s[denom_tbit],cbit_s[denom_cbit]);
  sprintf(detstr,"%s%s/%s%s",
    tbit_s[numer_tbit],cbit_s[numer_cbit],
    tbit_s[denom_tbit],cbit_s[denom_cbit]);
  printf("--------------------------\n\n");

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


  // check to see if matx tree is sorted by run index from hadd
  // -- if it's not, you'll need to write a sorting routine
  printf("checking matx...\n");
  Int_t index_tmp=0;
  matx->SetBranchAddress("i",&i);
  matx->SetBranchAddress("fi",&fi);
  for(Int_t ii=0; ii<matx->GetEntries(); ii++)
  {
    matx->GetEntry(ii);
    if(i != index_tmp)
    {
      if(i == index_tmp+1) index_tmp = i;
      else 
      {
        fprintf(stderr,"ERROR: hadded matx tree not correctly sorted\n");
        return;
      };
    };
  };


  // define ratio array "rat" and corresponding uncertainties "unc"
  // -- from documentation:  r^i = rat;  sigma_{r^i} = unc
  Int_t imax_tmp = matx->GetMaximum("i");
  const Int_t imax = imax_tmp;
  Double_t rat[imax][120]; // rat = num / den
  Double_t den[imax][120]; // denominator
  Double_t num[imax][120]; // numerator
  Double_t unc[imax][120]; // variances of rat
  printf("computing r^i for each run and bXing...\n");
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
  char rat_v_bx_n[imax][32];
  char rat_v_bx_t[imax][128];
  for(Int_t ii=0; ii<imax; ii++)
  {
    sprintf(rat_v_bx_n[ii],"rat_v_bx_run%d",ii);
    sprintf(rat_v_bx_t[ii],"%s vs. bXing for i=%d",detstr,ii);
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
    kicked_arr[i_c-1][bx_c]=kicked;
  };


  // set sums tree branch address and fill pattern number arrays
  Int_t i_s,pattern;
  Int_t pattern_arr[imax];
  sums_tr->SetBranchAddress("i",&i_s);
  sums_tr->SetBranchAddress("pattern",&pattern);
  for(Int_t ii=0; ii<sums_tr->GetEntries(); ii++)
  {
    sums_tr->GetEntry(ii);
    pattern_arr[i_s-1]=pattern;
  };

  
  // fill spin array H 
  // NOTE: H is automatically set to zero if kicked==true;
  //       and if H==0 for a bXing, it's not included in sigma functions
  //
  // -- see equ/bunch_fitting/H_table.pdf for H_a definitions
  //
  Int_t H[imax][10][120]; // [run] [asym] [bx]
  Int_t hb,hy;
  for(Int_t ii=0; ii<imax; ii++)
  {
    for(Int_t bb=0; bb<120; bb++)
    {
      hb = blue_arr[ii][bb];
      hy = yell_arr[ii][bb];
      if(kicked_arr[ii][bb]==0 && abs(hb*hy)==1)
      {
        H[ii][1][bb] = hy;
        H[ii][2][bb] = hb;
        H[ii][3][bb] = hb * hy;
        H[ii][4][bb] = (hb + hy) / 2;
        H[ii][5][bb] = (1 - hb) * hy / 2;
        H[ii][6][bb] = (1 - hy) * hb / 2;
        H[ii][7][bb] = (1 + hb) * hy / 2;
        H[ii][8][bb] = (hy - hb) / 2;
        H[ii][9][bb] = (1 + hy) * hb / 2;
      }
      else for(Int_t aa=1; aa<10; aa++) H[ii][aa][bb]=0;
    };
  };


  // H vs. bXing
  TH1D * H_v_bx[10][imax]; // [asym] [run index]
  TObjArray * H_v_bx_arr[10];
  char H_v_bx_n[10][imax][32];
  char H_v_bx_t[10][imax][128];
  for(Int_t aa=1; aa<10; aa++)
  {
    H_v_bx_arr[aa] = new TObjArray();
    for(Int_t ii=0; ii<imax; ii++)
    {
      sprintf(H_v_bx_n[aa][ii],"H%d_v_bx_%d",aa,ii);
      sprintf(H_v_bx_t[aa][ii],"H%d vs. bXing for i=%d",aa,ii);
      H_v_bx[aa][ii] = new TH1D(H_v_bx_n[aa][ii],H_v_bx_t[aa][ii],120,0,120);
      for(Int_t bb=0; bb<120; bb++)
      {
        H_v_bx[aa][ii]->SetBinContent(bb-1,H[ii][aa][bb]);
      };
      H_v_bx_arr[aa]->AddLast(H_v_bx[aa][ii]);
    };
  };


  // fill Hsum arrays (sum of H over bXings for each run & asymmetry number)
  Int_t Hsum[10][imax]; // [asym] [run]
  for(Int_t aa=1; aa<10; aa++)
  {
    for(Int_t ii=0; ii<imax; ii++)
    {
      Hsum[aa][ii] = 0;
      for(Int_t bb=0; bb<120; bb++)
      {
        if(kicked_arr[ii][bb]==0) Hsum[aa][ii] += H[ii][aa][bb];
      };
    };
  };


  // Hsum vs. run index
  TH1D * Hsum_dist[10]; //asym
  char Hsum_dist_n[10][32];
  char Hsum_dist_t[10][128];
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(Hsum_dist_n[aa],"Hsum_a%d",aa);
    sprintf(Hsum_dist_t[aa],"Sum of H_{%d} over filled bXings vs. run index",aa);
    Hsum_dist[aa] = new TH1D(Hsum_dist_n[aa],Hsum_dist_t[aa],imax,1,imax+1);
    for(Int_t ii=0; ii<imax; ii++) Hsum_dist[aa]->SetBinContent(ii+1,Hsum[aa][ii]);
  };


  // compute sigma functions for each run, summing over bXings
  // -- if specificPattern!=0, only let sigma function be nonzero for fills with pattern
  //    no. equal to specificPattern
  Double_t sigma_id[10][imax]; // [asym] [run]
  Double_t sigma_r[10][imax];
  Double_t sigma_h[10][imax];
  Double_t sigma_hr[10][imax];
  printf("computing sigma functions...\n");
  if(specificPattern!=0) printf(" [+] specificPattern=%d\n",specificPattern);
  for(Int_t aa=1; aa<10; aa++)
  {
    for(Int_t ii=0; ii<imax; ii++)
    {
      sigma_id[aa][ii]=0;
      sigma_r[aa][ii]=0;
      sigma_h[aa][ii]=0;
      sigma_hr[aa][ii]=0;
      if(specificPattern==0 || pattern_arr[ii]==specificPattern)
      {
        for(Int_t bb=0; bb<120; bb++)
        {
          if(kicked_arr[ii][bb]==0 && H[ii][aa][bb]!=0 && unc[ii][bb]>0)
          {
            sigma_id[aa][ii] += 1 / pow(unc[ii][bb],2);
            sigma_r[aa][ii] += rat[ii][bb] / pow(unc[ii][bb],2);
            sigma_h[aa][ii] += H[ii][aa][bb] / pow(unc[ii][bb],2);
            sigma_hr[aa][ii] += H[ii][aa][bb] * rat[ii][bb] / pow(unc[ii][bb],2);
          };
        };
      };
    };
  };


  // compute sigma_sum functions, which are sums of sigma functions over runs
  Double_t sigma_sum_id[10];
  Double_t sigma_sum_r[10];
  Double_t sigma_sum_h[10];
  Double_t sigma_sum_hr[10];
  for(Int_t aa=1; aa<10; aa++)
  {
    sigma_sum_id[aa] = 0;
    sigma_sum_r[aa] = 0;
    sigma_sum_h[aa] = 0;
    sigma_sum_hr[aa] = 0;
    for(Int_t ii=0; ii<imax; ii++)
    {
      sigma_sum_id[aa] += sigma_id[aa][ii];
      sigma_sum_r[aa] += sigma_r[aa][ii];
      sigma_sum_h[aa] += sigma_h[aa][ii];
      sigma_sum_hr[aa] += sigma_hr[aa][ii];
    };
  };


  // sigma function vs. run
  TH1D * sigma_id_dist[10]; // [asym]
  TH1D * sigma_r_dist[10];
  TH1D * sigma_h_dist[10];
  TH1D * sigma_hr_dist[10];
  char sigma_id_dist_n[10][16];
  char sigma_r_dist_n[10][16];
  char sigma_h_dist_n[10][16];
  char sigma_hr_dist_n[10][16];
  char sigma_id_dist_t[10][128];
  char sigma_r_dist_t[10][128];
  char sigma_h_dist_t[10][128];
  char sigma_hr_dist_t[10][128];
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(sigma_id_dist_n[aa],"sigma_id_a%d",aa);
    sprintf(sigma_r_dist_n[aa],"sigma_r_a%d",aa);
    sprintf(sigma_h_dist_n[aa],"sigma_h_a%d",aa);
    sprintf(sigma_hr_dist_n[aa],"sigma_hr_a%d",aa);
    sprintf(sigma_id_dist_t[aa],"#Sigma(1) vs. run // asym=%d // r=%s",aa,detstr);
    sprintf(sigma_r_dist_t[aa],"#Sigma(r) vs. run // asym=%d // r=%s",aa,detstr);
    sprintf(sigma_h_dist_t[aa],"#Sigma(H) vs. run // asym==%d // r=%s",aa,detstr);
    sprintf(sigma_hr_dist_t[aa],"#Sigma(Hr) vs. run // asym==%d // r=%s",aa,detstr);
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


  // constant vs. run index and epsilon vs. run index
  // -- if specificPattern!=0, only compute constant & epsilon for fills with pattern
  //    no. equal to specificPattern
  TH1D * cons_dist[10]; // [asym]
  char cons_dist_n[10][16];
  char cons_dist_t[10][64];
  TH1D * epsi_dist[10]; // [asym]
  char epsi_dist_n[10][16];
  char epsi_dist_t[10][64];
  Double_t Sid,Sr,Sh,Shr,cons,epsi;
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(cons_dist_n[aa],"cons_a%d_v_run",aa);
    sprintf(epsi_dist_n[aa],"epsi_a%d_v_run",aa);
    sprintf(cons_dist_t[aa],"c_{%d} vs. run index // r=%s",aa,detstr);
    sprintf(epsi_dist_t[aa],"#varepsilon_{%d} vs. run index // r=%s",aa,detstr);
    cons_dist[aa] = new TH1D(cons_dist_n[aa],cons_dist_t[aa],imax,1,imax+1);
    epsi_dist[aa] = new TH1D(epsi_dist_n[aa],epsi_dist_t[aa],imax,1,imax+1);
    for(Int_t ii=0; ii<imax; ii++)
    {
      if(specificPattern==0 || pattern_arr[ii]==specificPattern)
      {
        Sid = sigma_id[aa][ii];
        Sr = sigma_r[aa][ii];
        Sh = sigma_h[aa][ii];
        Shr = sigma_hr[aa][ii];
        cons = ( Sh * Shr - Sid * Sr) / ( pow(Sh,2) - pow(Sid,2) );
        epsi = ( Sh * Sr - Sid * Shr ) / ( Sh * Shr - Sid * Sr );
        cons_dist[aa]->SetBinContent(ii+1,cons);
        epsi_dist[aa]->SetBinContent(ii+1,epsi);
      };
    };
  };


  // compute statistical errors on constant and epsilon, vs. run index
  Double_t cons_err,epsi_err;
  for(Int_t aa=1; aa<10; aa++)
  {
    for(Int_t ii=0; ii<imax; ii++)
    {
      cons_err = 0; 
      epsi_err = 0;
      if(specificPattern==0 || pattern_arr[ii]==specificPattern)
      {
        // first sum the squares...
        for(Int_t bb=0; bb<120; bb++)
        {
          if(H[ii][aa][bb]!=0 && unc[ii][bb]>0)
          {
            epsi_err += pow(1/unc[ii][bb] * 
              ( ( pow(sigma_h[aa][ii],2) - pow(sigma_id[aa][ii],2) ) * 
                ( sigma_hr[aa][ii] - H[ii][aa][bb]*sigma_r[aa][ii] ) ) / 
              ( pow(sigma_hr[aa][ii]*sigma_h[aa][ii] - sigma_id[aa][ii]*sigma_r[aa][ii],2) ),2);
            cons_err += pow(1/unc[ii][bb] *
              ( H[ii][aa][bb]*sigma_h[aa][ii] - sigma_id[aa][ii] ) /
              ( pow(sigma_h[aa][ii],2) - pow(sigma_id[aa][ii],2) ),2);
          };
        };
        // ... then take the square root of the sum of squares
        cons_err = sqrt(cons_err);
        epsi_err = sqrt(epsi_err);
        cons_dist[aa]->SetBinError(ii+1,cons_err);
        epsi_dist[aa]->SetBinError(ii+1,epsi_err);
      };
    };
  };
      



  // binary search algorithm for parameter range such that
  Double_t optimal_cons[10]; // [asym]  // constant fit of cons vs. i 
  Double_t optimal_epsi[10]; // [asym]  // constant fit of epsi vs. i
  Double_t error_cons[10]; // [asym] // constant fit error of cons vs. i
  Double_t error_epsi[10]; // [asym] // constant fit error of epsi vs. i
  Double_t lb,ub; // binary search lower & upper bounds
  Double_t cx,ex,chi2_eval;
  for(Int_t aa=1; aa<10; aa++)
  {
    // don't do the constant fits if specificPattern!=0
    if(specificPattern==0)
    {
      cons_dist[aa]->Fit("pol0","Q","",1,imax+1);
      epsi_dist[aa]->Fit("pol0","Q","",1,imax+1);
      //optimal_cons[aa] = cons_dist[aa]->GetFunction("pol0")->GetParameter(0);
      //optimal_epsi[aa] = epsi_dist[aa]->GetFunction("pol0")->GetParameter(0);
      error_cons[aa] = cons_dist[aa]->GetFunction("pol0")->GetParError(0);
      error_epsi[aa] = epsi_dist[aa]->GetFunction("pol0")->GetParError(0);
    }
    else
    {
      error_cons[aa] = 1;  // set "fake error" if specificPattern!=0
      error_epsi[aa] = 1;
    };

    // determine overall chi2 minimum for all bXings and runs
    Sid = sigma_sum_id[aa];
    Sr = sigma_sum_r[aa];
    Sh = sigma_sum_h[aa];
    Shr = sigma_sum_hr[aa];
    optimal_cons[aa] = ( Sh * Shr - Sid * Sr) / ( pow(Sh,2) - pow(Sid,2) );
    optimal_epsi[aa] = ( Sh * Sr - Sid * Shr ) / ( Sh * Shr - Sid * Sr );
  };
  // determine minimum chi2 value
  Double_t chi2_min[10];
  for(Int_t aa=1; aa<10; aa++)
  {
    chi2_min[aa] = CalcChi2(aa,optimal_cons[aa],optimal_epsi[aa],imax,H,rat,unc);
  };
  Double_t ub_cons[10];  // upper / lower bounds for binary search 
  Double_t ub_epsi[10];
  Double_t lb_cons[10];
  Double_t lb_epsi[10];
  Double_t max_cons[10]; // max range of cons
  Double_t max_epsi[10]; // max range of epsi
  Double_t range_cons[10]; // max_cons - optimal_cons
  Double_t range_epsi[10]; // max_epsi - optimal_epsi
  Double_t max_delta = 5; // chi2 @ chi2 profile max range = chi2_min +/- max_delta
  Double_t max_delta_pad = 0.001; // how close chi2_check-chi2_min is to max_delta in
                                //   order to halt binary search
  Double_t chi2_check;
  Int_t NN;
  Double_t mp;
  if(evaluateChiSquare)
  {
    for(Int_t aa=1; aa<10; aa++)
    {
      // cons upper bound
      chi2_check=0;
      NN = 0;
      while(chi2_check - chi2_min[aa] < max_delta)
      {
        printf("[+] chi2_check-chi2_min = %f, max_delta = %f\n",
          chi2_check-chi2_min[aa],max_delta);
        NN++;
        chi2_check = CalcChi2(aa,
                              optimal_cons[aa] + NN * error_cons[aa],
                              optimal_epsi[aa],
                              imax,H,rat,unc);
      };
      printf("[+] chi2_check-chi2_min = %f, max_delta = %f\n",
        chi2_check-chi2_min[aa],max_delta);
      ub_cons[aa] = optimal_cons[aa] + NN * error_cons[aa];

      // cons binary search
      lb_cons[aa] = optimal_cons[aa];
      while(fabs(max_delta - (chi2_check - chi2_min[aa])) >= max_delta_pad)
      {
        mp = lb_cons[aa] + (ub_cons[aa] - lb_cons[aa])/2;
        chi2_check = CalcChi2(aa,mp,optimal_epsi[aa],imax,H,rat,unc);
        printf(" --> | max_delta - (chi2_check-chi2_min) | = %f\n",
          fabs(max_delta - (chi2_check-chi2_min[aa])));
        if(chi2_check - chi2_min[aa] > max_delta) ub_cons[aa] = mp;
        else lb_cons[aa] = mp;
      };
      max_cons[aa] = mp;
      printf(" ------> chi2_check-chi2_min = %f\n",chi2_check-chi2_min[aa]);
      printf("         optimal_cons=%f max_cons=%f\n",optimal_cons[aa],max_cons[aa]);
      

      // epsi upper bound
      chi2_check=0;
      NN = 0;
      while(chi2_check - chi2_min[aa] < max_delta)
      {
        NN++;
        chi2_check = CalcChi2(aa,
                              optimal_cons[aa],
                              optimal_epsi[aa] + NN * error_epsi[aa],
                              imax,H,rat,unc);
      };
      ub_epsi[aa] = optimal_epsi[aa] + NN * error_epsi[aa];

      // epsi binary search
      lb_epsi[aa] = optimal_epsi[aa];
      while(fabs(max_delta - (chi2_check - chi2_min[aa])) >= max_delta_pad)
      {
        mp = lb_epsi[aa] + (ub_epsi[aa] - lb_epsi[aa])/2;
        chi2_check = CalcChi2(aa,optimal_cons[aa],mp,imax,H,rat,unc);
        printf(" --> | max_delta - (chi2_check-chi2_min) | = %f\n",
          fabs(max_delta - (chi2_check-chi2_min[aa])));
        if(chi2_check - chi2_min[aa] > max_delta) ub_epsi[aa] = mp;
        else lb_epsi[aa] = mp;
      };
      max_epsi[aa] = mp;
      printf(" ------> chi2_check-chi2_min = %f\n",chi2_check-chi2_min[aa]);
      printf("         optimal_epsi=%f max_epsi=%f\n",optimal_epsi[aa],max_epsi[aa]);

      range_cons[aa] = fabs(max_cons[aa] - optimal_cons[aa]);
      range_epsi[aa] = fabs(max_epsi[aa] - optimal_epsi[aa]);

      printf("[***] max_cons - optimal_cons = %.15f\n",range_cons[aa]);
      printf("[***] max_epsi - optimal_epsi = %.15f\n",range_epsi[aa]);

    };
  };


  // evaluate chi2 by varying fit parameters
  const Int_t NC = 15; // number of bins around optimal to evaluate chi2
  TH1D * chi2_vs_cons[10];
  char chi2_vs_cons_n[10][32];
  char chi2_vs_cons_t[10][128];
  TH1D * chi2_vs_epsi[10];
  char chi2_vs_epsi_n[10][32];
  char chi2_vs_epsi_t[10][128];
  TH2D * chi2_vs_both[10];
  char chi2_vs_both_n[10][32];
  char chi2_vs_both_t[10][128];
  if(evaluateChiSquare)
  {
    for(Int_t aa=1; aa<10; aa++)
    {
      printf("evaluating chi2 profile for asym %d...\n",aa);
      sprintf(chi2_vs_cons_n[aa],"chi2_vs_cons_a%d",aa);
      sprintf(chi2_vs_epsi_n[aa],"chi2_vs_epsi_a%d",aa);
      sprintf(chi2_vs_both_n[aa],"chi2_vs_both_a%d",aa);
      sprintf(chi2_vs_cons_t[aa],
        "#chi^{2}_{%d} vs. c_{%d} // #varepsilon_{%d}=%f // r=%s;c_{%d}",
        aa,aa,aa,optimal_epsi[aa],detstr,aa);
      sprintf(chi2_vs_epsi_t[aa],
        "#chi^{2}_{%d} vs. #varepsilon_{%d} // c_{%d}=%f // r=%s;#varepsilon_{%d}",
        aa,aa,aa,optimal_cons[aa],detstr,aa);
      sprintf(chi2_vs_both_t[aa],
        "#chi^{2}_{%d} vs. c_{%d} vs. #varepsilon_{%d} // r=%s;#varepsilon_{%d};c_{%d}",
        aa,aa,aa,detstr,aa,aa);
      chi2_vs_cons[aa] = new TH1D(chi2_vs_cons_n[aa],chi2_vs_cons_t[aa],
        NC,optimal_cons[aa]-range_cons[aa],optimal_cons[aa]+range_cons[aa]);
      chi2_vs_epsi[aa] = new TH1D(chi2_vs_epsi_n[aa],chi2_vs_epsi_t[aa],
        NC,optimal_epsi[aa]-range_epsi[aa],optimal_epsi[aa]+range_epsi[aa]);
      chi2_vs_both[aa] = new TH2D(chi2_vs_both_n[aa],chi2_vs_both_t[aa],
        NC,optimal_epsi[aa]-range_epsi[aa],optimal_epsi[aa]+range_epsi[aa],
        NC,optimal_cons[aa]-range_cons[aa],optimal_cons[aa]+range_cons[aa]);


      // vary cons holding epsi fixed
      printf(" stage 1/3\n");
      ex = optimal_epsi[aa];
      for(Int_t bin=1; bin<=NC; bin++)
      {
        cx = chi2_vs_cons[aa]->GetBinCenter(bin);
        chi2_eval = CalcChi2(aa,cx,ex,imax,H,rat,unc);
        printf(" -- chi2_eval-chi2_min=%f\n",chi2_eval-chi2_min[aa]);
        chi2_vs_cons[aa]->SetBinContent(bin,chi2_eval);
      };

      // vary epsi holding cons fixed
      printf(" stage 2/3\n");
      cx = optimal_cons[aa];
      for(Int_t bin=1; bin<=NC; bin++)
      {
        ex = chi2_vs_epsi[aa]->GetBinCenter(bin);
        chi2_eval = CalcChi2(aa,cx,ex,imax,H,rat,unc);
        printf(" -- chi2_eval-chi2_min=%f\n",chi2_eval-chi2_min[aa]);
        chi2_vs_epsi[aa]->SetBinContent(bin,chi2_eval);
      };

      // vary both cons and epsi
      printf(" stage 3/3\n");
      for(Int_t binx=1; binx<=NC; binx++)
      {
        for(Int_t biny=1; biny<=NC; biny++)
        {
          ex = chi2_vs_both[aa]->GetXaxis()->GetBinCenter(binx);
          cx = chi2_vs_both[aa]->GetYaxis()->GetBinCenter(biny);
          chi2_eval = CalcChi2(aa,cx,ex,imax,H,rat,unc);
          printf(" -- chi2_eval-chi2_min=%f\n",chi2_eval-chi2_min[aa]);
          chi2_vs_both[aa]->SetBinContent(binx,biny,chi2_eval);
        };
      };
    };
  };
      


  // write output
  const char outfile_name[64];
  if(specificPattern==0)
    sprintf(outfile_name,"fit_result.%s%s.%s%s.root",
      tbit_s[numer_tbit],cbit_s[numer_cbit],
      tbit_s[denom_tbit],cbit_s[denom_cbit]);
  else
    sprintf(outfile_name,"pats/fit_result.%s%s.%s%s.pat%d.root",
      tbit_s[numer_tbit],cbit_s[numer_cbit],
      tbit_s[denom_tbit],cbit_s[denom_cbit],specificPattern);
  TFile * outfile = new TFile(outfile_name,"RECREATE");
  rat_v_bx_arr->Write("rat_v_bx_arr",TObject::kSingleKey);

  outfile->mkdir("H_v_bx_arrs");
  outfile->cd("/H_v_bx_arrs");
  char H_v_bx_arr_n[10][32];
  for(Int_t aa=1; aa<10; aa++)
  {
    sprintf(H_v_bx_arr_n[aa],"H_a%d_v_bx_arr",aa);
    H_v_bx_arr[aa]->Write(H_v_bx_arr_n[aa],TObject::kSingleKey);
  };

  outfile->mkdir("Hsum");
  outfile->mkdir("sigma_id");
  outfile->mkdir("sigma_r");
  outfile->mkdir("sigma_h");
  outfile->mkdir("sigma_hr");
  outfile->mkdir("constant");
  outfile->mkdir("epsilon");
  for(Int_t aa=1; aa<10; aa++) 
  {
    outfile->cd("/Hsum"); Hsum_dist[aa]->Write();
    outfile->cd("/sigma_id"); sigma_id_dist[aa]->Write();
    outfile->cd("/sigma_r"); sigma_r_dist[aa]->Write();
    outfile->cd("/sigma_h"); sigma_h_dist[aa]->Write();
    outfile->cd("/sigma_hr"); sigma_hr_dist[aa]->Write();
    outfile->cd("/constant"); cons_dist[aa]->Write();
    outfile->cd("/epsilon"); epsi_dist[aa]->Write();
  };
  if(evaluateChiSquare)
  {
    outfile->mkdir("chi2_profile_constant");
    outfile->mkdir("chi2_profile_epsilon");
    outfile->mkdir("chi2_profile_2d");
    for(Int_t aa=1; aa<10; aa++) 
    {
      outfile->cd("chi2_profile_constant"); chi2_vs_cons[aa]->Write();
      outfile->cd("chi2_profile_epsilon"); chi2_vs_epsi[aa]->Write();
      outfile->cd("chi2_profile_2d"); chi2_vs_both[aa]->Write();
    };
  };
  printf("%s created\n",outfile_name);
}


//Double_t CalcChi2(Double_t cx0, Double_t ex0, 
                  //Int_t imax0, Int_t * H0, Double_t * rat0, Double_t * unc0)
Double_t CalcChi2(Int_t aa0, Double_t cx0, Double_t ex0, 
                  Int_t imax0, 
                  Int_t H0[][10][120], 
                  Double_t rat0[][120], 
                  Double_t unc0[][120])
{
  Double_t value = 0;
  for(Int_t ii=0; ii<imax0; ii++)
  {
    for(Int_t bb=0; bb<120; bb++)
    {
      if(H0[aa0][ii][bb]!=0 && unc0[ii][bb]>0)
      {
        value += pow((cx0*(1+H0[aa0][ii][bb]*ex0)-rat0[ii][bb])/unc0[ii][bb], 2);
      };
    };
  };
  return value;
};
