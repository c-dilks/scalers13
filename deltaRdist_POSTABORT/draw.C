void draw() {
  TTree * t = new TTree();

  t->ReadFile("fitdata.dat","npa/I:id/I:mu/F:sigma/F");
  Int_t npa; // number of postabort bXings
  Int_t id; // 0: single gaus  1: left gaus  2: right gaus  -1: bad fit
  Float_t mu,sigma; // mean and stdev of gaus
  t->SetBranchAddress("npa",&npa);
  t->SetBranchAddress("id",&id);
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("sigma",&sigma);

  TGraphErrors * gr1 = new TGraphErrors(); // 1 gaus
  TGraphErrors * gr2 = new TGraphErrors(); // 2 gaus
  TGraphErrors * grBad = new TGraphErrors(); // untrustworthy fit
  Int_t gr1C=0;
  Int_t gr2C=0;
  Int_t grBadC=0;
  for(int x=0; x<t->GetEntries(); x++) {
    t->GetEntry(x);
    sigma = fabs(sigma);
    if(id==0) {
      gr1->SetPoint(gr1C,npa,mu);
      gr1->SetPointError(gr1C,0,sigma);
      gr1C++;
    } else if(id>0) {
      gr2->SetPoint(gr2C,npa,mu);
      gr2->SetPointError(gr2C,0,sigma);
      gr2C++;
    } else {
      grBad->SetPoint(grBadC,npa,mu);
      grBad->SetPointError(grBadC,0,sigma);
      grBadC++;
    };
  };

  gr1->SetMarkerStyle(kFullCircle);
  gr2->SetMarkerStyle(kFullCircle);
  grBad->SetMarkerStyle(kOpenCircle);

  gr1->SetMarkerColor(kGreen+1);
  gr2->SetMarkerColor(kBlue);
  grBad->SetMarkerColor(kMagenta);
  
  TMultiGraph * mgr = new TMultiGraph();
  mgr->Add(gr1);
  mgr->Add(gr2);
  mgr->Add(grBad);
  mgr->SetTitle(
   "#Delta R_{3} Gaussian Fit Means vs. N_{pa} Omitted Post-Abort bXings;N_{pa};#Delta R_{3} #mu #pm #sigma");
  
  TCanvas * c = new TCanvas("c","c",1000,600);
  c->SetGrid(1,1);
  mgr->Draw("APE");
  mgr->GetYaxis()->SetTitleOffset(1.3);
};
