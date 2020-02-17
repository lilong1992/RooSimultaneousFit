using namespace RooFit;

void gefittingfeb(TString filename){
  TRandom3 RA(0);
  TChain* tree = new TChain( "gerecoilTree" );
  Int_t nFilesAdded = tree->Add( filename );

  Int_t bdID;
  Double_t sigrec;
  Bool_t ifsig;

  tree->SetBranchAddress("bdID", &bdID);
  tree->SetBranchAddress("sigrec", &sigrec);
  tree->SetBranchAddress("ifsig", &ifsig);

  RooRealVar ph("ph","peakhigh/adc",1050.5,1850.5) ;
  RooDataSet* ds[27];
  TH1F* bghist[27];
  TH1F* meanhist[27];
  auto bghistin=new TH1F("bghistin","bghistin",5000,0.5,5000.5);
  auto meanhistout=new TH1F("meanhistout","meanhistout;keVee;AU",1300,0,1300);
  auto meanhistin=new TH1F("meanhistin","meanhistin;keVee;AU",1300,0,1300);
  for(int i=0;i<27;i++){
    TString s="ds";
    s.Append(to_string(i));
    ds[i]=new RooDataSet(s,s,RooArgSet(ph));
    s.Append("h");
    bghist[i]=new TH1F(s,s,1000,0.5,5000.5);
    s.Append("m");
    meanhist[i]=new TH1F(s,s,1300,0,1300);
  }

  Long64_t nEntries = tree->GetEntries();
  for(Long64_t entryNumber=0;entryNumber<nEntries;entryNumber++){
    tree->GetEntry(entryNumber);
    if (ifsig) {
      ph=sigrec;
      ds[bdID]->add(RooArgSet(ph));
    }
    else {
      //bghist[bdID]->Fill(sigrec);
      if (bdID<8){
        bghistin->Fill(sigrec);
        for (int j=0;j<8;j++) bghist[j]->Fill(sigrec);
      }
    }
  }

  RooRealVar sigmaout("sigmaout","width of gaussian",60,50,85) ;
  RooRealVar sigmain("sigmain","width of gaussian",35,30,55) ;
  RooCategory sampleout("sampleout","sampleout");
  RooCategory samplein("samplein","samplein");
  RooSimultaneous simPdfout("simPdfout","simultaneous pdf out",sampleout) ;
  RooSimultaneous simPdfin("simPdfin","simultaneous pdf in",samplein) ;

  RooRealVar* mean[27];
  RooGaussian* gs[27];
  RooUniform* uni[27];//uniform bg
  RooAddPdf* model[27];
  RooDataHist* bgh[27];
  RooHistPdf* bghpdf[27];//bghist to pdf
  RooRealVar* sigcounts[27];
  RooRealVar* bgcounts[27];
  RooHist* hresid[27];
  for (int i=0;i<27;i++){
    TString s1="mean";
    s1.Append(to_string(i));
    mean[i]=new RooRealVar(s1,s1,1450,1200,1700);
    TString s2="gauss";
    s2.Append(to_string(i));
    TString s3="det";
    s3.Append(to_string(i));
    TString s4="uni";
    s4.Append(to_string(i));
    uni[i]=new RooUniform(s4,s4,ph);
    TString s5="model";
    s5.Append(to_string(i));
    TString s6="sigcounts";
    s6.Append(to_string(i));
    sigcounts[i]=new RooRealVar(s6,s6,60,0,300);
    TString s7="bgcounts";
    s7.Append(to_string(i));
    bgcounts[i]=new RooRealVar(s7,s7,1,0,30000);
    TString s8="bgh";
    s8.Append(to_string(i));
    bgh[i]=new RooDataHist(s8,s8,ph,Import(*bghist[i]));
    s8.Append("pdf");
    bghpdf[i]=new RooHistPdf(s8,s8,ph,*bgh[i],1);

    if (i<8) {
      gs[i]=new RooGaussian(s2,s2,ph,*mean[i],sigmain);
      model[i]=new RooAddPdf(s5,s5,RooArgList(*gs[i],*bghpdf[i]),RooArgList(*sigcounts[i],*bgcounts[i]));
      samplein.defineType(s3);
      simPdfin.addPdf(*model[i],s3);
    }
    else {
      gs[i]=new RooGaussian(s2,s2,ph,*mean[i],sigmaout);
      model[i]=new RooAddPdf(s5,s5,RooArgList(*gs[i],*uni[i]),RooArgList(*sigcounts[i],*bgcounts[i]));
      sampleout.defineType(s3);
      simPdfout.addPdf(*model[i],s3);
    }
  }

  RooDataSet combDatain("combDatain","combined data in",ph,Index(samplein),Import("det0",*ds[0]));
  RooDataSet combDataout("combDataout","combined data out",ph,Index(sampleout),Import("det8",*ds[8]));
  for (int i=1;i<8;i++){
    TString s3="det";
    s3.Append(to_string(i));
    auto tmpds=new RooDataSet("tmpds","tmpds",ph,Index(samplein),Import(s3,*ds[i]));
    combDatain.append(*tmpds);
  }
  for (int i=9;i<27;i++){
    TString s3="det";
    s3.Append(to_string(i));
    auto tmpds=new RooDataSet("tmpds","tmpds",ph,Index(sampleout),Import(s3,*ds[i]));
    combDataout.append(*tmpds);
  }

  simPdfin.fitTo(combDatain,Range(1170,1500),NumCPU(6));
  auto c1=new TCanvas("c1","inner ring fitting",1600,900);
  c1->Divide(4,2);
  for(Int_t i=0;i<8;i++){
    c1->cd(i+1);
    //ph.setRange("signal",1251,1450) ;
    TString det=to_string(i+5); det.Append("th det");
    TString s3="det"; s3.Append(to_string(i));
    TString s8="bgh"; s8.Append(to_string(i)); s8.Append("pdf");
    RooPlot* frame = ph.frame(Title(det),Range(1150,1550)) ;
    //ds[i]->plotOn(frame,Binning(70));
    combDatain.plotOn(frame,Cut("samplein==samplein::"+s3),Binning(80)) ;
    simPdfin.plotOn(frame,Slice(samplein,s3),ProjWData(samplein,combDatain)) ;
    simPdfin.plotOn(frame,Slice(samplein,s3),Components(s8),ProjWData(samplein,combDatain),LineStyle(kDashed)) ;
    frame->GetYaxis()->SetRangeUser(-5,30) ;
    hresid[i] = frame->residHist() ;
    frame->Draw();
  }

  simPdfout.fitTo(combDataout,Range(1350,1800),NumCPU(6));
  auto c2=new TCanvas("c2","outer ring fitting",1600,900);
  c2->Divide(5,4);
  for(Int_t i=8;i<27;i++){
    c2->cd(i-7);
    //ph.setRange("signal",1350,1800) ;
    TString det=to_string(i+5); det.Append("th det");
    TString s3="det";    s3.Append(to_string(i));
    RooPlot* frame = ph.frame(Title(det),Range(1300,1850)) ;
    //ds[i]->plotOn(frame,Binning(35));
    combDataout.plotOn(frame,Cut("sampleout==sampleout::"+s3),Binning(40)) ;
    simPdfout.plotOn(frame,Slice(sampleout,s3),ProjWData(sampleout,combDataout)) ;

    frame->GetYaxis()->SetRangeUser(-5,30) ;
    hresid[i] = frame->residHist() ;
    frame->Draw();
  }

  Double_t anglesout[19]={165.9,180,194.1,30.9,45,59.1,75.9,90,104.1,255.9,270,284.1,300.9,315,329.1,345.9,0,14.1,149.1};
  Double_t meanout[19]={};
  Double_t meanerrout[19]={};
  Double_t countout[19]={};
  Double_t counterrout[19]={};
  Double_t anglesin[8]={45,90,135,180,225,270,315,0};
  Double_t meanin[8]={};
  Double_t meanerrin[8]={};
  Double_t countin[8]={};
  Double_t counterrin[8]={};
  Double_t totalcountsin=0,totalcountsout=0;

  for(int i=0;i<19;i++){
    meanout[i]=mean[i+8]->getValV();
    meanerrout[i]=mean[i+8]->getError();
    meanout[i]=meanout[i]*2.178-2513.5;
    meanerrout[i]=meanerrout[i]*2.178;
    countout[i]=sigcounts[i+8]->getValV();
    totalcountsout+=sigcounts[i+8]->getValV();
    counterrout[i]=sigcounts[i+8]->getError();
    for(int j=0;j<100000;j++) meanhist[i+8]->Fill(RA.Gaus(meanout[i],meanerrout[i]));
    meanhistout->Add(meanhist[i+8]);
  }
  for(int i=0;i<8;i++){
    meanin[i]=mean[i]->getValV();
    meanerrin[i]=mean[i]->getError();
    meanin[i]=meanin[i]*2.178-2513.5;
    meanerrin[i]=meanerrin[i]*2.178;
    countin[i]=sigcounts[i]->getValV();
    totalcountsin+=sigcounts[i]->getValV();
    counterrin[i]=sigcounts[i]->getError();
    for(int j=0;j<100000;j++) meanhist[i]->Fill(RA.Gaus(meanin[i],meanerrin[i]));
    meanhistin->Add(meanhist[i]);
  }

  auto c3=new TCanvas("c3","fitting results",1600,900);
  c3->Divide(2,2);
  c3->cd(1);
  auto gr=new TGraphErrors(19,anglesout,meanout,0,meanerrout);
  gr->SetTitle("Ge recoil energy at E_{nr}=4.93keV");
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(20);
  gr->GetXaxis()->SetTitle("Phi/degree");
  gr->GetYaxis()->SetTitle("Recoil energy/eVee");
  gr->GetYaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetLabelSize(0.03);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetRangeUser(900,1150);
  gr->GetYaxis()->SetTitleOffset(0.69);
  gr->GetXaxis()->SetLimits(-10,360);
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetXaxis()->SetLabelSize(0.03);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetNdivisions(515);
  gr->Draw("AP");
  c3->cd(2);
  auto gr2=new TGraphErrors(19,anglesout,countout,0,counterrout);
  gr2->SetTitle("Ge recoil counts vs angle Phi (outer)");
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(20);
  gr2->GetXaxis()->SetTitle("Phi/degree");
  gr2->GetYaxis()->SetTitle("Counts");
  gr2->GetYaxis()->SetTitleSize(0.04);
  gr2->GetYaxis()->SetLabelSize(0.03);
  gr2->GetYaxis()->CenterTitle();
  gr2->GetYaxis()->SetRangeUser(0,120);
  gr2->GetYaxis()->SetTitleOffset(0.69);
  gr2->GetXaxis()->SetLimits(-10,360);
  gr2->GetXaxis()->SetTitleSize(0.04);
  gr2->GetXaxis()->SetLabelSize(0.03);
  gr2->GetXaxis()->CenterTitle();
  gr2->GetXaxis()->SetNdivisions(515);
  gr2->Draw("AP");
  c3->cd(3);
  auto gri=new TGraphErrors(8,anglesin,meanin,0,meanerrin);
  gri->SetTitle("Ge recoil energy at E_{nr}=2.35keV");
  gri->SetMarkerColor(kRed);
  gri->SetMarkerStyle(20);
  gri->GetXaxis()->SetTitle("Phi/degree");
  gri->GetYaxis()->SetTitle("Recoil energy/eVee");
  gri->GetYaxis()->SetTitleSize(0.04);
  gri->GetYaxis()->SetLabelSize(0.03);
  gri->GetYaxis()->CenterTitle();
  gri->GetYaxis()->SetRangeUser(340,540);
  gri->GetYaxis()->SetTitleOffset(0.69);
  gri->GetXaxis()->SetLimits(-10,360);
  gri->GetXaxis()->SetTitleSize(0.04);
  gri->GetXaxis()->SetLabelSize(0.03);
  gri->GetXaxis()->CenterTitle();
  gri->GetXaxis()->SetNdivisions(515);
  gri->Draw("AP");
  c3->cd(4);
  auto gr2i=new TGraphErrors(8,anglesin,countin,0,counterrin);
  gr2i->SetTitle("Ge recoil counts vs angle Phi (inner)");
  gr2i->SetMarkerColor(kRed);
  gr2i->SetMarkerStyle(20);
  gr2i->GetXaxis()->SetTitle("Phi/degree");
  gr2i->GetYaxis()->SetTitle("Counts");
  gr2i->GetYaxis()->SetTitleSize(0.04);
  gr2i->GetYaxis()->SetLabelSize(0.03);
  gr2i->GetYaxis()->CenterTitle();
  gr2i->GetYaxis()->SetRangeUser(0,150);
  gr2i->GetYaxis()->SetTitleOffset(0.69);
  gr2i->GetXaxis()->SetLimits(-10,360);
  gr2i->GetXaxis()->SetTitleSize(0.04);
  gr2i->GetXaxis()->SetLabelSize(0.03);
  gr2i->GetXaxis()->CenterTitle();
  gr2i->GetXaxis()->SetNdivisions(515);
  gr2i->Draw("AP");

  auto c4=new TCanvas("c4","background histo",1600,900);
  c4->Divide(2,1);
  c4->cd(1);
  RooDataHist bghin("bghin","bghin",RooArgSet(ph),Import(*bghistin)) ;
  RooRealVar m0("m0","median of the lognormal",1200,0,10000) ;
  RooRealVar k("k","shape parameter of the lognormal",2,0,10) ;
  RooLognormal ln("ln","ln",ph,m0,k);
  RooRealVar bgmean("bgmean","mean",1200,0,10000) ;
  RooRealVar bgsigma("bgsigma","sigma",50,0,100) ;
  RooGaussian bggaus("bggaus","bggaus",ph,bgmean,bgsigma);

  //bggaus.fitTo(bghin,Range(1180,1240));
  ln.fitTo(bghin,Range(1140,1240));
  RooPlot* bghframe = ph.frame(Title("bghisto"),Range(1100,1300)) ;
  bghin.plotOn(bghframe);
  //bggaus.plotOn(bghframe);
  ln.plotOn(bghframe);
  RooHist* bghresid;
  bghresid=bghframe->residHist() ;
  bghframe->Draw();
  c4->cd(2);
  RooPlot* bghresframe = ph.frame(Title("bghresframe"),Range(1140,1240)) ;
  bghresframe->addPlotable(bghresid,"P") ;
  bghresframe->Draw();

  auto c5=new TCanvas("c5","inner residual",1600,900);
  c5->Divide(4,2);
  for (int i=0;i<8;i++){
    c5->cd(i+1);
    RooPlot* frame = ph.frame(Title("residual"),Range(1150,1550)) ;
    frame->addPlotable(hresid[i],"P") ;
    frame->Draw();
  }
  auto c6=new TCanvas("c6","outer residual",1600,900);
  c6->Divide(5,4);
  for(Int_t i=8;i<27;i++){
    c6->cd(i-7);
    RooPlot* frame = ph.frame(Title("residual"),Range(1300,1850)) ;
    frame->addPlotable(hresid[i],"P") ;
    frame->Draw();
  }

  auto c7=new TCanvas("c7","total mean gaussians",1600,900);
  c7->Divide(2,1);
  c7->cd(1);
  meanhistin->Draw();
  for(int i=0;i<8;i++){
    meanhist[i]->SetMarkerStyle(kFullCircle);
    meanhist[i]->SetLineColor(kRed);
    meanhist[i]->Draw("same");
  }
  c7->cd(2);
  meanhistout->Draw();
  for(int i=8;i<27;i++){
    meanhist[i]->SetMarkerStyle(kFullCircle);
    meanhist[i]->SetLineColor(kRed);
    meanhist[i]->Draw("same");
  }


  simPdfout.Print("t");
  simPdfin.Print("t");
  //ln.Print("t");
  Double_t integin=meanhistin->Integral(0,1299);
  Double_t integout=meanhistout->Integral(0,1299);

  for(int i=0;i<1300;i++){
    if (meanhistin->Integral(0,i)/integin > 0.16) {
      cout<<"lowin:"<<i<<endl;
      for (int j=i;j<1300;j++) {
        if (meanhistin->Integral(0,j)/integin > 0.84) {
          cout<<"highin:"<<j<<endl;
          break;
        }
      }
      break;
    }
  }

  for(int i=0;i<1300;i++){
    if (meanhistout->Integral(0,i)/integout > 0.16) {
      cout<<"lowout:"<<i<<endl;
      for (int j=i;j<1300;j++) {
        if (meanhistout->Integral(0,j)/integout > 0.84) {
          cout<<"highout:"<<j<<endl;
          break;
        }
      }
      break;
    }
  }
  cout<<"total counts inner:"<<totalcountsin<<endl;
  cout<<"total counts outer:"<<totalcountsout<<endl;
}


