//-----------------------------------------//
//root [0] .L MultiHistoOverlap.C++
//root [1] MultiHistoOverlap("BiasCheck_IDEAL.root=MC design,  BiasCheck_STARTUP.root=MC startup, BiasCheck_MP1073.root=Summer2011 Tk geometry",3)
//root [1] MultiHistoOverlap("BiasCheck_IDEAL.root=MC design,  BiasCheck_STARTUP.root=MC startup, BiasCheck_MP1073.root=Summer2011 Tk geometry, BiasCheck_MP0743.root=Tk Summer2011 (no mass constraint)",4)
//------------------------------------------//

#include <iostream>
#include "Gtypes.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "RooPlot.h"
#include "TAttMarker.h"

#include <iostream>
#include <sstream>
#include <algorithm> 

        
void MultiHistoOverlap(TString namesandlabels, const Int_t nOfFiles, const TString& outDir="./", const TString& legHeader="this validation", const TString& etaString="|#eta_{#mu}|<2.4"){

  Double_t minMplot = 90.;
  Double_t maxMplot = 92.;
 
  gROOT->Reset();

  const Color_t colorlist[7]={kBlack,kRed+1,kBlue+1,kRed+1,kMagenta,kViolet+5,kGreen+2}; //kBlue-4,
  const Int_t linestylelist[7]={1,1,1,1,1,1,1};
  const Int_t stylelist[7]={1,1,1,1,1,1,1};
  const Style_t markerstylelist[7]={kFullCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle};

//   // Thesis plots (MC)
//   const Color_t colorlist[7]={kRed+1,kRed+1,kBlue,kRed,kViolet+5,kOrange,kBlack}; //kBlue-4,
//   const Int_t linestylelist[7]={2,1,1,1,1,1,1};
//   const Int_t stylelist[7]={1,1,1,1,1,1,1};
//   const Style_t markerstylelist[7]={kOpenSquare,kFullSquare,kFullCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle};

//   // Thesis plots (DATA)
//   const Color_t colorlist[7]={kBlack,kBlack,kBlue,kRed,kViolet+5,kOrange,kBlack}; //kBlue-4,
//   const Int_t linestylelist[7]={2,1,1,1,1,1,1};
//   const Int_t stylelist[7]={1,1,1,1,1,1,1};
//   const Style_t markerstylelist[7]={kOpenCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle,kFullCircle};

  // TString text2 = "#sqrt{s} = 8 TeV, L=19.5 fb^{-1}";         
  // TString text2 = "#sqrt{s} = 8 TeV, Run2012D";         
  // TString text2 = "#sqrt{s} = 7 TeV, L=5.1 fb^{-1}";         
  // TString text2 = "7TeV Drell Yan Monte Carlo sample";         
  // TString text2 = "8TeV Drell Yan Monte Carlo sample";         
  TString text2 = "13 TeV Drell Yan Monte Carlo sample";         

 
 
 //  gSystem->Load("libRooFit");
 //  using namespace RooFit;
 // preamble
  TPaveText *cmsprel = new TPaveText(0.09, 0.95, 0.95, 0.99, "NDC");
  cmsprel->SetTextSize(0.035);
  cmsprel->SetTextFont(42);
  cmsprel->SetFillColor(0);
  cmsprel->SetBorderSize(0);
  cmsprel->SetMargin(0.01);
  cmsprel->SetTextAlign(12); // align left

  TString myEtaString;
  TString text = "CMS Preliminary";
  cmsprel->AddText(0.0, 0.5,text);  
  if (etaString=="BAR") {
    myEtaString.Clear();
    myEtaString.Append("|#eta_{#mu}|<0.9");
  } else if (etaString=="BWD") {
    myEtaString.Clear();
    myEtaString.Append("-2.4<#eta_{#mu}<-0.9");
  } else if (etaString=="FWD") {
    myEtaString.Clear();
    myEtaString.Append("0.9<#eta_{#mu}<2.4");
  } else if (etaString=="ALL") {
    myEtaString.Clear();
    myEtaString.Append("|#eta_{#mu}|<2.4");
  }
  else {
    myEtaString.Clear();
    myEtaString.Append(etaString);
  }



  cmsprel->AddText(0.35, 0.5, text2);
  cmsprel->AddText(0.85, 0.5, myEtaString);


  TList* FileList  = new TList();  
  TList* LabelList = new TList();    
  TObjArray *nameandlabelpairs = namesandlabels.Tokenize(",");  
  for (Int_t i = 0; i < nameandlabelpairs->GetEntries(); ++i) {    
    TObjArray *aFileLegPair = TString(nameandlabelpairs->At(i)->GetName()).Tokenize("=");       
    if(aFileLegPair->GetEntries() == 2) {      
      FileList->Add( TFile::Open(aFileLegPair->At(0)->GetName())  ); 
      LabelList->Add( aFileLegPair->At(1) );    
    } else {      
      std::cout << "Please give file name and legend entry in the following form:\n" 		<< " filename1=legendentry1,filename2=legendentry2\n";          
    }  
  }
  

 Int_t NOfFiles =  FileList->GetSize();  
 if ( NOfFiles!=nOfFiles ){
   std::cout<<"&MSG-e: NOfFiles = "<<nOfFiles<<std::endl;  
   return;
 }  
 

 TString LegLabels[nOfFiles];    
 for(Int_t j=0; j < nOfFiles; j++) {       
   TObjString* legend = (TObjString*)LabelList->At(j);    
   LegLabels[j] = legend->String();    
   std::cout<<"LegLabels["<<j<<"]"<<LegLabels[j]<<std::endl;  
 }

 TLegend *leg=0; 

 TCanvas* c0 = new TCanvas("c0", "c0",50, 20, 800,600);
 TCanvas* c1 = new TCanvas("c1", "c1",50, 20, 800,600);
 TCanvas* c2 = new TCanvas("c2", "c2",50, 20, 800,600);
 TCanvas* c3 = new TCanvas("c3", "c3",50, 20, 800,600);
 TCanvas* c4 = new TCanvas("c4", "c4",50, 20, 800,600);
 TCanvas* c5 = new TCanvas("c5", "c5",50, 20, 1200,800);
 TCanvas* c6 = new TCanvas("c6", "c6",50, 20, 1200,800);
 TCanvas* c7 = new TCanvas("c7", "c7",50, 20, 1200,800);

 TCanvas* c0s = new TCanvas("c0s", "c0s",50, 20, 800,600);
 TCanvas* c1s = new TCanvas("c1s", "c1s",50, 20, 800,600);
 TCanvas* c2s = new TCanvas("c2s", "c2s",50, 20, 800,600);
 TCanvas* c3s = new TCanvas("c3s", "c3s",50, 20, 800,600);

 TCanvas* cFit = new TCanvas("cFit", "cFit",50, 20, 1600, 800);

 //----------------- CANVAS C5 --------------//
 c5->SetFillColor(0);  
 c5->cd();
 c5->SetTopMargin(0.12);
 c5->SetRightMargin(0.16);

 leg = new TLegend(0.25,0.75,0.65,0.90);  
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);
 
 // Mass VS muon eta,phi plus -------------------------------
 TH2D *histoMassVsEtaPhiPlus[nOfFiles];

 TStyle *newStyle = new TStyle();
 newStyle->SetPalette(1);
 //newStyle->SetOptStat(0);
 // newStyle->SetOptTitle(1);

 const Int_t NRGBs = 5;
 const Int_t NCont = 255;
 
 Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
 Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
 Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
 Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 newStyle->SetNumberContours(NCont);


 Double_t zMin(86.);
 Double_t zMax(98.);
 Double_t maxDistLow(0.);
 Double_t maxDistUp(0.);
 Double_t maxDist(0.);
 Double_t ZMass = 91.1876;
 //Double_t ZMass = 91.09;

 // Fixed z-axis range
 maxDist=1.3;

//  // Choose the z-axis range looking at the largest deviation among all upwards/downwards deviations
//  for(Int_t j=0; j < nOfFiles; j++){
//    TFile *fin = (TFile*)FileList->At(j);    
//    if ( histoMassVsEtaPhiPlus[j] = (TH2D*)fin->Get("MassVsEtaPhiPlus/allHistos/meanHisto")){
//      zMin = std::max(0,histoMassVsEtaPhiPlus[j]->GetMinimum());
//      zMax = histoMassVsEtaPhiPlus[j]->GetMaximum();
//      std::cout << "### zMin = " << zMin << ", zMax = " << zMax << std::endl;
//      maxDistLow = fabs(ZMass-zMin);
//      maxDistUp = fabs(ZMass-zMax);
//      if (maxDistLow<maxDistUp && maxDistUp>maxDist) maxDist = maxDistUp;
//      else if (maxDistLow>maxDistUp && maxDistLow>maxDist) maxDist = maxDistLow;
//      else continue;
//    }
//  }
 
//  std::cout<< "#### maxDist = " << maxDist << std::endl;

 TPaveText *mapTitlePave = new TPaveText(0.1, 0.88, 0.7, 0.92, "NDC");
  mapTitlePave->SetTextSize(0.045);
  mapTitlePave->SetTextFont(42);
  mapTitlePave->SetFillColor(0);
  mapTitlePave->SetBorderSize(0);
  mapTitlePave->SetMargin(0.01);
  mapTitlePave->SetTextAlign(12); // align left
  gStyle->SetOptTitle(0);

  TString mapTitle = "";

 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsEtaPhiPlus[j] = (TH2D*)fin->Get("MassVsEtaPhiPlus/allHistos/meanHisto")){
     if ( j == 0 ) {
       histoMassVsEtaPhiPlus[j]->SetStats(false);
       histoMassVsEtaPhiPlus[j]->SetTitle(LegLabels[j]);
       histoMassVsEtaPhiPlus[j]->GetXaxis()->SetTitle("positive muon #phi (rad)");
       histoMassVsEtaPhiPlus[j]->GetYaxis()->SetTitle("positive muon #eta");
       histoMassVsEtaPhiPlus[j]->GetZaxis()->SetTitle("M^{fit}_{Z} (GeV)");
//        zMin = histoMassVsEtaPhiPlus[j]->GetMinimum();
//        zMax = histoMassVsEtaPhiPlus[j]->GetMaximum();
//        if (fabs(ZMass-zMin)<fabs(ZMass-zMax)) maxDist = fabs(ZMass-zMax);
//        else maxDist = fabs(ZMass-zMin);
//        zMin = ZMass-maxDist; zMax = ZMass+maxDist;
       histoMassVsEtaPhiPlus[j]->SetMinimum(ZMass-maxDist);
       histoMassVsEtaPhiPlus[j]->SetMaximum(ZMass+maxDist);
       // histoMassVsEtaPhiPlus[j]->GetYaxis()->SetMinimum(-2.4);
       // histoMassVsEtaPhiPlus[j]->GetYaxis()->SetMaximum(2.4);
       histoMassVsEtaPhiPlus[j]->Draw("COLZ");
       mapTitle.Append(LegLabels[j]);
       mapTitlePave->AddText(0.3, 0.5, mapTitle);
       cmsprel->Draw("same");
       mapTitlePave->Draw("same");
       c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file0.png"); c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file0.root"); c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file0.pdf"); c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file0.eps");
       mapTitle.Clear();
       mapTitlePave->Clear();
     } else {
       histoMassVsEtaPhiPlus[j]->SetTitle(LegLabels[j]);
       histoMassVsEtaPhiPlus[j]->SetStats(false);
       histoMassVsEtaPhiPlus[j]->GetXaxis()->SetTitle("positive muon #phi (rad)");
       histoMassVsEtaPhiPlus[j]->GetYaxis()->SetTitle("positive muon #eta");
       histoMassVsEtaPhiPlus[j]->GetZaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsEtaPhiPlus[j]->SetMinimum(ZMass-maxDist);
       histoMassVsEtaPhiPlus[j]->SetMaximum(ZMass+maxDist);

//        histoMassVsEtaPhiPlus[j]->SetMinimum(zMin);
//        histoMassVsEtaPhiPlus[j]->SetMaximum(zMax);
       // histoMassVsEtaPhiPlus[j]->GetYaxis()->SetMinimum(-2.4);
       // histoMassVsEtaPhiPlus[j]->GetYaxis()->SetMaximum(2.4);
       histoMassVsEtaPhiPlus[j]->Draw("COLZ");
       stringstream ssJ (stringstream::in | stringstream::out);
       ssJ << j;
       string j_string = ssJ.str();
       TString j_TString(j_string);
       mapTitle.Append(LegLabels[j]);
       mapTitlePave->AddText(0.3, 0.5, mapTitle);
       cmsprel->Draw("same"); 
       mapTitlePave->Draw("same");
       c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file"+j_TString+".png"); c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file"+j_TString+".root"); c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file"+j_TString+".pdf"); c5->SaveAs(outDir+"/MassVsEtaPhiPlus_file"+j_TString+".eps");
       mapTitle.Clear();
       mapTitlePave->Clear();
     }

   }
 }

 //----------------- CANVAS C6 --------------//
 c6->SetFillColor(0);  
 c6->cd();
 c6->SetTopMargin(0.12);
 c6->SetRightMargin(0.16);

 leg = new TLegend(0.25,0.75,0.65,0.90);  
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);


 // Mass VS muon eta,phi minus -------------------------------
 TH2D *histoMassVsEtaPhiMinus[nOfFiles];

 // // Choose the z-axis range looking at the largest deviation among all upwards/downwards deviations
//  for(Int_t j=0; j < nOfFiles; j++){
//    TFile *fin = (TFile*)FileList->At(j);    
//    if ( histoMassVsEtaPhiMinus[j] = (TH2D*)fin->Get("MassVsEtaPhiMinus/allHistos/meanHisto")){
//      zMin = histoMassVsEtaPhiMinus[j]->GetMinimum();
//      zMax = histoMassVsEtaPhiMinus[j]->GetMaximum();
//      maxDistLow = fabs(ZMass-zMin);
//      maxDistUp = fabs(ZMass-zMax);
//      if (maxDistLow<maxDistUp && maxDistUp>maxDist) maxDist = maxDistUp;
//      else if (maxDistLow>maxDistUp && maxDistLow>maxDist) maxDist = maxDistLow;
//      else continue;
//    }
//  }
 

 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsEtaPhiMinus[j] = (TH2D*)fin->Get("MassVsEtaPhiMinus/allHistos/meanHisto")){
     if ( j == 0 ) {
       histoMassVsEtaPhiMinus[j]->SetStats(false);
       histoMassVsEtaPhiMinus[j]->SetTitle(LegLabels[j]);
       histoMassVsEtaPhiMinus[j]->GetXaxis()->SetTitle("negative muon #phi (rad)");
       histoMassVsEtaPhiMinus[j]->GetYaxis()->SetTitle("negative muon #eta");
       histoMassVsEtaPhiMinus[j]->GetZaxis()->SetTitle("M^{fit}_{Z} (GeV)");
//        zMin = histoMassVsEtaPhiMinus[j]->GetMinimum();
//        zMax = histoMassVsEtaPhiMinus[j]->GetMaximum();
//        if (fabs(ZMass-zMin)<fabs(ZMass-zMax)) maxDist = fabs(ZMass-zMax);
//        else maxDist = fabs(ZMass-zMin);
//        zMin = ZMass-maxDist; zMax = ZMass+maxDist;
       histoMassVsEtaPhiMinus[j]->SetMinimum(ZMass-maxDist);
       histoMassVsEtaPhiMinus[j]->SetMaximum(ZMass+maxDist);
       // histoMassVsEtaPhiMinus[j]->GetYaxis()->SetMinimum(-2.4);
       // histoMassVsEtaPhiMinus[j]->GetYaxis()->SetMaximum(2.4);
       histoMassVsEtaPhiMinus[j]->Draw("COLZ");
       mapTitle.Append(LegLabels[j]);
       mapTitlePave->AddText(0.3, 0.5, mapTitle);
       cmsprel->Draw("same");
       mapTitlePave->Draw("same");
       c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file0.png"); c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file0.root"); c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file0.pdf"); c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file0.eps");
       mapTitle.Clear();
       mapTitlePave->Clear();
     } else {
       histoMassVsEtaPhiMinus[j]->SetTitle(LegLabels[j]);
       histoMassVsEtaPhiMinus[j]->SetStats(false);
       histoMassVsEtaPhiMinus[j]->GetXaxis()->SetTitle("negative muon #phi (rad)");
       histoMassVsEtaPhiMinus[j]->GetYaxis()->SetTitle("negative muon #eta");
       histoMassVsEtaPhiMinus[j]->GetZaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsEtaPhiMinus[j]->SetMinimum(ZMass-maxDist);
       histoMassVsEtaPhiMinus[j]->SetMaximum(ZMass+maxDist);

//        histoMassVsEtaPhiMinus[j]->SetMinimum(zMin);
//        histoMassVsEtaPhiMinus[j]->SetMaximum(zMax);
       // histoMassVsEtaPhiMinus[j]->GetYaxis()->SetMinimum(-2.4);
       // histoMassVsEtaPhiMinus[j]->GetYaxis()->SetMaximum(2.4);
       histoMassVsEtaPhiMinus[j]->Draw("COLZ");
       stringstream ssJ (stringstream::in | stringstream::out);
       ssJ << j;
       string j_string = ssJ.str();
       TString j_TString(j_string);
       mapTitle.Append(LegLabels[j]);
       mapTitlePave->AddText(0.3, 0.5, mapTitle);
       cmsprel->Draw("same"); 
       mapTitlePave->Draw("same");
       c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file"+j_TString+".png"); c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file"+j_TString+".root"); c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file"+j_TString+".pdf"); c6->SaveAs(outDir+"/MassVsEtaPhiMinus_file"+j_TString+".eps");
       mapTitle.Clear();
       mapTitlePave->Clear();
     }

   }
 }
 




 gROOT->ProcessLine(".L tdrstyle.C"); 
 gROOT->ProcessLine("setTDRStyle()");


 //----------------- CANVAS C0 --------------//
 c0->SetFillColor(0);  
 c0->cd();

 leg = new TLegend(0.15,0.15,0.6,0.3);  
 leg->SetHeader(legHeader);
 ////leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);
 
// Mass VS muon phi plus -------------------------------
 TH1D *histoMassVsPhiPlus[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsPhiPlus[j] = (TH1D*)fin->Get("MassVsPhiPlus/allHistos/meanHisto")){
     histoMassVsPhiPlus[j]->SetLineStyle(linestylelist[j]);
     histoMassVsPhiPlus[j]->SetLineWidth(2);
     histoMassVsPhiPlus[j]->SetMarkerColor(colorlist[j]);
     histoMassVsPhiPlus[j]->SetLineColor(colorlist[j]);
     histoMassVsPhiPlus[j]->SetMarkerStyle(markerstylelist[j]); 
//      histoMassVsPhiPlus[j]->GetXaxis()->SetTitleFont(40);
//      histoMassVsPhiPlus[j]->GetYaxis()->SetTitleFont(40);
     //     histoMassVsPhiPlus[j]->SetMarkerSize(0.75);
     if ( j == 0 ) {
       histoMassVsPhiPlus[j]->GetXaxis()->SetTitle("positive muon #phi (rad)");
       histoMassVsPhiPlus[j]->GetXaxis()->SetTitleOffset(1.2);
       histoMassVsPhiPlus[j]->GetYaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsPhiPlus[j]->GetYaxis()->SetTitleOffset(1.2);
       //       histoMassVsPhiPlus[j]->GetYaxis()->SetRangeUser(88.5,93.5);
       histoMassVsPhiPlus[j]->GetYaxis()->SetRangeUser(minMplot, maxMplot);
       histoMassVsPhiPlus[j]->GetXaxis()->SetRangeUser(-3.14,3.14);
       histoMassVsPhiPlus[j]->Draw();
     } else {
       histoMassVsPhiPlus[j]->Draw("SAME");
     }
     leg->AddEntry(histoMassVsPhiPlus[j],LegLabels[j],"PL");  
   }
 }
 cmsprel->Draw("same");
 leg->Draw("same");
 c0->SaveAs(outDir+"/MassVsPhiPlus.png"); c0->SaveAs(outDir+"/MassVsPhiPlus.root"); c0->SaveAs(outDir+"/MassVsPhiPlus.pdf"); c0->SaveAs(outDir+"/MassVsPhiPlus.eps"); 


 //----------------- CANVAS C1 --------------//
 c1->SetFillColor(0);  
 c1->cd();

 leg = new TLegend(0.15,0.15,0.6,0.3);
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);
 
// Mass VS muon eta plus -------------------------------
 TH1D *histoMassVsEtaPlus[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsEtaPlus[j] = (TH1D*)fin->Get("MassVsEtaPlus/allHistos/meanHisto")){
     histoMassVsEtaPlus[j]->SetLineStyle(linestylelist[j]);
     histoMassVsEtaPlus[j]->SetLineWidth(2);
     histoMassVsEtaPlus[j]->SetMarkerColor(colorlist[j]);
     histoMassVsEtaPlus[j]->SetLineColor(colorlist[j]);
     histoMassVsEtaPlus[j]->SetMarkerStyle(markerstylelist[j]); 
     //     histoMassVsEtaPlus[j]->SetMarkerSize(0.75);
     if ( j == 0 ) {
       histoMassVsEtaPlus[j]->GetXaxis()->SetTitle("positive muon #eta");
       histoMassVsEtaPlus[j]->GetXaxis()->SetTitleOffset(1.2);
       histoMassVsEtaPlus[j]->GetYaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsEtaPlus[j]->GetYaxis()->SetTitleOffset(1.2);
       //       histoMassVsEtaPlus[j]->GetYaxis()->SetRangeUser(88.5,93.5);
       histoMassVsEtaPlus[j]->GetYaxis()->SetRangeUser(minMplot, maxMplot);
       histoMassVsEtaPlus[j]->GetXaxis()->SetRangeUser(-2.41,2.41);
       histoMassVsEtaPlus[j]->Draw();
     } else {
       histoMassVsEtaPlus[j]->Draw("SAME");
     }
     leg->AddEntry(histoMassVsEtaPlus[j],LegLabels[j],"PL");  
   }
 }
 cmsprel->Draw("same");
 leg->Draw("same");
 c1->SaveAs(outDir+"/MassVsEtaPlus.png"); c1->SaveAs(outDir+"/MassVsEtaPlus.root"); c1->SaveAs(outDir+"/MassVsEtaPlus.pdf"); c1->SaveAs(outDir+"/MassVsEtaPlus.eps"); 


 //----------------- CANVAS C2 --------------//
 c2->SetFillColor(0);  
 c2->cd();

 leg = new TLegend(0.15,0.15,0.6,0.3);  
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);
 
// Mass VS muon eta plus - eta minus  -------------------------------
 TH1D *histoMassVsEtaPlusMinusDiff[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsEtaPlusMinusDiff[j] = (TH1D*)fin->Get("MassVsEtaPlusMinusDiff/allHistos/meanHisto")){
     histoMassVsEtaPlusMinusDiff[j]->SetLineStyle(linestylelist[j]);
     histoMassVsEtaPlusMinusDiff[j]->SetLineWidth(2);
     histoMassVsEtaPlusMinusDiff[j]->SetMarkerColor(colorlist[j]);
     histoMassVsEtaPlusMinusDiff[j]->SetLineColor(colorlist[j]);
     histoMassVsEtaPlusMinusDiff[j]->SetMarkerStyle(markerstylelist[j]); 
     //     histoMassVsEtaPlusMinusDiff[j]->SetMarkerSize(0.75);
     if ( j == 0 ) {
       histoMassVsEtaPlusMinusDiff[j]->GetXaxis()->SetTitle("#eta pos. muon -  #eta neg. muon");
       histoMassVsEtaPlusMinusDiff[j]->GetXaxis()->SetTitleOffset(1.2);
       histoMassVsEtaPlusMinusDiff[j]->GetYaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsEtaPlusMinusDiff[j]->GetYaxis()->SetTitleOffset(1.2);
       //       histoMassVsEtaPlusMinusDiff[j]->GetYaxis()->SetRangeUser(88.0,96.0);
       histoMassVsEtaPlusMinusDiff[j]->GetYaxis()->SetRangeUser(minMplot, maxMplot);
       histoMassVsEtaPlusMinusDiff[j]->GetXaxis()->SetRangeUser(-3,3);
       histoMassVsEtaPlusMinusDiff[j]->Draw();
     } else {
       histoMassVsEtaPlusMinusDiff[j]->Draw("SAME");
     }
     leg->AddEntry(histoMassVsEtaPlusMinusDiff[j],LegLabels[j],"PL");  
   }
 }
 cmsprel->Draw("same");
 leg->Draw("same");
 c2->SaveAs(outDir+"/MassVsEtaPlusMinusDiff.png"); c2->SaveAs(outDir+"/MassVsEtaPlusMinusDiff.root"); c2->SaveAs(outDir+"/MassVsEtaPlusMinusDiff.pdf"); c2->SaveAs(outDir+"/MassVsEtaPlusMinusDiff.eps");

 //----------------- CANVAS C3 --------------//
 c3->SetFillColor(0);  
 c3->cd();

 leg = new TLegend(0.15,0.15,0.6,0.3);  
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);
 
// Mass VS muon phi minus -------------------------------
 TH1D *histoMassVsPhiMinus[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsPhiMinus[j] = (TH1D*)fin->Get("MassVsPhiMinus/allHistos/meanHisto")){
     histoMassVsPhiMinus[j]->SetLineStyle(linestylelist[j]);
     histoMassVsPhiMinus[j]->SetLineWidth(2);
     histoMassVsPhiMinus[j]->SetMarkerColor(colorlist[j]);
     histoMassVsPhiMinus[j]->SetLineColor(colorlist[j]);
     histoMassVsPhiMinus[j]->SetMarkerStyle(markerstylelist[j]); 
     //     histoMassVsPhiMinus[j]->SetMarkerSize(0.75);
     if ( j == 0 ) {
       histoMassVsPhiMinus[j]->GetXaxis()->SetTitle("negative muon #phi (rad)");
       histoMassVsPhiMinus[j]->GetXaxis()->SetTitleOffset(1.2);
       histoMassVsPhiMinus[j]->GetYaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsPhiMinus[j]->GetYaxis()->SetTitleOffset(1.2);
       //       histoMassVsPhiMinus[j]->GetYaxis()->SetRangeUser(88.5,93.5);
       histoMassVsPhiMinus[j]->GetYaxis()->SetRangeUser(minMplot, maxMplot);
       histoMassVsPhiMinus[j]->GetXaxis()->SetRangeUser(-3.14,3.14);
       histoMassVsPhiMinus[j]->Draw();
     } else {
       histoMassVsPhiMinus[j]->Draw("SAME");
     }
     leg->AddEntry(histoMassVsPhiMinus[j],LegLabels[j],"PL");  
   }
 }
 cmsprel->Draw("same");
 leg->Draw("same");
 c3->SaveAs(outDir+"/MassVsPhiMinus.png"); c3->SaveAs(outDir+"/MassVsPhiMinus.root"); c3->SaveAs(outDir+"/MassVsPhiMinus.pdf"); c3->SaveAs(outDir+"/MassVsPhiMinus.eps"); 


 //----------------- CANVAS C4 --------------//
 c4->SetFillColor(0);  
 c4->cd();

 leg = new TLegend(0.15,0.15,0.6,0.3);  
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);
 
// Mass VS muon eta minus -------------------------------
 TH1D *histoMassVsEtaMinus[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoMassVsEtaMinus[j] = (TH1D*)fin->Get("MassVsEtaMinus/allHistos/meanHisto")){
     histoMassVsEtaMinus[j]->SetLineStyle(linestylelist[j]);
     histoMassVsEtaMinus[j]->SetLineWidth(2);
     histoMassVsEtaMinus[j]->SetMarkerColor(colorlist[j]);
     histoMassVsEtaMinus[j]->SetLineColor(colorlist[j]);
     histoMassVsEtaMinus[j]->SetMarkerStyle(markerstylelist[j]); 
     //     histoMassVsEtaMinus[j]->SetMarkerSize(0.75);
     if ( j == 0 ) {
       histoMassVsEtaMinus[j]->GetXaxis()->SetTitle("negative muon #eta");
       histoMassVsEtaMinus[j]->GetXaxis()->SetTitleOffset(1.2);
       histoMassVsEtaMinus[j]->GetYaxis()->SetTitle("M^{fit}_{Z} (GeV)");
       histoMassVsEtaMinus[j]->GetYaxis()->SetTitleOffset(1.2);
       //       histoMassVsEtaMinus[j]->GetYaxis()->SetRangeUser(88.5,93.5);
       histoMassVsEtaMinus[j]->GetYaxis()->SetRangeUser(minMplot, maxMplot);
       histoMassVsEtaMinus[j]->GetXaxis()->SetRangeUser(-2.41,2.41);
       histoMassVsEtaMinus[j]->Draw();
     } else {
       histoMassVsEtaMinus[j]->Draw("SAME");
     }
     leg->AddEntry(histoMassVsEtaMinus[j],LegLabels[j],"PL");  
   }
 }
 cmsprel->Draw("same");
 leg->Draw("same");
 c4->SaveAs(outDir+"/MassVsEtaMinus.png"); c4->SaveAs(outDir+"/MassVsEtaMinus.root"); c4->SaveAs(outDir+"/MassVsEtaMinus.pdf"); c4->SaveAs(outDir+"/MassVsEtaMinus.eps"); 


//  //----------------- CANVAS C7 --------------//
//  c7->SetFillColor(0);  
//  c7->cd();

//  leg = new TLegend(0.20,0.2,0.6,0.35);  
//  leg->SetHeader(legHeader);
//  //leg->SetBorderSize(1);
//  leg->SetFillColor(0);
//  leg->SetTextFont(42); leg->SetTextSize(0.04);
 
// // Mass VS muon eta minus -------------------------------
//  TH1D *histoMassVsPt[nOfFiles];
//  for(Int_t j=0; j < nOfFiles; j++) {     
   
//    TFile *fin = (TFile*)FileList->At(j);    
//    if ( histoMassVsPt[j] = (TH1D*)fin->Get("MassVsPt/allHistos/meanHisto")){
//      histoMassVsPt[j]->SetLineStyle(linestylelist[j]);
//      histoMassVsPt[j]->SetLineWidth(2);
//      histoMassVsPt[j]->SetMarkerColor(colorlist[j]);
//      histoMassVsPt[j]->SetLineColor(colorlist[j]);
//      histoMassVsPt[j]->SetMarkerStyle(markerstylelist[j]); 
//      //     histoMassVsPt[j]->SetMarkerSize(0.75);
//      if ( j == 0 ) {
//        histoMassVsPt[j]->GetXaxis()->SetTitle("muon p_{T}");
//        histoMassVsPt[j]->GetYaxis()->SetTitle("M^{fit}_{Z} (GeV)");
//        //       histoMassVsPt[j]->GetYaxis()->SetRangeUser(88.5,93.5);
//        histoMassVsPt[j]->GetYaxis()->SetRangeUser(minMplot, maxMplot);
//        histoMassVsPt[j]->GetXaxis()->SetRangeUser(20.,200.);
//        histoMassVsPt[j]->Draw();
//      } else {
//        histoMassVsPt[j]->Draw("SAME");
//      }
//      leg->AddEntry(histoMassVsPt[j],LegLabels[j],"PL");  
//    }
//  }
//  //cmsprel->Draw("same");
//  leg->Draw("same");
//  c7->SaveAs(outDir+"/MassVsPt.png"); 



//  //----------------- CANVAS C0S --------------//
//  c0s->SetFillColor(0);  
//  c0s->cd();

//  leg = new TLegend(0.15,0.15,0.6,0.3);  
//   leg->SetHeader(legHeader);
//  //leg->SetBorderSize(1);
//  leg->SetFillColor(0);
//  leg->SetTextFont(42); leg->SetTextSize(0.04);

// // Sigma VS muon phi plus -------------------------------
//  TH1D *histoSigmaVsPhiPlus[nOfFiles];
//  for(Int_t j=0; j < nOfFiles; j++) {     
   
//    TFile *fin = (TFile*)FileList->At(j);    
//    if ( histoSigmaVsPhiPlus[j] = (TH1D*)fin->Get("MassVsPhiPlus/allHistos/sigmaHisto")){
//      histoSigmaVsPhiPlus[j]->SetLineStyle(linestylelist[j]);
//      histoSigmaVsPhiPlus[j]->SetMarkerColor(colorlist[j]);
//      histoSigmaVsPhiPlus[j]->SetLineColor(colorlist[j]);
//      histoSigmaVsPhiPlus[j]->SetMarkerStyle(markerstylelist[j]); 
//      //     histoSigmaVsPhiPlus[j]->SetMarkerSize(0.75);
//      if ( j == 0 ) {
//        histoSigmaVsPhiPlus[j]->GetXaxis()->SetTitle("positive muon #phi (rad)");
//        histoSigmaVsPhiPlus[j]->GetYaxis()->SetTitle("#sigma(M^{fit}_{Z}) (GeV)");
//        //       histoSigmaVsPhiPlus[j]->GetYaxis()->SetRangeUser(88.5,93.5);
//        histoSigmaVsPhiPlus[j]->GetYaxis()->SetRangeUser(0.,3.);
//        histoSigmaVsPhiPlus[j]->GetXaxis()->SetRangeUser(-3.14,3.14);
//        histoSigmaVsPhiPlus[j]->Draw();
//      } else {
//        histoSigmaVsPhiPlus[j]->Draw("SAME");
//      }
//      leg->AddEntry(histoSigmaVsPhiPlus[j],LegLabels[j],"PL");  
//    }
//  }
//  //cmsprel->Draw("same");
//  leg->Draw("same");
//  c0s->SaveAs(outDir+"/SigmaVsPhiPlus.png"); 


 //----------------- CANVAS C1S --------------//
 c1s->SetFillColor(0);  
 c1s->cd();

 leg = new TLegend(0.25,0.7,0.7,0.85);  
 leg->SetHeader(legHeader);
 //leg->SetBorderSize(1);
 leg->SetFillColor(0);
 leg->SetTextFont(42); leg->SetTextSize(0.04);

 
// Sigma VS muon eta plus -------------------------------
 TH1D *histoSigmaVsEtaPlus[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoSigmaVsEtaPlus[j] = (TH1D*)fin->Get("MassVsEtaPlus/allHistos/sigmaHisto")){
     histoSigmaVsEtaPlus[j]->SetLineStyle(linestylelist[j]);
     histoSigmaVsEtaPlus[j]->SetLineWidth(2);
     histoSigmaVsEtaPlus[j]->SetMarkerColor(colorlist[j]);
     histoSigmaVsEtaPlus[j]->SetLineColor(colorlist[j]);
     histoSigmaVsEtaPlus[j]->SetMarkerStyle(markerstylelist[j]); 
     //     histoSigmaVsEtaPlus[j]->SetMarkerSize(0.75);
     if ( j == 0 ) {
       histoSigmaVsEtaPlus[j]->GetXaxis()->SetTitle("positive muon #eta");
       histoSigmaVsEtaPlus[j]->GetXaxis()->SetTitleOffset(1.2);
       histoSigmaVsEtaPlus[j]->GetYaxis()->SetTitle("#sigma(M_{#mu#mu}) (GeV)");
       histoSigmaVsEtaPlus[j]->GetXaxis()->SetTitleOffset(1.2);
       //       histoSigmaVsEtaPlus[j]->GetYaxis()->SetRangeUser(88.5,93.5);
       histoSigmaVsEtaPlus[j]->GetYaxis()->SetRangeUser(0.7,3.5);
       histoSigmaVsEtaPlus[j]->GetXaxis()->SetRangeUser(-2.41,2.41);
       histoSigmaVsEtaPlus[j]->Draw();
     } else {
       histoSigmaVsEtaPlus[j]->Draw("SAME");
     }
     leg->AddEntry(histoSigmaVsEtaPlus[j],LegLabels[j],"PL");  
   }
 }
 cmsprel->Draw("same");
 leg->Draw("same");
 c1s->SaveAs(outDir+"/SigmaVsEtaPlus.png"); c1s->SaveAs(outDir+"/SigmaVsEtaPlus.root"); c1s->SaveAs(outDir+"/SigmaVsEtaPlus.pdf"); c1s->SaveAs(outDir+"/SigmaVsEtaPlus.eps");


//  //----------------- CANVAS C2S --------------//
//  c2s->SetFillColor(0);  
//  c2s->cd();

//   leg = new TLegend(0.15,0.15,0.6,0.3);  
//  leg->SetHeader(legHeader);
//  //leg->SetBorderSize(1);
//  leg->SetFillColor(0);
//  leg->SetTextFont(42); leg->SetTextSize(0.04);

// // Sigma VS muon eta plus - eta minus  -------------------------------
//  TH1D *histoSigmaVsEtaPlusMinusDiff[nOfFiles];
//  for(Int_t j=0; j < nOfFiles; j++) {     
   
//    TFile *fin = (TFile*)FileList->At(j);    
//    if ( histoSigmaVsEtaPlusMinusDiff[j] = (TH1D*)fin->Get("MassVsEtaPlusMinusDiff/allHistos/sigmaHisto")){
//      histoSigmaVsEtaPlusMinusDiff[j]->SetLineStyle(linestylelist[j]);
//      histoSigmaVsEtaPlusMinusDiff[j]->SetMarkerColor(colorlist[j]);
//      histoSigmaVsEtaPlusMinusDiff[j]->SetLineColor(colorlist[j]);
//      histoSigmaVsEtaPlusMinusDiff[j]->SetMarkerStyle(markerstylelist[j]); 
//      //     histoSigmaVsEtaPlusMinusDiff[j]->SetMarkerSize(0.75);
//      if ( j == 0 ) {
//        histoSigmaVsEtaPlusMinusDiff[j]->GetXaxis()->SetTitle("#eta pos. muon - #eta neg. muon");
//        histoSigmaVsEtaPlusMinusDiff[j]->GetYaxis()->SetTitle("#sigma(M^{fit}_{Z}) (GeV)");
//        //       histoSigmaVsEtaPlusMinusDiff[j]->GetYaxis()->SetRangeUser(88.0,96.0);
//        histoSigmaVsEtaPlusMinusDiff[j]->GetYaxis()->SetRangeUser(0.,3.);
//        //histoSigmaVsEtaPlusMinusDiff[j]->GetYaxis()->SetRangeUser(90.60,90.75);
//        histoSigmaVsEtaPlusMinusDiff[j]->GetXaxis()->SetRangeUser(-3.2,3.2);
//        histoSigmaVsEtaPlusMinusDiff[j]->Draw();
//      } else {
//        histoSigmaVsEtaPlusMinusDiff[j]->Draw("SAME");
//      }
//      leg->AddEntry(histoSigmaVsEtaPlusMinusDiff[j],LegLabels[j],"PL");  
//    }
//  }
//  //cmsprel->Draw("same");
//  leg->Draw("same");
//  c2s->SaveAs(outDir+"/SigmaVsEtaPlusMinusDiff.png"); 


//  //----------------- CANVAS C3S --------------//
//   leg = new TLegend(0.15,0.15,0.6,0.3);  
//  leg->SetHeader(legHeader);
//  //leg->SetBorderSize(1);
//  leg->SetFillColor(0);
//  leg->SetTextFont(42); leg->SetTextSize(0.04);

//  c3s->SetFillColor(0);  
//  c3s->cd();

 
// // Sigma VS muon pT  -------------------------------
//  TH1D *histoSigmaVsPt[nOfFiles];
//  for(Int_t j=0; j < nOfFiles; j++) {     
   
//    TFile *fin = (TFile*)FileList->At(j);    
//    if ( histoSigmaVsPt[j] = (TH1D*)fin->Get("MassVsPt/allHistos/sigmaHisto")){
//      histoSigmaVsPt[j]->SetLineStyle(linestylelist[j]);
//      histoSigmaVsPt[j]->SetMarkerColor(colorlist[j]);
//      histoSigmaVsPt[j]->SetLineColor(colorlist[j]);
//      histoSigmaVsPt[j]->SetMarkerStyle(markerstylelist[j]); 
//      //     histoSigmaVsPt[j]->SetMarkerSize(0.75);
//      if ( j == 0 ) {
//        histoSigmaVsPt[j]->GetXaxis()->SetTitle("muon p_{T} (GeV)");
//        histoSigmaVsPt[j]->GetYaxis()->SetTitle("#sigma(M^{fit}_{Z}) (GeV)");
//        //       histoSigmaVsPt[j]->GetYaxis()->SetRangeUser(88.0,96.0);
//        histoSigmaVsPt[j]->GetYaxis()->SetRangeUser(0.,3.);
//        //histoSigmaVsPt[j]->GetYaxis()->SetRangeUser(90.60,90.75);
//        histoSigmaVsPt[j]->GetXaxis()->SetRangeUser(15.,105.);
//        histoSigmaVsPt[j]->Draw();
//      } else {
//        histoSigmaVsPt[j]->Draw("SAME");
//      }
//      leg->AddEntry(histoSigmaVsPt[j],LegLabels[j],"PL");  
//    }
//  }
//  //cmsprel->Draw("same");
//  leg->Draw("same");
//  c3s->SaveAs(outDir+"/SigmaVsPt.png"); 

 //----------------- CANVAS CFIT --------------//
 cFit->SetFillColor(0);  
 cFit->cd();
 Float_t nN = TMath::Sqrt(nOfFiles);
 Int_t nX = (Int_t)nN;
 if ( nN-nX > 0.5 ) nX++;
 Int_t nY = (Int_t)(nOfFiles/nX);
 std::cout << nX << " ," << nY << std::endl;
 cFit->Divide(nOfFiles,1);
 
// Mass VS muon phi plus -------------------------------
 TFile *ZFitFile = new TFile("ZFitFile.root","RECREATE");
 RooPlot *histoLineShape[nOfFiles];
 for(Int_t j=0; j < nOfFiles; j++) {     
   
   TFile *fin = (TFile*)FileList->At(j);    
   if ( histoLineShape[j] = (RooPlot*)fin->Get("hRecBestResAllEvents_Mass_frame")){
     std::cout<<"Writing fit histogrem file n. "<<j<<std::endl;
     histoLineShape[j]->Write();
     cFit->cd(j+1);
     histoLineShape[j]->SetTitle(LegLabels[j]);
     histoLineShape[j]->Draw();
     histoLineShape[j]->GetXaxis()->SetTitle("M^{fit}_{Z} (GeV)");
//      TPaveText *cmsprel2 = new TPaveText(0.19, 0.95, 0.95, 0.99, "NDC");
//      cmsprel2->SetTextSize(0.03);
//      cmsprel2->SetTextFont(42);
//      cmsprel2->SetFillColor(0);
//      cmsprel2->SetBorderSize(0);
//      cmsprel2->SetMargin(0.01);
//      cmsprel2->SetTextAlign(12); // align left
//      cmsprel2->AddText(0.666666, 0.5, LegLabels[j]);

   }
 }
 ZFitFile->Close();
 // cmsprel2->Draw("same");
 cFit->SaveAs(outDir+"/ZFitFile.root");
 //cFit->SaveAs(outDir+"/ZFitFile.png");


 
 return; 
};

