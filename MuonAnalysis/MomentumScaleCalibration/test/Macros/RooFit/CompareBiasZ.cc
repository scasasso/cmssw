#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TSystem.h"
#include "RooMsgService.h"

#include "FitMassSlices.cc"
#include "FitMass1D.cc"
#include "Legend.h"


class CompareBiasZ
{
public:
  CompareBiasZ(const TString& inputFile="0_zmumuHisto.root", const TString& outputFile="BiasCheck_0.root", const TString& sig="relBreitWignerTimesCB", const TString& bkg="exponential", const double& Mmin = 75., const double& Mmax = 105.) 
  {

    gROOT->SetStyle("Plain");    

    TString inputFileName(inputFile);
    TString outputFileName(outputFile);

    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
    FitMassSlices fitter;

    fitter.useChi2 = false;

    // #######################
    // fitter (slicer and fitter)
    // #######################     
    // Further rebin
    fitter.rebinX = 2;
    // fitter.rebinZ = 20;
    // fitter.rebinY = 2;
    fitter.sigma2 = 1.;
    fitter.fit(inputFileName, outputFileName, sig, bkg, 91, Mmin, Mmax, 1.2, 0., 5.);

    // #######################
    // 1D fitter (all events)
    // #######################
    FitMass1D fitMass1D;
    //fitMass1D.rebinX = 10;
    fitMass1D.fitter()->initMean(91, Mmin, Mmax);
    fitMass1D.fitter()->initGamma( 2.4952, 0., 10.);
    fitMass1D.fitter()->gamma()->setConstant(kTRUE);
    fitMass1D.fitter()->initMean2(0., -20., 20.);
    fitMass1D.fitter()->mean2()->setConstant(kTRUE);
    fitMass1D.fitter()->initSigma(1.2, 0., 5.);
    fitMass1D.fitter()->initAlpha(3., 0., 30.);
    fitMass1D.fitter()->initN(2, 0., 50.);
    fitMass1D.fitter()->initExpCoeffA0( -0.3, -10., 10. );
    fitMass1D.fitter()->initFsig(0.9, 0., 1.);
    fitMass1D.fitter()->initA0(0., -10., 10.);
    fitMass1D.fitter()->initA1(0., -10., 10.);
    fitMass1D.fitter()->initA2(0., -10., 10.);
    fitMass1D.fitter()->initA3(0., -10., 10.);
    fitMass1D.fitter()->initA4(0., -10., 10.);
    fitMass1D.fitter()->initA5(0., -10., 10.);
    fitMass1D.fitter()->initA6(0., -10., 10.);
    fitMass1D.fitter()->initExpCoeffA0(-1, -10., 10.);
    fitMass1D.fitter()->initExpCoeffA1(0., -10., 10.);
    fitMass1D.fitter()->initExpCoeffA2(0., -2., 2.);

    /// Let's fit
    fitMass1D.fit(inputFileName, outputFileName, "UPDATE", Mmin, Mmax, sig, bkg);

  }
protected:  

  TFile * file_;
};
