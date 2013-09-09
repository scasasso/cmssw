#ifndef MuonPair_h
#define MuonPair_h
 
#include <TObject.h>
#include "MuonAnalysis/MomentumScaleCalibration/interface/Muon.h"

/**
 * Simple class used to save the muon pairs in a root tree. <br>
 * Includes the information on the run and event number.
 */
 
class MuonPair : public TObject
{
 public:
  MuonPair() :
    //initialize 2 object of class muon
    mu1(lorentzVector(0,0,0,0),-1),
    mu2(lorentzVector(0,0,0,0),1),
    run(0),
    event(0)
      {}
    
    MuonPair(const MuScleFitMuon & initMu1, const MuScleFitMuon & initMu2, const UInt_t & initRun=0,  const UInt_t & initEvent=0) :
      mu1(initMu1),
      mu2(initMu2),
      run(initRun),
      event(initEvent)
	{}

      
      /// Used to copy the content of another MuonPair
      void copy(const MuonPair & copyPair)
      {
	mu1 = copyPair.mu1;
	mu2 = copyPair.mu2;
	run= copyPair.run;
	event= copyPair.event;
      }
      
      MuScleFitMuon mu1;
      MuScleFitMuon mu2;
      UInt_t run;
      UInt_t event;
      
      ClassDef(MuonPair, 3)
	};
ClassImp(MuonPair)
   
#endif
