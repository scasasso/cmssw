import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# usage w/ command line options
# cmsRun Zmumu_validator.py  sampleName=TEST 

### helper class definition ###

####################
class SetParameters:
####################
    def __init__(self):
        self.parFix = cms.vint32()
        self.parOrder = cms.vint32()
        self.par = cms.vdouble()
        self.parStep = cms.untracked.vdouble()
        self.parMin  = cms.untracked.vdouble()
        self.parMax  = cms.untracked.vdouble()

    def set(self, fix, order, value, step, min, max):
        self.parFix.append(fix)
        self.parOrder.append(order)
        self.par.append(value)
        self.parStep.append(step)
        self.parMin.append(min)
        self.parMax.append(max)
### end of the  definition ###
        
### fit scale

scaleSetter = SetParameters()
#fix, order, value, step, min, max

##V6 (change initialization of 2nd order harmonics)
scaleSetter.set( 1 ,0, 0.           , 0.0001,   -0.01  , 0.01 )    # (1+par[0])

# VBWD
scaleSetter.set( 0 ,0, 0.0005       , 0.00015,   0.    , 0.0015  )  # par[1] BWD  par[1]*sin(phi+par[2])+par[3]*eta
scaleSetter.set( 0 ,0, 0.           , 0.4,      -3.1416, 3.1416  )  # par[2]
scaleSetter.set( 0 ,0, 0.           , 0.0003,   -0.0015, 0.0015 )  # par[3]
#
scaleSetter.set( 1 ,0,-2.1          , 0.1,      -2.21   ,-1.99 )      # par[4] VBWD/BWD boundary
# BWD
scaleSetter.set( 0 ,0, 0.0005       , 0.00015,   0.    , 0.0015 )  # par[5] BAR  par[5]*sin(phi+par[6])+par[7]*eta
scaleSetter.set( 0 ,0, 0.           , 0.4,      -3.1416, 3.1416 )  # par[6]
scaleSetter.set( 0 ,0, 0.           , 0.0002,   -0.001 , 0.001  )  # par[7]
#
scaleSetter.set( 1 ,0,-1.5          , 0.1,      -1.61  ,-1.39 )      # par[8] BWD/BAR boundary (NB min!=max otherwise MINUIT overwrites the settings)
# BAR
scaleSetter.set( 0 ,0, 0.0005       , 0.00015,   0.    , 0.0015 )  # par[9] FWD  par[9]*sin(phi+par[10])+par[11]*eta
scaleSetter.set( 0 ,0, 0            , 0.4,      -3.1416, 3.1416 )  # par[10]
scaleSetter.set( 0 ,0, 0            , 0.0002,   -0.001,  0.001  )  # par[11]
#
scaleSetter.set( 1 ,0, 1.5          , 0.1,       1.39    , 1.61 )     # par[12] BAR/FWD boundary (NB min!=max otherwise MINUIT overwrites the settings)
# FWD
scaleSetter.set( 0 ,0, 0.0005        , 0.00015,  0.    , 0.0015 )  # par[13] FWD  par[13]*sin(phi+par[14])+par[15]*eta
scaleSetter.set( 0 ,0, 0             , 0.4,     -3.1416, 3.1416 )  # par[14]
scaleSetter.set( 0 ,0, 0.            , 0.0002,  -0.001 , 0.001  )  # par[15]
#
scaleSetter.set( 1, 0, 2.1           , 0.1,       1.99  ,  2.21 )     # par[16] FWD/VFWD boundary
# VFWD
scaleSetter.set( 0 ,0, 0.0005        , 0.00015,   0.    , 0.0015 )  # par[17] FWD  par[17]*sin(phi+par[18])+par[19]*eta
scaleSetter.set( 0 ,0, 0.            , 0.4,      -3.1416, 3.1416 )  # par[18]
scaleSetter.set( 0 ,0, 0.            , 0.0003,   -0.0015 ,0.0015 )  # par[19]

scaleSetter.set( 0 ,0, 5e-05         , 0.00005,  -0.00015 , 0.00015 )  # par[20] charge dep. bias

#VBWD
scaleSetter.set( 0 ,0, 0.0005        , 0.00015,   0.      , 0.0015  )  #par[21]*sin(par[22]*phi+par[23])
scaleSetter.set( 1 ,0, 2.            , 0.1,       1.      , 5.      )  
scaleSetter.set( 0 ,0, 0.            , 0.4,      -3.1416  , 3.1416  )  

#VFWD
scaleSetter.set( 0 ,0, 0.0005        , 0.00015,   0.      , 0.0015  )  #par[24]*sin(par[25]*phi+par[26])
scaleSetter.set( 1 ,0, 2.            , 0.1,       1.      , 5.      )  
scaleSetter.set( 0 ,0, 0.            , 0.4,      -3.1416  , 3.1416  )  



### fit resolution
resolSetter = SetParameters()
#fix, order, value, step, min, max
resolSetter.set( 0 ,0, 0.00025    , 0.0001,   0.   , 0.10  ) # par[0]
resolSetter.set( 0 ,0, 0.030    , 0.005,  -0.05 , 0.05  ) # par[1] 
resolSetter.set( 0 ,0, 0.015    , 0.005,  -0.05 , 0.05  ) # par[2]
resolSetter.set( 0 ,0, 0.010    , 0.005,  -0.05 , 0.05  ) # par[3] 
resolSetter.set( 0 ,0, 0.010    , 0.005,  -0.05 , 0.05  ) # par[4]
resolSetter.set( 0 ,0, 0.010    , 0.005,  -0.05 , 0.05  ) # par[5] 
resolSetter.set( 0 ,0, 0.005    , 0.005,  -0.05 , 0.05  ) # par[6] 
resolSetter.set( 0 ,0, 0.005    , 0.005,  -0.05 , 0.05  ) # par[7] 
resolSetter.set( 0 ,0, 0.010    , 0.005,  -0.05 , 0.05  ) # par[8]
resolSetter.set( 0 ,0, 0.010    , 0.005,  -0.05 , 0.05  ) # par[9] 
resolSetter.set( 0 ,0, 0.010    , 0.005,  -0.05 , 0.05  ) # par[10]
resolSetter.set( 0 ,0, 0.015    , 0.005,  -0.05 , 0.05  ) # par[11]
resolSetter.set( 0 ,0, 0.030    , 0.005,  -0.05 , 0.05  ) # par[12]   

process = cms.Process("TEST")

# Messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


#setup 'standard'  options
options = VarParsing.VarParsing()
# setup any defaults you want

options.register('sampleName',
                 "TEST", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "sample name")

options.register('appendName',
                 "ALL", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string,         # string, int, or float
                 "additional name for the output")

## eta ranges steerable
options.register('etaMax1',
                 2.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "eta max (muon1)")

options.register('etaMin1',
                 -2.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "eta min (muon1)")

options.register('etaMax2',
                 2.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "eta max (muon2)")

options.register('etaMin2',
                 -2.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "eta min (muon2)")



options.parseArguments()


#### INPUT Tree: to be changed in case the input is the Tree ############################


#### SETTING THE VARIABLES #######################


EOSBASEDIR="root://eoscms//eos/cms/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/MC_13TeV" 
TREEINPUTNAME=EOSBASEDIR+"/zmumuTree_DYToMuMu_13TeV_Phys14DR-PU20bx25.root" #PHYS14

TREEOUTPUTNAME=""

print "Running on sample _%s_ , selecting muons in the range: eta_mu1=[%f,%f], eta_mu2=[%f,%f]" % (options.sampleName,options.etaMin1,options.etaMax1,options.etaMin2,options.etaMax2)

# Use this when running on a tree
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(0))

    
process.looper = cms.Looper(
    "MuScleFit",
    # Only used when reading events from a root tree
    MaxEventsFromRootTree = cms.int32(10),

    # Specify a file if you want to read events from a root tree in a local file.
    # In this case the input source should be an empty source with 0 events.
    
    InputRootTreeFileName = cms.string(TREEINPUTNAME),
    
    # Specify the file name where you want to save a root tree with the muon pairs.
    # Leave empty if no file should be written.
    
    OutputRootTreeFileName = cms.string(TREEOUTPUTNAME),
    

    # Choose the kind of muons you want to run on
    # -------------------------------------------
    MuonLabel = cms.InputTag("TrackRefitter"),
    MuonType = cms.int32(2),

    # This line allows to switch to PAT muons. Default is false.
    # Note that the onia selection works only with onia patTuples.
    PATmuons = cms.untracked.bool(False),

    # ---------------- #
    # Select resonance #
    # ---------------- #
    # The resonances are to be specified in this order:
    # Z0, Y(3S), Y(2S), Y(1S), Psi(2S), J/Psi
    # -------------------------------------------------
    resfind = cms.vint32(1, 0, 0, 0, 0, 0),

    # Likelihood settings
    # -------------------
    # step #1
    maxLoopNumber = cms.untracked.int32(4),
    # Select which fits to do in which loop (0 = do not, 1 = do)
    doResolFit =        cms.vint32(1,0,0,0),
    doScaleFit =        cms.vint32(0,1,1,0),
    doBackgroundFit =   cms.vint32(0,0,0,0),
    doCrossSectionFit = cms.vint32(0,0,0,0),

    # Use the probability file or not. If not it will perform a simpler selection taking the muon pair with
    # invariant mass closer to the pdf value and will crash if some fit is attempted.
    UseProbsFile = cms.untracked.bool(True),

    # False = use also MC information
    speedup = cms.bool(True),
    # Set this to false if you do not want to use simTracks.
    # (Note that this is skipped anyway if speedup == True).
    compareToSimTracks = cms.bool(False),

    # Output settings
    # ---------------
    OutputFileName = cms.untracked.string("zmumuHisto_"+options.sampleName+"_"+options.appendName+".root"),

    # BiasType=0 means no bias to muon momenta
    # ----------------------------------------
    BiasType = cms.int32(0),
    parBias = cms.vdouble(),

    # SmearType=0 means no smearing applied to muon momenta
    # -----------------------------------------------------
    SmearType = cms.int32(0),
    parSmear = cms.vdouble(),

    ### taken from AN 2011-166 #########################

    # ------------------------- #
    # Resolution fit parameters #
    # ------------------------- #

    ResolFitType = cms.int32(45),
    parResolFix   = resolSetter.parFix,
    parResolOrder = resolSetter.parOrder,
    parResol      = resolSetter.par,
    parResolStep  = resolSetter.parStep,
    parResolMin   = resolSetter.parMin,
    parResolMax   = resolSetter.parMax,


    # -------------------- #
    # Scale fit parameters #
    # -------------------- #

    # -----------------------------------------------------------------------------------
    ScaleFitType  = cms.int32(50),
    parScaleOrder = scaleSetter.parOrder,
    parScaleFix   = scaleSetter.parFix,
    parScale      = scaleSetter.par,
    parScaleStep  = scaleSetter.parStep,
    parScaleMin   = scaleSetter.parMin,
    parScaleMax   = scaleSetter.parMax,
    

    # ---------------------------- #
    # Cross section fit parameters #
    # ---------------------------- #
    # Note that the cross section fit works differently than the others, it
    # fits ratios of parameters. Fix and Order should not be used as is, they
    # are there mainly for compatibility.
    parCrossSectionOrder = cms.vint32(0, 0, 0, 0, 0, 0),
    parCrossSectionFix   = cms.vint32(0, 0, 0, 0, 0, 0),
    parCrossSection      = cms.vdouble(1.233, 2.07, 6.33, 13.9, 2.169, 127.2),

    # ------------------------- #
    # Background fit parameters #
    # ------------------------- #
    # Window factors for: Z, Upsilons and (J/Psi,Psi2S) regions
    LeftWindowBorder  = cms.vdouble(75., 8.7, 1.391495),
    RightWindowBorder = cms.vdouble(100., 9.7, 5.391495),

    # The two parameters of BgrFitType=2 are respectively:
    # bgr fraction, (negative of) bgr exp. slope, bgr constant
    # --------------------------------------------------------
    # The function types for resonances in a region must be the same
    
    BgrFitType = cms.vint32(2, 2, 2), 
    # These empty parameters should be used when there is no background
    parBgr = cms.vdouble(0., 0.,   0., 0.,   0., 0.,
                         0., 0.,   0., 0.,   0., 0.,   0., 0.,   0., 0.,   0., 0.),
    parBgrFix = cms.vint32(0, 0,   0, 0,   0, 0,
                           # The rest of the parameters is used for the resonance regions. They are automatically fixed in the code
                           # because they are never used to fit the background, but only after the rescaling.
                           1, 1,   1, 1,   1, 1,   1, 1,   1, 1,   1, 1),
    parBgrOrder = cms.vint32(0, 0,   0, 0,   0, 0,
                             0, 0,   0, 0,   0, 0,   0, 0,   0, 0,   0, 0),


    # ----------------------- #

    # Set Minuit fit strategy
    FitStrategy = cms.int32(2),

    # Fit accuracy and debug parameters
    StartWithSimplex   = cms.bool(False),
    ComputeMinosErrors = cms.bool(False),
    MinimumShapePlots  = cms.bool(False),


    ########## TO BE ENABLED ################################
    # Set the cuts on muons to be used in the fit
    SeparateRanges = cms.untracked.bool(True),
    MinMuonEtaFirstRange = cms.untracked.double(options.etaMin1),
    MaxMuonEtaFirstRange = cms.untracked.double(options.etaMax1),
    MinMuonEtaSecondRange = cms.untracked.double(options.etaMin2),
    MaxMuonEtaSecondRange = cms.untracked.double(options.etaMax2),
    PileUpSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
    PrimaryVertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
    
    #ProbabilitiesFileInPath = cms.untracked.string("MuonAnalysis/MomentumScaleCalibration/test/Probs_merge.root"),

    RapidityBinsForZ = cms.untracked.bool(False),

    ProbabilitiesFile = cms.untracked.string("Probs_Z_1001_8TeV_v2.root"), # 8 Tev

    # The following parameters can be used to filter events
    TriggerResultsLabel = cms.untracked.string("TriggerResults"),
    TriggerResultsProcess = cms.untracked.string("HLT"),
    TriggerPath = cms.untracked.vstring(""),
    # Negate the result of the trigger
    NegateTrigger = cms.untracked.bool(False),
    debug = cms.untracked.int32(0),
)
