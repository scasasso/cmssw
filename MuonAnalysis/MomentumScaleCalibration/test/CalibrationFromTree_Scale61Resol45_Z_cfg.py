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

# ### V0 (starting from V1 of Scale56 and adding initialization of the additionale 3 parameters) -> GOOD for data 2011 Legacy
# scaleSetter.set(1 ,0, 0. , 0.0001,   -0.01  , 0.01 )    # (1+par[0])
# scaleSetter.set(0, 0, -0.05, 0.01, -0.11, 0.01) #par[1]
# scaleSetter.set(0, 0, -0.07, 0.01, -0.15, 0.01) #par[2]
# scaleSetter.set(0, 0, -0.10, 0.02, -0.22, 0.02) #par[3]
# scaleSetter.set(0, 0, -0.13, 0.025, -0.25, 0.02) #par[4]
# scaleSetter.set(1, 0,   0.9,  0.1,   0.5,  1.1) #par[5]
# scaleSetter.set(1, 0,   1.6,  0.1,   1.3,  1.8) #par[6]
# scaleSetter.set(1, 0,   2.1,  0.1,   1.9,  2.3) #par[7]


### V1 (same as V0 but range of parameters is revisited for the MC)
scaleSetter.set(1 ,0, 0. , 0.0001,   -0.01  , 0.01 )    # (1+par[0])
scaleSetter.set(0, 0, 0., 0.02, -0.1, 0.1) #par[1]
scaleSetter.set(0, 0, 0., 0.03, -0.15, 0.15) #par[2]
scaleSetter.set(0, 0, 0., 0.05, -0.25, 0.25) #par[3]
scaleSetter.set(0, 0, 0., 0.05, -0.25, 0.25) #par[4]
scaleSetter.set(1, 0,   0.9,  0.1,   0.5,  1.1) #par[5]
scaleSetter.set(1, 0,   1.6,  0.1,   1.3,  1.8) #par[6]
scaleSetter.set(1, 0,   2.1,  0.1,   1.9,  2.3) #par[7]


# ### V2 (same as V1 but let p0 floating: weird results, need to investigate)
# scaleSetter.set(0 ,0, 0. , 0.002,   -0.01  , 0.01 )    # (1+par[0])
# scaleSetter.set(0, 0, 0., 0.02, -0.1, 0.1) #par[1]
# scaleSetter.set(0, 0, 0., 0.03, -0.15, 0.15) #par[2]
# scaleSetter.set(0, 0, 0., 0.05, -0.25, 0.25) #par[3]
# scaleSetter.set(0, 0, 0., 0.05, -0.25, 0.25) #par[4]
# scaleSetter.set(1, 0,   0.9,  0.1,   0.5,  1.1) #par[5]
# scaleSetter.set(1, 0,   1.6,  0.1,   1.3,  1.8) #par[6]
# scaleSetter.set(1, 0,   2.1,  0.1,   1.9,  2.3) #par[7]

# ### V3 (same as V0 but let p0 floating: weird results, need to investigate)
# scaleSetter.set(0 ,0, 0. , 0.002,   -0.01  , 0.01 )    # (1+par[0])
# scaleSetter.set(0, 0, -0.05, 0.01, -0.11, 0.01) #par[1]
# scaleSetter.set(0, 0, -0.07, 0.01, -0.15, 0.01) #par[2]
# scaleSetter.set(0, 0, -0.10, 0.02, -0.22, 0.02) #par[3]
# scaleSetter.set(0, 0, -0.13, 0.025, -0.25, 0.02) #par[4]
# scaleSetter.set(1, 0,   0.9,  0.1,   0.5,  1.1) #par[5]
# scaleSetter.set(1, 0,   1.6,  0.1,   1.3,  1.8) #par[6]
# scaleSetter.set(1, 0,   2.1,  0.1,   1.9,  2.3) #par[7]

# ### V4 (same as V2 but re-adjust ranges)
# scaleSetter.set(0 ,0, 0. , 0.002,   -0.01  , 0.01 )    # (1+par[0])
# scaleSetter.set(0, 0, 0., 0.06, -0.3, 0.3) #par[1]
# scaleSetter.set(0, 0, 0., 0.08, -0.4, 0.4) #par[2]
# scaleSetter.set(0, 0, 0., 0.1, -0.5, 0.5) #par[3]
# scaleSetter.set(0, 0, 0., 0.14, -0.7, 0.7) #par[4]
# scaleSetter.set(1, 0,   0.9,  0.1,   0.5,  1.1) #par[5]
# scaleSetter.set(1, 0,   1.6,  0.1,   1.3,  1.8) #par[6]
# scaleSetter.set(1, 0,   2.1,  0.1,   1.9,  2.3) #par[7]







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

options.register('tolerance',
                 0.1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Minuit tolerance")



options.parseArguments()


#### INPUT Tree: to be changed in case the input is the Tree ############################


#### SETTING THE VARIABLES #######################

#TREEINPUTNAME="root://eoscms//eos/cms/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/MC_8TeV/"
#TREEINPUTNAME="root://eoscms//eos/cms/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/Data2011_7TeV/"
TREEINPUTNAME="root://eoscms//eos/cms/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/MC_7TeV/"

TREEOUTPUTNAME=""

print "Running on sample _%s_ , selecting muons in the range: eta_mu1=[%f,%f], eta_mu2=[%f,%f]" % (options.sampleName,options.etaMin1,options.etaMax1,options.etaMin2,options.etaMax2)

# Use this when running on a tree
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(0))
TREEINPUTNAME+="zmumuTree_"+options.sampleName+".root"

    
process.looper = cms.Looper(
    "MuScleFit",
    # Only used when reading events from a root tree
    MaxEventsFromRootTree = cms.int32(-1),

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
    maxLoopNumber = cms.untracked.int32(5),
    # Select which fits to do in which loop (0 = do not, 1 = do)
    doResolFit =        cms.vint32(1,0,0,0,0),
    doScaleFit =        cms.vint32(0,1,1,1,0),
    doBackgroundFit =   cms.vint32(0,0,0,0,0),
    doCrossSectionFit = cms.vint32(0,0,0,0,0),

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
    ScaleFitType  = cms.int32(61),
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
    
    # type11: Z DATA 2012 Run2012C ReReco
    # type12: Z DATA 2011 Run2011AB 12Oct2013 ReReco (Legacy)
    
    BgrFitType = cms.vint32(12, 12, 12), 
    # These empty parameters should be used when there is no background
    parBgr = cms.vdouble(1., 1.,   1., 1.,   1., 1.,
                         1., 1.,   1., 1.,   1., 1.,   1., 1.,   1., 1.,   1., 1.),
    parBgrFix = cms.vint32(0, 0,   0, 0,   0, 0,
                           # The rest of the parameters is used for the resonance regions. They are automatically fixed in the code
                           # because they are never used to fit the background, but only after the rescaling.
                           1, 1,   1, 1,   1, 1,   1, 1,   1, 1,   1, 1),
    parBgrOrder = cms.vint32(0, 0,   0, 0,   0, 0,
                             0, 0,   0, 0,   0, 0,   0, 0,   0, 0,   0, 0),


    # ----------------------- #

    # Set Minuit fit strategy
    FitStrategy = cms.int32(2),
    MinuitTolerance = cms.untracked.double(options.tolerance),    

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
    
    RapidityBinsForZ = cms.untracked.bool(False),

    ProbabilitiesFile = cms.untracked.string("Probs_Z_1001_7TeV_v2.root"), # 7 Tev
    #ProbabilitiesFile = cms.untracked.string("Probs_Z_1001_8TeV_v2.root"), # 8 Tev

    # The following parameters can be used to filter events
    TriggerResultsLabel = cms.untracked.string("TriggerResults"),
    TriggerResultsProcess = cms.untracked.string("HLT"),
    TriggerPath = cms.untracked.vstring(""),
    # Negate the result of the trigger
    NegateTrigger = cms.untracked.bool(False),
    debug = cms.untracked.int32(0),
)
