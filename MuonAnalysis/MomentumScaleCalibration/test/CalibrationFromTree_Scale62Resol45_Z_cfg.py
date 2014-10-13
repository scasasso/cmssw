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
### V0 (beta version)
scaleSetter = SetParameters()
#fix, order, value, step, min, max
scaleSetter.set(1, 0, 0., 0.0001, -0.01, 0.01) #(1+par[0])
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[1]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[2]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[3]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[4]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[5]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[6]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[7]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[8]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[9]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[10]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[11]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[12]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[13]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[14]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[15]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[16]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[17]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[18]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[19]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[20]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[21]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[22]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[23]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[24]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[25]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[26]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[27]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[28]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[29]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[30]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[31]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[32]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[33]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[34]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[35]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[36]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[37]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[38]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[39]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[40]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[41]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[42]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[43]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[44]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[45]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[46]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[47]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[48]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[49]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[50]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[51]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[52]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[53]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[54]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[55]
scaleSetter.set(0, 0, 0., 0.001, -0.005, 0.005) #par[56]


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
TREEINPUTNAME="root://eoscms//eos/cms/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/MC_7TeV/"
#TREEINPUTNAME="root://eoscms//eos/cms/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/Data2011_7TeV/"

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
    ScaleFitType  = cms.int32(62),
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
    # type13: Z MC 2011 Summer11LegacyDR ReReco (Legacy)
    
    BgrFitType = cms.vint32(13, 13, 13), 
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
