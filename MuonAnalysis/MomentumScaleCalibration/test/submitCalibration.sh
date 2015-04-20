#!/bin/bash 

cfgName=CalibrationFromTree_Scale50Resol45_Z_cfg.py
jobName=ZMC_Scale50Resol45
dataset=MC_13TeV
sampleName=DYToMuMu_13TeV_Phys14DR-PU20bx25
etaMax=2.4
append=""

echo "Running job ${jobName}:"
echo "  cfg: ${cfgName}"
echo "  sampleName: ${sampleName}"

# Setup variables
LSF_DIR=/afs/cern.ch/work/s/scasasso/public/jobs/lsf
LOG_DIR=/afs/cern.ch/work/s/scasasso/public/jobs/log

if [ "${4}" == "" ]; then 
    append=$4
else append=_${4}
fi
    
OUT_DIR=/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/Results/Calibration/${dataset}/${jobName}/${sampleName}
echo "Output directory will be ${OUT_DIR}"
cmsMkdir $OUT_DIR


lsfname=${LSF_DIR}/${jobName}_${sampleName}.lsf
logname=${jobName}_${sampleName}.log

stdOutName=${LOG_DIR}/${jobName}_${sampleName}${append}.out.txt
stdErrName=${stdOutName/out/err}

# 
# -- Clean lsf files
if [ -f $lsfname ]; then 
    rm  $lsfname
fi

if [ -f  ${stdOutName} ]; then rm ${stdOutName}; fi
if [ -f ${stdErrName} ]; then rm ${stdErrName}; fi


# 
# -- Prepare the script to be submitted
touch $lsfname
echo '#!/bin/bash ' >> $lsfname
echo '#BSUB -L /bin/bash'    >> $lsfname
echo '#BSUB -u scasasso@cern.ch'    >> $lsfname
echo '#BSUB -J ' ${jobName}   >> $lsfname
echo '#BSUB -q 1nh'    >> $lsfname
echo 'cfgName='$cfgName >> $lsfname
echo 'jobName='$jobName >> $lsfname
echo 'sampleName='$sampleName >> $lsfname
echo 'logname='$logname >> $lsfname
echo 'etaMax='$etaMax >> $lsfname
echo 'attend='$attend >> $lsfname
echo 'OUT_DIR='$OUT_DIR >> $lsfname

cat >>$lsfname <<"LSF"
echo  -----------------------
echo  Job started at `date`
echo  -----------------------

CMSSW_DIR=/afs/cern.ch/work/s/scasasso/private/MuScleFit_tutorial/CMSSW_7_2_4/src
TEST_DIR=$CMSSW_DIR/MuonAnalysis/MomentumScaleCalibration/test

LXBATCH_DIR=`pwd` 

cd ${CMSSW_DIR}
eval `scram r -sh`

cd ${TEST_DIR}

cmsRun ${cfgName} etaMin1=-${etaMax} etaMax1=${etaMax} etaMin2=-${etaMax} etaMax2=${etaMax} sampleName=${sampleName} 2>&1 | tee ${logname}


for RootOutputFile in $(ls *mumuHisto*root ); do cmsStage -f ${RootOutputFile} $OUT_DIR; done
cmsStage -f ${cfgName} ${OUT_DIR} 
for file in $(ls *log *txt); do cmsStage -f ${file} $OUT_DIR/; done
for file in $(ls *out); do gzip ${file}; cmsStage -f ${file}.gz $OUT_DIR; done


echo  -----------------------
echo  Job ended at `date`
echo  -----------------------    
LSF

    
chmod u+x $lsfname    
bsub -o ${stdOutName} -e ${stdErrName} < $lsfname 	
