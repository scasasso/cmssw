#!/bin/bash 

dataset=DoubleMu_Run2011AB-12Oct2013-v1

#jobName as first option
jobName=$1
etaMax=$2
treeName=$3

SampleName=${treeName}

# Setup variables
LSF_DIR=/afs/cern.ch/work/s/scasasso/private/jobs/lsf
LOG_DIR=/afs/cern.ch/work/s/scasasso/private/jobs/log


OUT_DIR=/store/caf/user/scasasso/Alignment/ZMuMu/MuScleFit2.0/Results/Calibration/Data2011_7TeV/${jobName}/${dataset}
cmsMkdir $OUT_DIR


lsfname=${LSF_DIR}/${jobName}_${SampleName}.lsf
logname=${jobName}_${SampleName}.log

# 
# -- Clean lsf files
if [ -f $lsfname ]; then 
    rm  $lsfname
fi

if [ -f ${LOG_DIR}/${jobName}_${SampleName}.out.txt ]; then rm ${LOG_DIR}/${jobName}_${SampleName}.out.txt; fi
if [ -f ${LOG_DIR}/${jobName}_${SampleName}.err.txt ]; then rm ${LOG_DIR}/${jobName}_${SampleName}.err.txt; fi


# 
# -- Prepare the script to be submitted
touch $lsfname
echo '#!/bin/sh ' >> $lsfname
echo '#BSUB -L /bin/sh'    >> $lsfname
echo '#BSUB -R "select[type==SLC6_64]"'   >>$lsfname
echo '#BSUB -u stefano.casasso@cern.ch'    >> $lsfname
echo '#BSUB -J ' ${jobName}   >> $lsfname
# echo '#BSUB -oo' ${LOG_DIR}/${jobName}_${SampleName}.log >> $lsfname
echo '#BSUB -o ' ${LOG_DIR}/${jobName}_${SampleName}.out.txt   >> $lsfname
echo '#BSUB -e ' ${LOG_DIR}/${jobName}_${SampleName}.err.txt   >> $lsfname
echo '#BSUB -q cmscaf1nw'    >> $lsfname
echo 'jobName='$jobName >> $lsfname
echo 'SampleName='$SampleName >> $lsfname
echo 'logname='$logname >> $lsfname
echo 'etaMax='$etaMax >> $lsfname
echo 'treeName='$treeName >> $lsfname
echo 'OUT_DIR='$OUT_DIR >> $lsfname

cat >>$lsfname <<"LSF"
echo  -----------------------
echo  Job started at `date`
echo  -----------------------

#Some variable definitions
echo "Defining variables ..."
cmssw_ver=CMSSW_5_3_9
CMSSW_DIR=/afs/cern.ch/work/s/scasasso/private/MuScleFit_dev4/${cmssw_ver}/src
TEST_DIR=$CMSSW_DIR/MuonAnalysis/MomentumScaleCalibration/test

# Save current dir on lxbatch machine
LXBATCH_DIR=`pwd` 

# Setup variables
echo "Setup CAF environment (is this really needed??) ..."
. /afs/cern.ch/cms/caf/setup.sh


echo "Setup CMSSW release ..."
scram p CMSSW ${cmssw_ver}
cd ${cmssw_ver}/src/
eval `scram r -sh`

echo "Inspecting ${LXBATCH_DIR}/${cmssw_ver}/src:"
ls -lart

echo "Gymnastics to make Git work on lxbatch ..."
#Some gymnastics to make Git work on lxbatch
git config --global user.name "Stefano Casasso"
git config --global user.email stefano.casasso@cern.ch
git config --global user.github scasasso
git config --global http.sslVerify false

#git clone git://github.com/fwyzard/cms-git-tools
#PATH=${PWD}/cms-git-tools:$PATH 
git cms-init --ssh


echo "Checkout MuScleFit from repo ..."
#Co the package
git cms-addpkg MuonAnalysis/MomentumScaleCalibration 
#git clone https://github.com/scasasso/cmssw
git remote add upstream https://github.com/scasasso/cmssw.git
git fetch upstream
git checkout upstream/test_binned_function


# echo "Update Functions.* files with the latest version (not in the repo) ..."
# cp $CMSSW_DIR/MuonAnalysis/MomentumScaleCalibration/src/Functions.cc MuonAnalysis/MomentumScaleCalibration/src
# cp $CMSSW_DIR/MuonAnalysis/MomentumScaleCalibration/interface/Functions.h MuonAnalysis/MomentumScaleCalibration/interface
grep -A2 -i "scalefunctiontype53" MuonAnalysis/MomentumScaleCalibration/src/Functions.cc
grep -A2 -i "scalefunctiontype53" MuonAnalysis/MomentumScaleCalibration/interface/Functions.h


# #Compile
echo "Compile ..."
scram b -j 20

echo "Get into test directory ..."
#Go to the test directory and prepare to start the job
cd MuonAnalysis/MomentumScaleCalibration/test

echo "Copy cfg file ..."
cp ${TEST_DIR}/CalibrationFromTree_${jobName}_cfg.py ./dummy_cfg.py 
ls -larth dummy_cfg.py

echo "Copying propbability ROOT files ..."
cp ${TEST_DIR}/Probs_Z_1001_8TeV_v2.root .
cp ${TEST_DIR}/Probs_Z_1001_7TeV_v2.root .

echo 'Content of test directory just before running cmsRun:'
ls -lh .

echo "Issue cmsRun ..."
cmsRun dummy_cfg.py etaMin1=-${etaMax} etaMax1=${etaMax} etaMin2=-${etaMax} etaMax2=${etaMax} sampleName=${treeName} | tee ${logname}

echo "Inspecting test directory after cmsRun ..."
ls -lh . 


# Copy output files to EOS 
echo "Copy output files to EOS ..."
for RootOutputFile in $(ls *mumuHisto*root ); do 
    echo " ... ${RootOutputFile}"
    cmsStage -f ${RootOutputFile} $OUT_DIR/; 
done

echo "... dummy_cfg.py"
cmsStage -f dummy_cfg.py ${OUT_DIR}/CalibrationFromTree_${jobName}_cfg.py 

for file in $(ls *log *txt); do 
    echo "... ${file}"
    cmsStage -f ${file} $OUT_DIR/; 
done

for file in $(ls *out); do
    gzip ${file}
    echo "... ${file}"
    cmsStage -f ${file}.gz $OUT_DIR/
done


echo  -----------------------
echo  Job ended at `date`
echo  -----------------------    
LSF

    
chmod u+x $lsfname    
bsub < $lsfname 	
