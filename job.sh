#!/bin/bash
#Source
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97/Geant4/10.06.p01/x86_64-centos7-gcc8-opt/Geant4-env.sh

source /cvmfs/sft.cern.ch/lcg/releases/LCG_97/Geant4/10.06.p01/x86_64-centos7-gcc8-opt/bin/geant4.sh
# 3 = Energie, 4 = targtick , 5 = absthick, 6 = number of e-/bunch 

# create the mac file 
sed "s/Energy/${3}/g" run_batch.tmp > run_${1}.mac
sed -i "s/targthick/${4}/g"  run_${1}.mac
sed -i "s/absthick/${5}/g"   run_${1}.mac
sed -i "s/enum/${6}/g"       run_${1}.mac 


./leap_sims -m run_${1}.mac -r results_${2}_${1}.root

# move the root files to results/ and the run_*.mac to macfiles/


 mv results_${2}_${1}.root results 
 mv run_${1}.mac macfiles 



