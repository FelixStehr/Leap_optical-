#!/bin/bash
#Source
source /cvmfs/sft.cern.ch/lcg/releases/LCG_97/Geant4/10.06.p01/x86_64-centos7-gcc8-opt/Geant4-env.sh

source /cvmfs/sft.cern.ch/lcg/releases/LCG_97/Geant4/10.06.p01/x86_64-centos7-gcc8-opt/bin/geant4.sh


# Simulation 
cd /afs/desy.de/group/flc/pool/fstehr/Simulation_Felix/leap_sims_Felix_neu/ActionClasses/leap_sims_from_basic_ActionClassBatch/build

./leap_sims -m run_batch.mac -r hallo1.root

# move the root files to results/
 cd /afs/desy.de/group/flc/pool/fstehr/Simulation_Felix/leap_sims_Felix_neu/ActionClasses/leap_sims_from_basic_ActionClassBatch/build/
 mv hallo1.root results  
 cd results 

# merge the root files to one and delete the old 
 hadd HalloWelt.root hallo1.root
rm hallo1.root

# 
