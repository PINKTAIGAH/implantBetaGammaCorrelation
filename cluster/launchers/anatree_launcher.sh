#!/bin/bash

export FAIRROOTPATH=/cvmfs/fairsoft.gsi.de/debian11/fairroot/v18.8.0_nov22p1
export SIMPATH=/cvmfs/fairsoft.gsi.de/debian11/fairsoft/nov22p1
export UCESB_DIR=/lustre/gamma/gbrunic/G302/ucesb
export UCESB_BASE_DIR=/lustre/gamma/gbrunic/G302/ucesb

# Set up environment variables for ROOT, FAIRROOT, and other dependencies
export ROOTSYS="${FAIRROOTPATH}"
export PATH="${FAIRROOTPATH}/bin:${SIMPATH}/bin:${UCESB_DIR}/bin:${PATH}"
export LD_LIBRARY_PATH="${FAIRROOTPATH}/lib:${SIMPATH}/lib:${UCESB_DIR}/lib:${LD_LIBRARY_PATH}"

# Source setup scripts
. "${FAIRROOTPATH}/bin/FairRootConfig.sh"
. "${SIMPATH}/bin/thisroot.sh"
. "/lustre/gamma/gbrunic/G302/build/config.sh"

#for some reasion thisroot.sh seems to unset some of these guys:
export FAIRROOTPATH=/cvmfs/fairsoft.gsi.de/debian11/fairroot/v18.8.0_nov22p1
export SIMPATH=/cvmfs/fairsoft.gsi.de/debian11/fairsoft/nov22p1
export UCESB_DIR=/lustre/gamma/gbrunic/G302/ucesb
export UCESB_BASE_DIR=/lustre/gamma/gbrunic/G302/ucesb

LISTFILE="/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/filelists/tree2anatree_file_list.txt"
NFILES=$(cat ${LISTFILE} | wc -l)
echo "Analysing" $NFILES " files"

declare -a size
while IFS= read -r line
do
    size+=($line)
done < "$LISTFILE"

index=$(($SLURM_ARRAY_TASK_ID - 1))
infile_dir="${size[index]}"
infile_basename=$(basename ${infile_dir})
outfile_dir="/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/anatrees/${infile_basename//tree.root/anatree.root}"

echo $infile_dir
echo $outfile_dir

echo "Running c4! <3"
root -b -l <<EOF

gSystem->AddIncludePath("${FAIRROOTPATH}/include");
gSystem->AddIncludePath("${SIMPATH}/include");

gSystem->AddLinkedLibs("-L/lustre/gamma/gbrunic/G302/build/lib -llibc4source.so");
gSystem->AddLinkedLibs("-L/lustre/gamma/gbrunic/G302/build/lib -llibc4Analysis.so");  
gSystem->AddLinkedLibs("-L/lustre/gamma/gbrunic/G302/build/lib -llibc4Data.so"); 
gSystem->AddLinkedLibs("-L/lustre/gamma/gbrunic/G302/build/lib -llibc4MacroCompiler.so"); 
gSystem->AddLinkedLibs("-L/lustre/gamma/gbrunic/G302/build/lib -llibc4Base.so"); 

.x /lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/scripts/makeAnatrees.C("$infile_dir", "$outfile_dir")
EOF

unset size
