#!/bin/bash
LISTFILE="/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/filelists/ionbeta_file_list.txt"
NFILES=$(cat ${LISTFILE} | wc -l)

declare -a size
while IFS= read -r line
do
    size+=($line)
    echo $line
done < "$LISTFILE"

echo "Making trees from " $NFILES " lmd files."

sbatch -J ionbeta \
--cpus-per-task=2 \
--mem-per-cpu=8G \
--array=1-$NFILES \
-o /lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/cluster/logs/ionbeta_%A_%a.out.log \
-e /lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/cluster/logs/ionbeta_%A_%a.err.log \
-- /lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/cluster/launchers/ionbeta_launcher.sh

unset size
