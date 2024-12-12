#!/usr/bin/dash
# Third script of imputation process
# share a comparison matrix to ppm server
# must be executed from parent directory `/usr/bin/dash scripts/compare_process.sh`
CORES=16
mkdir -p outbox

cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-1-user-3-compare.tar.gz
tarx -f pack-2-reference-3-compare.tar.gz

cd ..
parallel --max-procs $CORES Rscript scripts/compare_process.R :::: inbox/1-user-3-compare-chunks.txt

cd outbox
touch 3-compare-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-3-compare-4-ppm.tar.gz 3-compare-4-ppm-*
