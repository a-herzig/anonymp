#!/usr/bin/dash
# Third script of imputation process
# share a comparison matrix to ppm server
# must be executed from parent directory `/usr/bin/dash scripts/compare_process.sh`
CORES="$1"
mkdir -p outbox

parallel --max-procs $CORES Rscript scripts/compare_process.R :::: inbox/1-user-3-compare-chunks.txt

touch outbox/3-compare-*
