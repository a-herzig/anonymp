#!/usr/bin/dash
# Eighth script of imputation process
# decode imputation and compute score
# must be executed from parent directory `/usr/bin/dash scripts/user_process_final.sh`
CORES="$1"
mkdir -p tmp outbox


parallel --max-procs $CORES Rscript scripts/summation_process.R :::: inbox/1-user-6-summation-chunks.txt

cd outbox
touch 6-summation-*
