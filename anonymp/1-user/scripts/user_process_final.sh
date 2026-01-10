#!/usr/bin/dash
# Eighth script of imputation process
# decode imputation and compute score
# must be executed from parent directory `/usr/bin/dash scripts/user_process_final.sh`
CORES="$1"

parallel --max-procs $CORES Rscript scripts/user_process_final.R :::: inbox/1-user-1-user-chunks.txt
