#!/usr/bin/dash
# Eighth script of imputation process
# decode imputation and compute score
# must be executed from parent directory `/usr/bin/dash scripts/user_process_final.sh`
CORES=16

cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-1-user-1-user.tar.gz
tarx -f pack-2-reference-1-user.tar.gz
tarx -f pack-4-ppm-1-user.tar.gz
tarx -f pack-5-product-1-user.tar.gz

cd ..
parallel --max-procs $CORES Rscript scripts/user_process_final.R :::: inbox/1-user-1-user-chunks.txt
