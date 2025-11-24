#!/usr/bin/dash
# Eighth script of imputation process
# decode imputation and compute score
# must be executed from parent directory `/usr/bin/dash scripts/user_process_final.sh`
CORES="$1"
mkdir -p tmp outbox


cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-1-user-6-summation.tar.gz
tarx -f pack-4-ppm-6-summation.tar.gz
tarx -f pack-5-product-6-summation.tar.gz

cd ..
parallel --max-procs $CORES Rscript scripts/summation_process.R :::: inbox/1-user-6-summation-chunks.txt

cd outbox
touch 6-summation-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-6-summation-1-user.tar.gz 6-summation-1-user-*