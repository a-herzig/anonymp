#!/usr/bin/dash
# Fifth script of imputation process
# share an encoded reference haplotype to imputation product server
# must be executed from parent directory `/usr/bin/dash scripts/reference_process_final.sh`

CORES="$1"

mkdir -p outbox tmp

cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-2-reference-2-reference.tar.gz
tarx -f pack-4-ppm-2-reference.tar.gz

cd ..
echo process reference haplotypes
parallel --max-procs $CORES Rscript scripts/reference_process_final.R :::: inbox/1-user-2-reference-chunks.txt

cd outbox
touch 2-reference-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-2-reference-5-product.tar.gz 2-reference-5-product-*
tarc -f pack-2-reference-1-user.tar.gz 2-reference-1-user-*
