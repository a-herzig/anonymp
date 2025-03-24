#!/usr/bin/dash
# Sixth script of imputation process
# must be executed from parent directory `/usr/bin/dash scripts/product_process.sh`

CORES="$1"

mkdir -p outbox tmp

cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-1-user-5-product.tar.gz
tarx -f pack-4-ppm-5-product.tar.gz
tarx -f pack-2-reference-5-product.tar.gz

cd ..
echo process imputation product
parallel --max-procs $CORES Rscript scripts/product_process.R :::: inbox/1-user-5-product-chunks.txt

cd outbox
touch 5-product-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-5-product-1-user.tar.gz 5-product-1-user-*
