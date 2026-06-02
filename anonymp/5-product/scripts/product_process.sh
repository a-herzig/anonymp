#!/usr/bin/dash
# Sixth script of imputation process
# must be executed from parent directory `/usr/bin/dash scripts/product_process.sh`

CORES="$1"
ENCRYPTED_MESSAGES="$3"
SECURE_SUMMATION="$4"

mkdir -p outbox tmp

echo process imputation product
parallel --max-procs $CORES Rscript scripts/product_process.R :::: inbox/1-user-5-product-chunks.txt

touch outbox/5-product-*
