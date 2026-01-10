#!/usr/bin/dash
# Fifth script of imputation process
# share an encoded reference haplotype to imputation product server
# must be executed from parent directory `/usr/bin/dash scripts/reference_process_final.sh`

CORES="$1"

mkdir -p outbox tmp

echo process reference haplotypes
parallel --max-procs $CORES Rscript scripts/reference_process_final.R :::: inbox/1-user-2-reference-chunks.txt

touch outbox/2-reference-*
