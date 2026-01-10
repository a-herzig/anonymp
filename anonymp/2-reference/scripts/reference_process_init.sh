#!/usr/bin/dash
# Second script of imputation process
# share an encoded reference haplotype to compare server
# must be executed from parent directory `/usr/bin/dash scripts/reference_process_init.sh`

CORES="$1"

mkdir -p outbox tmp

echo extract positions from reference haplotypes
zcat res/chr15_5popSim_4B11_Ref.vcf.gz \
  | gawk 'FNR > 10 { print $2 }' > tmp/positions.txt

echo process reference haplotypes
parallel --max-procs $CORES Rscript scripts/reference_process_init.R :::: inbox/1-user-2-reference-chunks.txt

cd outbox
touch 2-reference-*
