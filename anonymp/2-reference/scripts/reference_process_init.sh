#!/usr/bin/dash
# Second script of imputation process
# share an encoded reference haplotype to compare server
# must be executed from parent directory `/usr/bin/dash scripts/reference_process_init.sh`

CORES=16

mkdir -p outbox tmp

cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-1-user-2-reference.tar.gz

cd ..
echo extract positions from reference haplotypes
zcat res/chr15_5popSim_4B11_Ref.vcf.gz \
  | gawk 'FNR > 10 { print $2 }' > tmp/positions.txt

echo process reference haplotypes
parallel --max-procs $CORES Rscript scripts/reference_process_init.R :::: inbox/1-user-2-reference-chunks.txt

cd outbox
touch 2-reference-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-2-reference-3-compare.tar.gz 2-reference-3-compare-*
tarc -f pack-2-reference-4-ppm.tar.gz 2-reference-4-ppm-*
tarc -f pack-2-reference-2-reference.tar.gz 2-reference-2-reference-*
tarc -f pack-2-reference-1-user.tar.gz 2-reference-1-user-*
