#!/usr/bin/dash
# Fourth script of imputation process
# compute ppm matrix and share it to imputation server
# must be executed from parent directory `/usr/bin/dash scripts/ppm_process.sh`
CORES="$1"
CHUNKS_PATH=inbox/1-user-4-ppm-chunks.txt
mkdir -p tmp outbox

rand_int32() {
  rand_uint32=$(echo "ibase=16; $(openssl rand --hex 4 | tr [:lower:] [:upper:])" | bc)
  if [ $rand_uint32 -ge $(echo "2^31" | bc) ]
  then
    output=$(echo "$rand_uint32 - 2^32" | bc)
  else
    output=$rand_uint32
  fi
  echo $output
}

cd inbox
alias tarx='tar --extract --use-compress-program="pigz -3"'
tarx -f pack-1-user-4-ppm.tar.gz
tarx -f pack-2-reference-4-ppm.tar.gz
tarx -f pack-3-compare-4-ppm.tar.gz

cd ..
for chunk_name in $(cat $CHUNKS_PATH)
do
  echo "$(rand_int32)\n$(rand_int32)" > tmp/rand-chunk${chunk_name}.txt
  # share seed to compute full_shuffle_key to 2-reference
  cp tmp/rand-chunk${chunk_name}.txt outbox/4-ppm-2-reference-rand-chunk${chunk_name}.txt
done

parallel --max-procs $CORES Rscript scripts/ppm_process.R :::: $CHUNKS_PATH

cd outbox
touch 4-ppm-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-4-ppm-2-reference.tar.gz 4-ppm-2-reference-*
tarc -f pack-4-ppm-5-product.tar.gz 4-ppm-5-product-*
tarc -f pack-4-ppm-1-user.tar.gz 4-ppm-1-user-*
