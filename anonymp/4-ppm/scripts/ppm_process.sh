#!/usr/bin/dash
# Fourth script of imputation process
# compute ppm matrix and share it to imputation server
# must be executed from parent directory `/usr/bin/dash scripts/ppm_process.sh`
CORES="$1"
ENCRYPTED_MESSAGES="$3"
SECURE_SUMMATION="$4"

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

for chunk_name in $(cat $CHUNKS_PATH)
do
  echo "$(rand_int32)" > tmp/rand-chunk${chunk_name}.txt
  echo "$(rand_int32)" >> tmp/rand-chunk${chunk_name}.txt
  # share seed to compute full_shuffle_key to 2-reference
  cp tmp/rand-chunk${chunk_name}.txt outbox/4-ppm-2-reference-rand-chunk${chunk_name}.txt
done

parallel --max-procs $CORES Rscript scripts/ppm_process.R :::: $CHUNKS_PATH

touch outbox/4-ppm-*
