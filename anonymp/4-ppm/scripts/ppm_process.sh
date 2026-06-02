#!/usr/bin/bash
# Fourth script of imputation process
# compute ppm matrix and share it to imputation server
# must be executed from parent directory `/usr/bin/dash scripts/ppm_process.sh`
CORES="$1"
ENCRYPTED_MESSAGES="$3"
SECURE_SUMMATION="$4"

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    USER_PPM_PWD="user_ppm"
    REFERENCE_PPM_PWD="reference_ppm"
    COMPARE_PPM_PWD="compare_ppm"

    USER_PPM_PREFIX="1-user-4-ppm"
    REFERENCE_PPM_PREFIX="2-reference-4-ppm"
    COMPARE_PPM_PREFIX="3-compare-4-ppm"

    for dest in USER_PPM REFERENCE_PPM COMPARE_PPM
    do
	pwd_var=${dest}_PWD
	prefix_var=${dest}_PREFIX
	for file in $(find inbox -name "${!prefix_var}*.aes128")
	do
	    cat $file | openssl enc -aes128 -pbkdf2 -a -d -k "${!pwd_var}" > "${file%.aes128}"
	done
    done
fi

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

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    PPM_USER_PWD="ppm_user"
    PPM_REFERENCE_PWD="ppm_reference"
    PPM_PRODUCT_PWD="ppm_product"

    PPM_USER_PREFIX="4-ppm-1-user"
    PPM_REFERENCE_PREFIX="4-ppm-2-reference"
    PPM_PRODUCT_PREFIX="4-ppm-5-product"

    for dest in PPM_USER PPM_REFERENCE PPM_PRODUCT
    do
	pwd_var=${dest}_PWD
	prefix_var=${dest}_PREFIX
	for file in $(find outbox -name "${!prefix_var}*")
	do
	    cat $file | openssl enc -aes128 -pbkdf2 -a -e -k "${!pwd_var}" > "${file}.aes128"
	    rm $file
	done
    done
fi

touch outbox/4-ppm-*
