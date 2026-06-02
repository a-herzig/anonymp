#!/usr/bin/bash
# First script of imputation process
# share an encoded genotype to compare server
# share encode keys to reference server
# must be executed from parent directory `/usr/bin/dash scripts/user_process_init.sh`
CORES="$1"
targets="$2"
ENCRYPTED_MESSAGES="$3"
SECURE_SUMMATION="$4"

chunks=$(gawk '{print $1}' res/coordinates.I51.2_15.txt)

mkdir -p tmp outbox
for target in $targets
do
  for chunk in $chunks
  do
    code=$(od /dev/urandom --address-radix=n --format=x4 --read-bytes=4 | tr -d '[:blank:]')
    echo $code >> tmp/chunks_anon.txt
    echo $code $target $chunk >> tmp/chunks.txt
  done
done
sort tmp/chunks_anon.txt -o tmp/chunks_anon.txt
sort tmp/chunks.txt -o tmp/chunks.txt

cp tmp/chunks_anon.txt outbox/1-user-1-user-chunks.txt
cp tmp/chunks_anon.txt outbox/1-user-2-reference-chunks.txt
cp tmp/chunks_anon.txt outbox/1-user-3-compare-chunks.txt
cp tmp/chunks_anon.txt outbox/1-user-4-ppm-chunks.txt
cp tmp/chunks_anon.txt outbox/1-user-5-product-chunks.txt

if [ -n "$(uniq tmp/chunks_anon.txt --repeated)" ]
then
  echo "FATAL : same random name generated several times"
  echo "This script can be run again safely. "
  exit 1
fi

parallel --max-procs $CORES --colsep ' ' Rscript scripts/user_process_init.R :::: tmp/chunks.txt

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    USER_REFERENCE_PWD="user_reference"
    USER_COMPARE_PWD="user_compare"
    USER_PPM_PWD="user_ppm"
    USER_PRODUCT_PWD="user_product"

    USER_REFERENCE_PREFIX="1-user-2-reference"
    USER_COMPARE_PREFIX="1-user-3-compare"
    USER_PPM_PREFIX="1-user-4-ppm"
    USER_PRODUCT_PREFIX="1-user-5-product"

    for dest in USER_REFERENCE USER_COMPARE USER_PPM USER_PRODUCT
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

touch outbox/1-user-*
