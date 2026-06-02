#!/usr/bin/bash
# Fifth script of imputation process
# share an encoded reference haplotype to imputation product server
# must be executed from parent directory `/usr/bin/dash scripts/reference_process_final.sh`

CORES="$1"
ENCRYPTED_MESSAGES="$3"

mkdir -p outbox tmp

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    PPM_REFERENCE_PWD="ppm_reference"
    PPM_REFERENCE_PREFIX="4-ppm-2-reference"

    for dest in PPM_REFERENCE
    do
	pwd_var=${dest}_PWD
	prefix_var=${dest}_PREFIX
	for file in $(find inbox -name "${!prefix_var}*.aes128")
	do
	    cat $file | openssl enc -aes128 -pbkdf2 -a -d -k "${!pwd_var}" > "${file%.aes128}"
	done
    done
fi

echo process reference haplotypes
parallel --max-procs $CORES Rscript scripts/reference_process_final.R :::: inbox/1-user-2-reference-chunks.txt

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    REFERENCE_USER_PWD="reference_user"
    REFERENCE_PRODUCT_PWD="reference_product"

    REFERENCE_USER_PREFIX="2-reference-1-user"
    REFERENCE_PRODUCT_PREFIX="2-reference-5-product"

    for dest in REFERENCE_USER REFERENCE_PRODUCT
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

touch outbox/2-reference-*
