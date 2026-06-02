#!/usr/bin/bash
# Sixth script of imputation process
# must be executed from parent directory `/usr/bin/dash scripts/product_process.sh`

CORES="$1"
ENCRYPTED_MESSAGES="$3"
SECURE_SUMMATION="$4"

mkdir -p outbox tmp

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    USER_PRODUCT_PWD="user_product"
    REFERENCE_PRODUCT_PWD="reference_product"
    PPM_PRODUCT_PWD="ppm_product"

    USER_PRODUCT_PREFIX="1-user-5-product"
    REFERENCE_PRODUCT_PREFIX="2-reference-5-product"
    PPM_PRODUCT_PREFIX="4-ppm-5-product"

    for dest in USER_PRODUCT REFERENCE_PRODUCT PPM_PRODUCT
    do
	for file in $(find inbox -name "${!prefix_var}*.aes128")
	do
	    cat $file | openssl enc -aes128 -pbkdf2 -a -d -k "${!pwd_var}" > "${file%.aes128}"
	done
    done
fi

echo process imputation product
parallel --max-procs $CORES Rscript scripts/product_process.R :::: inbox/1-user-5-product-chunks.txt

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    PRODUCT_USER_PWD="product_user"
    PRODUCT_USER_PREFIX="5-product-1-user"
    dest=PRODUCT_USER
    pwd_var=${dest}_PWD
    prefix_var=${dest}_PREFIX

    for file in $(find outbox -name "${!prefix_var}*")
    do
	cat $file | openssl enc -aes128 -pbkdf2 -a -e -k "${!pwd_var}" > "${file}.aes128"
	rm $file
    done
fi

touch outbox/5-product-*
