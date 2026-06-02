#!/usr/bin/bash
# Eighth script of imputation process
# decode imputation and compute score
# must be executed from parent directory `/usr/bin/dash scripts/user_process_final.sh`
CORES="$1"
ENCRYPTED_MESSAGES="$3"
SECURE_SUMMATION="$4"

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    REFERENCE_USER_PWD="reference_user"
    PPM_USER_PWD="ppm_user"
    PRODUCT_USER_PWD="product_user"

    REFERENCE_USER_PREFIX="2-reference-1-user"
    PPM_USER_PREFIX="4-ppm-1-user"
    PRODUCT_USER_PREFIX="5-product-1-user"

    for dest in REFERENCE_USER PPM_USER PRODUCT_USER
    do
	pwd_var=${dest}_PWD
	prefix_var=${dest}_PREFIX
	for file in $(find inbox -name "${!prefix_var}*.aes128")
	do
	    cat $file | openssl enc -aes128 -pbkdf2 -a -d -k "${!pwd_var}" > "${file%.aes128}"
	done
    done
fi

parallel --max-procs $CORES Rscript scripts/user_process_final.R :::: inbox/1-user-1-user-chunks.txt
