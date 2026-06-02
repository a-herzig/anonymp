#!/usr/bin/bash
# Third script of imputation process
# share a comparison matrix to ppm server
# must be executed from parent directory `/usr/bin/dash scripts/compare_process.sh`
CORES="$1"
ENCRYPTED_MESSAGES="$3"

mkdir -p outbox

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    USER_COMPARE_PWD="user_compare"
    REFERENCE_COMPARE_PWD="reference_compare"

    USER_COMPARE_PREFIX="1-user-3-compare"
    REFERENCE_COMPARE_PREFIX="2-reference-3-compare"

    for dest in USER_COMPARE REFERENCE_COMPARE
    do
	pwd_var=${dest}_PWD
	prefix_var=${dest}_PREFIX
	for file in $(find inbox -name "${!prefix_var}*.aes128")
	do
	    cat $file | openssl enc -aes128 -pbkdf2 -a -d -k "${!pwd_var}" > "${file%.aes128}"
	done
    done
fi

parallel --max-procs $CORES Rscript scripts/compare_process.R :::: inbox/1-user-3-compare-chunks.txt

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    COMPARE_PPM_PWD="compare_ppm"
    COMPARE_PPM_PREFIX="3-compare-4-ppm"
    dest=COMPARE_PPM
    pwd_var=${dest}_PWD
    prefix_var=${dest}_PREFIX
    for file in $(find outbox -name "${!prefix_var}*")
    do
	cat $file | openssl enc -aes128 -pbkdf2 -a -e -k "${!pwd_var}" > "${file}.aes128"
	rm $file
    done
fi

touch outbox/3-compare-*
