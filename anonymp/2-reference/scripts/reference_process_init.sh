#!/usr/bin/bash
# Second script of imputation process
# share an encoded reference haplotype to compare server
# must be executed from parent directory `/usr/bin/dash scripts/reference_process_init.sh`

CORES="$1"
ENCRYPTED_MESSAGES="$3"

mkdir -p outbox tmp

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    USER_REFERENCE_PWD="user_reference"
    USER_REFERENCE_PREFIX="1-user-2-reference"
    dest=USER_REFERENCE
    pwd_var=${dest}_PWD
    prefix_var=${dest}_PREFIX
    for file in $(find inbox -name "${!prefix_var}*.aes128")
    do
	cat $file | openssl enc -aes128 -pbkdf2 -a -d -k "${!pwd_var}" > "${file%.aes128}"
    done
fi

echo extract positions from reference haplotypes
zcat res/chr15_5popSim_4B11_Ref.vcf.gz \
  | gawk 'FNR > 10 { print $2 }' > tmp/positions.txt

echo process reference haplotypes
parallel --max-procs $CORES Rscript scripts/reference_process_init.R :::: inbox/1-user-2-reference-chunks.txt

if [ "$ENCRYPTED_MESSAGES" = true ] ; then
    REFERENCE_USER_PWD="reference_user"
    REFERENCE_COMPARE_PWD="reference_compare"
    REFERENCE_PPM_PWD="reference_ppm"

    REFERENCE_USER_PREFIX="2-reference-1-user"
    REFERENCE_COMPARE_PREFIX="2-reference-3-compare"
    REFERENCE_PPM_PREFIX="2-reference-4-ppm"

    for dest in REFERENCE_USER REFERENCE_COMPARE REFERENCE_PPM
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
