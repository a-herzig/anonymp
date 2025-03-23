#!/usr/bin/dash
# First script of imputation process
# share an encoded genotype to compare server
# share encode keys to reference server
# must be executed from parent directory `/usr/bin/dash scripts/user_process_init.sh`
CORES="$1"
targets="$2"
chunks=$(gawk '{print $2}' res/coordinates.I51.2_15.txt)

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


cd outbox
touch 1-user-*
alias tarc='tar --create --use-compress-program="pigz -3" --remove-files'
tarc -f pack-1-user-2-reference.tar.gz 1-user-2-reference-*
tarc -f pack-1-user-3-compare.tar.gz 1-user-3-compare-*
tarc -f pack-1-user-4-ppm.tar.gz 1-user-4-ppm-*
tarc -f pack-1-user-5-product.tar.gz 1-user-5-product-*
tarc -f pack-1-user-1-user.tar.gz 1-user-1-user-*
