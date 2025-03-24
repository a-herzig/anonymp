#!/usr/bin/dash

basedir=$(pwd)
mkdir -p "$basedir/tmp"
prefix="$basedir/tmp/$(date '+%Y%m%d-%H%M%S-')"
if [ "$#" -eq 0 ]
then
  echo "pipeline require target index to impute"
  exit 2
fi
cores="$1"
targets="$2"

step_total_duration_file="${prefix}duration_step_total.txt"

runstep () {
  echo process $1
  start=$(date +%s.%N)
  $2 $cores "$3"
  duration=$(echo "$(date +%s.%N) - $start" | bc)
  echo "$1 $duration" >> $step_total_duration_file
}

dispatch () {
  for dest in $2
  do
    mkdir -p $basedir/$dest/inbox/
    mv --verbose $basedir/$1/outbox/pack-$1-$dest.tar.gz $basedir/$dest/inbox/
  done
}

./scripts/clean.sh

ACTOR="1-user"
cd $basedir/$ACTOR
runstep user_init ./scripts/user_process_init.sh "$targets"
dispatch $ACTOR "1-user 2-reference 3-compare 4-ppm 5-product"

ACTOR="2-reference"
cd $basedir/$ACTOR
runstep reference_init ./scripts/reference_process_init.sh
dispatch $ACTOR "1-user 2-reference 3-compare 4-ppm"

ACTOR="3-compare"
cd $basedir/$ACTOR
runstep compare ./scripts/compare_process.sh
dispatch $ACTOR "4-ppm"

ACTOR="4-ppm"
cd $basedir/$ACTOR
runstep ppm ./scripts/ppm_process.sh
dispatch $ACTOR "1-user 2-reference 5-product"

ACTOR=2-reference
cd $basedir/$ACTOR
runstep reference_final ./scripts/reference_process_final.sh
dispatch $ACTOR "5-product"

ACTOR=5-product
cd $basedir/$ACTOR
runstep product ./scripts/product_process.sh
dispatch $ACTOR "1-user"

ACTOR=1-user
cd $basedir/$ACTOR
runstep user_final ./scripts/user_process_final.sh

cd $basedir
duration_stats_file="${prefix}duration_stats.txt"
sum_tmp_file="${prefix}duration_sum.txt"
stats_file="${prefix}stats.csv"
stats_plot_file="${prefix}score.png"
ppm_guess_genotype_score_plot_file="${prefix}4_ppm_guess_genotye_score_plot.pdf"
mv 1-user/outbox/1-user-stats.csv $stats_file

SCORE=$(gawk --field-separator ',' '{s+=$6;u+=$7}END{print "scale=3;" 100*s "/" u}' $stats_file | bc)
echo final imputation score is ${SCORE}%

gawk --field-separator ',' '{ d[$1]+=$3; n[$1]++ }END{ for (i in d) { print i, d[i] } }' */outbox/*-duration_record.csv > $sum_tmp_file
gawk 'NR==FNR{ d[$1]=$2; next } { print $1, d[$1], $2 }' $sum_tmp_file $step_total_duration_file > $duration_stats_file
column -t $duration_stats_file
rm $sum_tmp_file $step_total_duration_file

Rscript scripts/score_plot.R $stats_file $stats_plot_file
Rscript scripts/4_ppm_guess_genotype_score_plot.R $ppm_guess_genotype_score_plot_file
