#!/bin/bash

source ../env.sh

echo "@ = $@"
echo "pwd -P: $(pwd -P)"
echo "TMPDIR = $TMPDIR"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"

# while [ -a "~/.condor_cp_mutex" ]
# do
#   echo waiting
#   sleep 1
# done

# sleep $(bc <<< "0.1*$1")
sleep $(python -c "print 0.5*(${1}%25)")

../../bin/analyses/hist_Hjets ../../cards/Hjets.json $2 --tmp-dir=$TMPDIR

exit_code=$?
if [ $exit_code -ne 0 ]; then
  printf "\033[31;1m"
else
  printf "\033[32;1m"
  mv -v $2 $(sed -r 's/(_card)(\.json)$/\1_done\2/' <<< "$2")
fi
printf "Analysis finished\033[0m\n"

rm -rfv $TMPDIR/*

