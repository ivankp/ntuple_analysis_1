#!/bin/bash

source ../env.sh

echo "$@"
echo "pwd -P: $(pwd -P)"
echo "TMPDIR = $TMPDIR"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"

# echo "ls from job.sh"
# ls -lh /home/ivanp/.condor_cp_mutex
# ls -lh .condor_cp_mutex
while [ -a "~/.condor_cp_mutex" ]
do
  echo waiting
  sleep 1
done

# touch ~/.condor_cp_mutex

../../bin/analyses/test ../../cards/test.json $1 --tmp-dir=$TMPDIR

exit_code=$?
if [ $exit_code -ne 0 ]; then
  printf "\033[31;1m"
else
  printf "\033[32;1m"
fi
printf "Analysis finished\033[0m\n"

rm -fv $TMPDIR/*

