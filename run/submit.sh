#!/bin/bash

rm -f ~/.condor_cp_mutex

for dir in $@; do
  echo $dir
  cd $dir
  condor_submit ../condor
done
