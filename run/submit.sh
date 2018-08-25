#!/bin/bash

for dir in $@; do
  echo $dir
  cd $dir
  condor_submit ../condor
done
