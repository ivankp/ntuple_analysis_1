#!/bin/bash

dir1=out
dir2=merged

mkdir -p $dir2

for x in $(
  find $dir1/ -name '*.root' |
  sed 's:.*/\(.*\)_[0-9]\+\.root$:\1:' |
  sort -u
); do
  echo $x
  rm -f $dir2/${x}.root
  hadd tmp.root $dir1/${x}_*.root
  ../../bin/merge_root $dir2/${x}.root tmp.root
  rm tmp.root
done

for f in $(
  find merged -name '*.root' |
  sed -r 's/(B|RS|I|V)_/NLO_/' |
  sort -u
); do
  hadd $f $(ls $(sed 's/NLO_/*_/' <<< $f))
done

