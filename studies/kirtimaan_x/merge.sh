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

