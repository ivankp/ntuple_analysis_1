#!/bin/bash

dir1=out
dir2=merged

mkdir -p $dir2

find $dir1 -type f -name '*.root' |
sed -r 's:.*/(.*)_[0-9]+\.root$:\1:' |
sort -u |
while read x
do
  echo $x
  rm -f $dir2/${x}.root
  hadd tmp.root $dir1/${x}_*.root
  ../../bin/merge_root $dir2/${x}.root tmp.root
  rm tmp.root
done

find $dir2 -type f -name '*.root' |
sed -nr 's/([^_])(B|RS|I|V)_/\1NLO_/p' |
sort -u |
while read x
do
  parts=$(sed -r 's/(.*)NLO(_.*)/\1B\2\n\1RS\2\n\1I\2\n\1V\2/' <<< $x)
  if stat --printf='' $parts 2>/dev/null
  then
    ../../bin/merge_root $x $parts
  fi
done
