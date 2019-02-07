#!/bin/bash

dir=out

for x in $(
  find $dir -name '*.root' |
  sed 's:.*/\(.*\)_[0-9]\+\.root$:\1:' |
  sort -u
); do
  echo $x
  rm -f ${x}.root
  hadd tmp.root $dir/${x}_*.root
  /home/ivanp/work/ntuple_analysis/bin/merge_root ${x}.root tmp.root
  rm tmp.root
done

for x in $(
  ls *.root |
  sed 's/_\(eft\|mtop\)_/__/' |
  sort -u
); do
  hrat \
    $(sed 's/__/_mtop_/' <<< $x) \
    $(sed 's/__/_eft_/' <<< $x) \
    $(sed 's/__/_rat_/' <<< $x)
done

