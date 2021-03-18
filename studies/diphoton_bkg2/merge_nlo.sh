#!/bin/bash

find merged -type f -name '*.root' |
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
