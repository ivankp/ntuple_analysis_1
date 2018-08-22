#!/bin/bash

source ../env.sh

echo "$@"
echo "pwd -P: $(pwd -P)"
echo "TMPDIR = $TMPDIR"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"

ofname=$(python -c "
import json
with open('$1','r') as f:
  j = json.load(f)
  print j['output']
")
echo "ofname = $ofname"

path='/home/ivanp/work/ntuple_analysis'
${path}/bin/analyses/test ${path}/runcards/test.json $1

# xz $ofname

