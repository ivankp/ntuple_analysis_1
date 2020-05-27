#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o Hj_basic.db \
  --hist-regex='^.*' \
  -l proc energy type diag jet weight photon_cuts isp var1

