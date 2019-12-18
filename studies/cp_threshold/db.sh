#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o cp_threshold.db \
  -l proc type diag jet weight isp photon_cuts T24_sqrtS \
     var1 var2

