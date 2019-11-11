#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o cp_threshold.db \
  -l proc type jet weight photon_cuts isp \
     var1 var2

