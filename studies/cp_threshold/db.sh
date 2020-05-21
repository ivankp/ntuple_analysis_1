#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o cp_threshold.db \
  -l proc energy type diag jet weight \
     $(sed -n 's/^#define CATEGORIES //p' analysis.cc | sed 's/[()]\+/ /g') \
     var1 var2

