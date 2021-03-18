#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o mass_vs_pT.db \
  -l proc energy type diag jet weight \
     $(sed -n 's/^#define CATEGORIES //p' analysis.cc | sed 's/[()]\+/ /g') \
     var1 var2 var3

  # --hist-regex '((Hj[12]|jj)_mass)_(H|j[12])_pT'

