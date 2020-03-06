#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o kirtimaan_x.db \
  --hist-regex='.+' \
  -l proc type diag jet weight \
     $(sed -n 's/^#define CATEGORIES //p' analysis.cc | sed 's/[()]\+/ /g') \
     var1

sqlite3 kirtimaan_x.db << SQL
  delete from hist where
  (swap_45="all" and ( var1 like "x\\_%" ESCAPE "\\"
    or var1 in ("s34","s45","t15")
  )) or
  (swap_45!="all" and ( var1 like "Njets\\_%" ESCAPE "\\"
    or var1 in ("s12","t23","x1","x2")
  ))
SQL

