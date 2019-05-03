#!/bin/bash

if [ $# -ne 1 ]; then
  echo "usage: $0 info_string"
  exit 1
fi

sqlite3 /home/ivanp/work/ntuple_analysis/sql/ntuples.db "
select count(id),sum(nentries)/1e6,info,dir
from ntuples
where
  info=\"$1\"
  and energy=13
  and part=\"B\"
  and njets=2
"

