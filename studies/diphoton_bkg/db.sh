#!/bin/bash

rm -fv diphoton_bkg.db

../../bin/root2sql \
  -i merged/*.root \
  -o diphoton_bkg.db \
  --hist-regex='[^_]+(?:_yy|_zoom|_incl|_excl|\[[^)]+\))*' \
  -l proc type diag jet weight photon_cuts isp var1 var2

sqlite3 diphoton_bkg.db << SQL
  CREATE TABLE tmp(
    proc TEXT,
    jet TEXT,
    photon_cuts TEXT,
    isp TEXT,
    var1 TEXT,
    var2 TEXT,
    axis INTEGER,
    bins TEXT
  );
  INSERT INTO tmp
  SELECT
    proc, jet, photon_cuts, isp, var1, var2, axis, bins
  FROM hist;
  DROP TABLE hist;
  ALTER TABLE tmp RENAME TO hist;
SQL

