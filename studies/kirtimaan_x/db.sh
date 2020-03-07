#!/bin/bash

../../bin/root2sql \
  -i merged/*.root \
  -o kirtimaan_x.db \
  --hist-regex='.+' \
  -l proc type diag jet weight \
     $(sed -n 's/^#define CATEGORIES //p' analysis.cc | sed 's/[()]\+/ /g') \
     var1

sqlite3 kirtimaan_x.db << SQL
  DELETE FROM hist WHERE
    (swap_45="all" AND ( var1 LIKE "x\\_%" ESCAPE "\\"
      OR var1 IN ("s34","s45","t15")
    )) OR
    (swap_45!="all" AND ( var1 LIKE "Njets\\_%" ESCAPE "\\"
      OR var1 IN ("s12","t23","x1","x2")
    ));

  CREATE TEMPORARY TABLE temp AS
    SELECT * FROM hist WHERE swap_45="all";
  UPDATE temp SET swap_45="no";
  INSERT INTO hist SELECT * FROM temp;
  DROP TABLE temp;

  CREATE TEMPORARY TABLE temp AS
    SELECT * FROM hist WHERE swap_45="all";
  UPDATE temp SET swap_45="yes";
  INSERT INTO hist SELECT * FROM temp;
  DROP TABLE temp;

  DELETE FROM hist WHERE swap_45="all";
SQL
