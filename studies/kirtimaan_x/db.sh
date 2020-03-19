#!/bin/bash

rm -fv kirtimaan_x.db

../../bin/root2sql \
  -i merged/*.root \
  -o kirtimaan_x.db \
  --hist-regex='.+(?=_zoom\d*$)|zoom\d*|(?!_zoom\d*).+' \
  -l proc type diag jet weight \
     $(sed -n 's/^#define CATEGORIES //p' analysis.cc | sed 's/[()]\+/ /g') \
     var1 zoom

# https://regex101.com/r/yN4tJ6/343

sqlite3 kirtimaan_x.db << SQL
  DELETE FROM hist WHERE
    (swap_45="all" AND ( var1 LIKE "x\\_%" ESCAPE "\\"
      OR var1 IN ("s34","s35","s45",
            "t15","t34","t35","t45")
      OR var1 IN ("sqrt_s34","sqrt_s35","sqrt_s45",
       "sqrt_t15","sqrt_t34","sqrt_t35","sqrt_t45")
    )) OR
    (swap_45!="all" AND ( var1 LIKE "Njets\\_%" ESCAPE "\\"
      OR var1 IN ("s12","t23","x1","x2")
      OR var1 IN ("sqrt_s12","sqrt_t23","x1","x2")
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
