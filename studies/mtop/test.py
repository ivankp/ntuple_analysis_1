#!/usr/bin/env python

import sqlite3, struct

db = sqlite3.connect('eft_mtop_H2jB.db')

hs = db.execute('''
select var1, var2, photon_cuts
from hist
where
type="eft" and
jet="antikt4" and
isp="all" and
photon_cuts="all" and
higgs_y_cut="central_higgs" and
nsubjets="all" and
flavor1="all" and
flavor2="all" and
var1="Hj_mass" and
var2=""
''').fetchall()

print hs

