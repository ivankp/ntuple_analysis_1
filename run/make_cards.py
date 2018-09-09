#!/usr/bin/env python

import os, sqlite3, json

def mkdir(d):
    if not os.path.isdir(d): os.makedirs(d)

path = '/home/ivanp/work/ntuple_analysis'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

odir = path+'/run/Hjets/'
print odir
mkdir(odir)

# sets = set()

for x in cur.execute('''
SELECT dir,file,id,particle,njets,part
FROM ntuples
WHERE dir="/msu/data/t3work4/luisonig/H1jets_ggf/NTuplesFiles"
  and njets=1 and energy=13 and part="B"
'''):
    name = '{3}{4}j_{5}'.format(*x)
    # sets.add(name)
    name = '{}_{}'.format(x[2],name)
    print name
    name = odir+name
    with open(name+'_card.json','w') as card:
        json.dump({
            "input": [{ "files": [ x[0]+'/'+x[1] ] }],
            "output": name+'_hist.json.xz'
        },card,separators=(',',':'))

# for s in sets:
#     with open(odir+'jobs.dag','w') as dag:
#         dag.write('''\
# JOB {0} ../analysis.condor
# VARS {0} match_cards="*_{0}_card.json"
# SCRIPT POST {1}/bin/merge -o merged_{0}_hist.json.xz *_{0}_hist.json.xz
#
# '''.format(s,path))
