#!/usr/bin/env python

import os, sqlite3, json

def mkdir(d):
    if not os.path.isdir(d): os.makedirs(d)

path = '/home/ivanp/work/ntuple_analysis'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

odir = path+'/run/Hjets'
print odir
mkdir(odir)

for x in cur.execute('''
SELECT dir,file,id,particle,njets,part
FROM ntuples
WHERE dir="/msu/data/t3work4/luisonig/H1jets_ggf/NTuplesFiles"
  and njets=1 and energy=13 and part="B"
'''):
    name = '{2}_{3}{4}j_{5}'.format(*x)
    print name
    name = odir+'/'+name
    with open(name+'_card.json','w') as ofile:
        json.dump({
            "input": [{ "files": [ x[0]+'/'+x[1] ] }],
            "output": name+'_hist.json.xz'
        },ofile,separators=(',',':'))

