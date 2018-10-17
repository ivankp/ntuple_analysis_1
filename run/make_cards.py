#!/usr/bin/env python

import os, sqlite3, json

def mkdir(d):
    if not os.path.isdir(d): os.makedirs(d)

path = '/home/ivanp/work/ntuple_analysis'
db  = sqlite3.connect(path+'/sql/ntuples.db')
cur = db.cursor()

# odir = path+'/run/Hjets/'
odir = '/msu/data/t3work9/ivanp/ntuples_hists/rel/raw/'
if odir[-1]!='/': odir += '/'
print odir
mkdir(odir)

# sets = set()

for x in cur.execute('''
SELECT dir,file,id,particle,njets,part,info
FROM ntuples
WHERE energy=13
  and not instr(info,"GGFHTLS")
  and not instr(info,"VBFHT4VBF_INCL")
  and not instr(info,"VBFHT4VBF_PTH200")
'''):
    name = '{3}{4}j_{5}'.format(*x)
    # sets.add(name)
    name = '{}_{}'.format(x[2],name)
    info = x[-1].split()
    if 'ED' in info: name += '_ED'
    if 'mtop' in info: name += '_mtop'
    if 'VBFHT4VBF' in info: name += '_VBF'
    print name
    name = odir+name
    with open(name+'_card.json','w') as card:
        json.dump({
            "input": [{ "files": [ x[0]+'/'+x[1] ] }],
            "output": name+'_hist.json.xz',
            "analysis": { "jets": { "njets_born": x[4] } }
        },card,separators=(',',':'))

# for s in sets:
#     with open(odir+'jobs.dag','w') as dag:
#         dag.write('''\
# JOB {0} ../analysis.condor
# VARS {0} match_cards="*_{0}_card.json"
# SCRIPT POST {1}/bin/merge -o merged_{0}_hist.json.xz *_{0}_hist.json.xz
#
# '''.format(s,path))
