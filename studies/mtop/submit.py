#!/usr/bin/env python

import sys, os, sqlite3, json
from subprocess import Popen, PIPE
from collections import defaultdict

def mkdirs(*ds):
    for d in ds:
        if not os.path.isdir(d):
            os.makedirs(d)

project_dir = '/home/ivanp/work/ntuple_analysis'
loc = project_dir+'/studies/mtop'
exe = loc + '/analysis2'

jetR = float(sys.argv[1]) if len(sys.argv)>1 else 4
njets = 2

chunk_size = 20e6

mkdirs(loc+'/condor',loc+'/out')

db  = sqlite3.connect(project_dir+'/sql/ntuples.db')
cur = db.cursor()

LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH']

subcount = defaultdict(lambda:0)

def get(info,nj):
    cur.execute('''
SELECT dir,file,particle,njets,part,info,nentries
FROM ntuples
WHERE info="{}" and njets={} and particle="H" and energy=13 and part="B"
'''.format(info,nj))
    fs = [ ( x[-1], x[0]+'/'+x[1], '{}{}j{}_{}_antikt{:g}'.format(
        x[2], x[3], x[4], ('mtop' if ('mtop' in x[5]) else 'eft'), jetR
    )) for x in cur.fetchall() ]
    pref = set([x[2] for x in fs])
    if len(pref) != 1:
        raise Exception('multiple types in single selection: '+' '.join(pref))
    pref = pref.pop()

    chunks = [ ]
    n = 0
    for x in fs:
        if n == 0:
            subcount[pref] += 1
            chunks.append(('{}_{:0>3d}'.format(pref,subcount[pref]),[]))
        chunks[-1][1].append(x[1])
        n += x[0]
        if n >= chunk_size:
            n = 0

    return chunks

def condor(chunk):
    script = loc+'/condor/'+chunk[0]+'.sh'
    with open(script,'w') as f:
        f.write('''\
#!/bin/bash
export LD_LIBRARY_PATH={lib}\n
{exe} - << CARD
{card}
CARD
'''.format(
    lib = LD_LIBRARY_PATH,
    exe = exe,
    card = json.dumps({
'input': [{ 'files': chunk[1] }],
'analysis': {
  'jets': {
    "cuts": { "pT": 30, "eta": 4.4 },
    "alg": [ "antikt", jetR*0.1 ],
    "njets_min": 1
  },
  'binning': exe+'.bins'
  },
'output': loc+'/out/'+chunk[0]+'.root'
}, indent=2, separators=(',',': '))
))

    os.chmod(script,0o775)
    return '''\
Universe   = vanilla
Executable = {0}.sh
Output     = {0}.out
Error      = {0}.err
Log        = {0}.log
getenv = True
Queue
'''.format(chunk[0])

def jobs(*infos):
    for info in infos:
        for chunk in get(info,njets):
            print chunk[0]
            yield condor(chunk)

os.chdir(loc+'/condor')

for job in jobs(
'ED GGFHT pt25.0 eta4.5',
'mtop GGFHT pt25.0 eta4.5',
'mtop amegic_GGFHT pt25.0 eta4.5'
):
    p = Popen(('condor_submit','-'), stdin=PIPE, stdout=PIPE)
    p.stdin.write(job)
    p.communicate()
    p.stdin.close()

