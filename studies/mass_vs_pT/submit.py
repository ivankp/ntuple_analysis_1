#!/usr/bin/env python

import sys, os, sqlite3, json, re
from subprocess import Popen, PIPE
from collections import defaultdict
from itertools import product

def mkdirs(*ds):
    for d in ds:
        if not os.path.isdir(d):
            os.makedirs(d)

loc = os.path.dirname(os.path.realpath(__file__))
exe = loc + '/analysis'

project_dir = loc
while os.path.basename(project_dir) != 'ntuple_analysis':
    project_dir = os.path.dirname(project_dir)

jetR = 4

chunk_size = 20e6

mkdirs(loc+'/condor',loc+'/out')

db = sqlite3.connect(project_dir+'/sql/ntuples.db')

LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH']

subcount = defaultdict(lambda:0)

def diagram(info):
    for m in re.finditer(' (tr|bx|pn)(?: |$)',info):
        return m.group(1)
    return 'all'

def get(names,vals):
    fs = [ ( x[-1], x[0]+'/'+x[1], x[3], x[5],
        '{}{}j{}_{:g}TeV_{}_{}_antikt{:g}'.format(
            x[2], x[3], x[4],
            x[5],
            ('mtop' if ('mtop' in x[6]) else 'eft'),
            diagram(x[6]),
            jetR
        )
    ) for x in db.execute('''
SELECT dir,file,particle,njets,part,energy,info,nentries
FROM ntuples
WHERE
'''+' and '.join(a+'=?' for a in names),vals).fetchall() ]

    pref = set([f[-1] for f in fs])
    if len(pref) > 1:
        raise Exception('multiple types in single selection: '+' '.join(pref))
    elif len(pref)==0:
        return [ ]
    pref = pref.pop()

    chunks = [ ]
    n = 0
    for f in fs:
        if n == 0:
            subcount[pref] += 1
            chunks.append((
                '{}_{:0>3d}'.format(pref,subcount[pref]), [], f[3], f[2]
            ))
        chunks[-1][1].append(f[1])
        n += f[0]
        if n >= chunk_size:
            n = 0

    return chunks

def make_job(chunk):
    script = chunk[0]+'.sh'
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
            'rootS': chunk[2],
            'jets': {
                "cuts": { "pT": 30, "eta": 4.4 },
                "alg": [ "antikt", jetR*0.1 ],
                "njets_min": chunk[3]
            },
            'binning': loc+'/analysis.bins'
        },
        'output': loc+'/out/'+chunk[0]+'.root'
    }, indent=2, separators=(',',': ')) ))

    os.chmod(script,0o775)

os.chdir(loc+'/condor')

params = list(zip(
    ('njets',(2,)),
    ('part',('B',)),
    ('particle',('H','AA')),
    ('energy',(13,100)),
    ('info',(
        'Nico',
        'ED GGFHT pt25.0 eta4.5',
        'mtop GGFHT pt25.0 eta4.5',
        'mtop amegic_GGFHT pt25.0 eta4.5',
        'amegic_GGFHT pt25.0 eta4.5 mtop tr',
        'amegic_GGFHT pt25.0 eta4.5 mtop bx',
        'amegic_GGFHT pt25.0 eta4.5 mtop pn',
        'eos GGFHT pt25.0 eta4.5',
        'eos GGFHT pt25.0 eta10.0',
        'eos amegic_GGFHT pt25.0 eta10.0',
        'eos mtop GGFHT_FCC pt25.0 eta10.0'
    ))
))

with open('finish.sh','w') as f:
    f.write('''\
#!/bin/bash
export LD_LIBRARY_PATH={lib}
cd ..
./merge.sh
./db.sh
'''.format(lib = LD_LIBRARY_PATH))
os.chmod('finish.sh',0o775)

with open('job.sub','w') as f:
    f.write('''\
Universe   = vanilla
Executable = $(name).sh
Output     = $(name).out
Error      = $(name).err
Log        = $(name).log
getenv = True
queue
''')

with open('jobs.dag','w') as f:
    n = 0
    for vals in product(*params[1]):
        for chunk in get(params[0],vals):
            print(chunk[0])
            make_job(chunk)
            f.write('JOB j{0} job.sub\nVARS j{0} name="{1}"\n\n'.format(
                n, chunk[0]
            ))
            n += 1

    f.write('''\
JOB finish job.sub
VARS finish name="finish"

PARENT '''+(' '.join('j{}'.format(i) for i in range(n)))+' CHILD finish\n')

Popen(('condor_submit_dag','jobs.dag')).communicate()

