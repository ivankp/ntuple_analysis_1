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

chunk_size = 50e6

mkdirs(loc+'/condor',loc+'/out')

db = sqlite3.connect(project_dir+'/sql/ntuples.db')

LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH']

subcount = defaultdict(lambda:0)

def diagram(info):
    for m in re.finditer(' (tr|bx|pn)(?: |$)',info):
        return m.group(1)
    return 'all'

def get(names,vals):
    fs = [ ( x[-1], x[0]+'/'+x[1], x[3],
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
    for x in fs:
        if n == 0:
            subcount[pref] += 1
            chunks.append(('{}_{:0>3d}'.format(pref,subcount[pref]),[],x[2]))
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
                "njets_min": chunk[2]
            },
            'binning': loc+'/analysis.bins'
        },
        'output': loc+'/out/'+chunk[0]+'.root'
    }, indent=2, separators=(',',': ')) ))

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

os.chdir(loc+'/condor')

params = zip(
    ('njets',(1,2,3)),
    ('part',('B','RS','I','V')),
    ('particle',('H','AA')),
    ('energy',(8,13,100))
)

infos = [
    ({ 'particle': 'H', 'energy': 8 },
     [ 'GGFHT pt25.0 eta4.5' ]),

    ({ 'njets': 1, 'particle': 'H', 'energy': 13 },
     [ 'GGFHT pt25.0 eta4.5', 'mtop GGFHT pt25.0 eta4.5',
       'amegic_GGFHT pt25.0 eta4.5 mtop tr',
       'amegic_GGFHT pt25.0 eta4.5 mtop bx' ]),

    ({ 'njets': 2, 'particle': 'H', 'energy': 13 },
     [ 'ED GGFHT pt25.0 eta4.5', 'mtop GGFHT pt25.0 eta4.5',
       'amegic_GGFHT pt25.0 eta4.5 mtop tr',
       'amegic_GGFHT pt25.0 eta4.5 mtop bx',
       'amegic_GGFHT pt25.0 eta4.5 mtop pn' ]),

    ({ 'njets': 3, 'particle': 'H', 'energy': 13 },
     [ 'ED GGFHT pt25.0 eta4.5', 'mtop GGFHT pt25.0 eta4.5' ]),

    ({ 'particle': 'H', 'energy': 100, 'njets': 1 },
     [ 'eos amegic_GGFHT pt25.0 eta10.0', 'eos GGFHT pt25.0 eta10.0',
       'eos mtop amegic_GGFHT pt25.0 eta10.0' ]),

    ({ 'particle': 'H', 'energy': 100, 'njets': 2 },
     [ 'eos amegic_GGFHT pt25.0 eta10.0', 'eos GGFHT pt25.0 eta10.0',
       'eos mtop GGFHT_FCC pt25.0 eta10.0' ]),

    ({ 'particle': 'H', 'energy': 100, 'njets': 2 },
     [ 'eos GGFHT pt25.0 eta4.5' ]),

    ({ 'particle': 'AA', 'energy': 13 },
     [ 'Nico' ])
]

for vals in product(*params[1]):
    for info_key, info in infos:
        matched = True
        for key, val in zip(params[0],vals):
            val2 = info_key.get(key)
            if val2 is not None and val2 != val:
                matched = False
                break
        if matched:
            for x in info:
                for chunk in get(params[0]+('info',),vals+(x,)):
                    outf = loc + '/condor/' + chunk[0] + '.out'
                    if os.path.exists(outf):
                        continue

                    print chunk[0]
                    job = condor(chunk)

                    p = Popen(('condor_submit','-'), stdin=PIPE, stdout=PIPE)
                    p.stdin.write(job)
                    p.communicate()
                    p.stdin.close()

