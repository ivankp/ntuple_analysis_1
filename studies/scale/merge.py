#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', help='force output recreation', action='store_true')
args = parser.parse_args()

import os, re, errno
from collections import defaultdict

def mkdirs(*ds):
    for d in ds:
        try: os.makedirs(d)
        except OSError as e:
            if e.errno == errno.EEXIST: continue
            else: raise

def files(d):
    for (path, dirs, files) in os.walk(d):
        return sorted(files)

def groups(d,r):
    r = re.compile(r)
    gs = defaultdict(list)
    for f in files(d):
        m = r.match(f)
        if m: gs['_'.join(m.groups())].append(os.path.join(d,f))
    return gs.items()

merger = '/home/ivanp/work/ntuple_analysis/bin/merge_sql'
mkdirs('merged')

for (g,fs) in groups('out',r'^(.+)_\d+\.db$'):
    print '\n', g
    ofname = 'merged/'+g+'.db'
    if args.f or not os.path.isfile(ofname):
        os.system(merger+' -n -o '+ofname+' '+' '.join(fs))

re_nlo = re.compile(r'^([^\d]+\d+j)(B|RS|I|V)(.*)\.db$')
nlo = defaultdict(list)
for f in files('merged'):
    m = re_nlo.match(f)
    if m:
        g = m.groups()
        nlo[(g[0],g[2])].append((g[1],f))
for (g,fs) in nlo.items():
    if reduce(
        lambda a, p: a and filter(lambda x: x[0]==p, fs),
        ['B','RS','I','V'],
        len(fs)==4
    ):
        ofname = g[0]+'NLO'+g[1]
        print '\n', ofname
        ofname = 'merged/'+ofname+'.db'
        if args.f or not os.path.isfile(ofname):
            os.system(merger+' -x -o '+ofname+' '+
                ' '.join('merged/'+f[1] for f in fs))

ofname = 'merged/all.db'
print '\n', ofname
if args.f or not os.path.isfile(ofname):
    os.system('./join.py '+ofname+' merged/*.db')

