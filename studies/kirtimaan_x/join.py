#!/usr/bin/env python

import argparse, sys, os, re, sqlite3

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='ifnames', nargs='+', required=True,
    help='input files')
parser.add_argument('-o', dest='ofname', required=True,
    help='output file')
parser.add_argument('-t', dest='tags', nargs='+', required=True,
    help='regex for tokenizing input file names')
args = parser.parse_args()

tags_names = args.tags[1:]
ntags = len(tags_names)
tags_re = re.compile(args.tags[0])

tags = [ ]
for f in args.ifnames:
    m = tags_re.match(f)
    if not m:
        print "input file names must match tags expression"
        sys.exit(1)
    m = list(m.groups())
    print m
    if len(m)!=ntags:
        print "number of match groups must equal the number of tag names"
        sys.exit(1)
    tags.append(m)

def all_same(xs, f=lambda x: x):
    x0 = None
    for x in xs:
        if x0 is None: x0 = f(x)
        elif f(x) != x0: return False
    return True

tagi = filter(lambda i: not all_same(t[i] for t in tags), xrange(ntags))
tags = [ [ t[i] for i in tagi ] for t in tags ]
tags_names = [ tags_names[i] for i in tagi ]
ntags = len(tags_names)
print tags_names
print tags

curs = [ sqlite3.connect(x).cursor() for x in args.ifnames ]

tabs = None
for cur in curs:
    xs = cur.execute(
        'SELECT name FROM sqlite_master WHERE type="table" order by name')
    if not tabs: tabs = xs.fetchall()
    else:
        for i,x in enumerate(xs):
            if x != tabs[i]:
                print "all input databases must have the same tables"
                sys.exit(1)
tabs = [ t[0] for t in tabs ]
tabs.remove('hist')
print tabs

for tab in tabs:
    for rows in zip(cur.execute('SELECT * FROM '+tab) for cur in curs):
        x = None
        for row in rows:
            if not x: x = row
            elif x != row:
                print "unequal auxiliary table"
                sys.exit(1)

cols = None
for cur in curs:
    xs = [ (c[1],c[2]) for c in cur.execute('PRAGMA table_info(hist)') ]
    if not cols: cols = xs
    elif cols != xs:
        print "different columns in table \"hist\""
        sys.exit(1)
print cols

print args.ofname
if os.path.isfile(args.ofname):
    os.remove(args.ofname)

out_db = sqlite3.connect(args.ofname)
out_cur = out_db.cursor()
def exe(sql):
    print sql
    out_cur.execute(sql)

for i,f in enumerate(args.ifnames):
    exe('ATTACH DATABASE "{}" AS db{}'.format(f,i))
    if i==0:
        for tab in tabs:
            exe('CREATE TABLE {0} AS SELECT * FROM db0.{0}'.format(tab))
        exe('CREATE TABLE hist(' +
            ','.join('\n {} TEXT'.format(tag) for tag in tags_names) +
            ',' +
            ','.join('\n {} {}'.format(*col) for col in cols) +
            '\n)')
        out_db.commit()
    exe('INSERT INTO hist SELECT {}, * FROM db{}.hist'.format(
        ', '.join('"{}"'.format(tag) for tag in tags[i]), i))
    out_db.commit()
    exe('DETACH DATABASE db{}'.format(i))

