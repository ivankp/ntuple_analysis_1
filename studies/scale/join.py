#!/usr/bin/env python

import sys, os, re, sqlite3

if len(sys.argv)<3:
    print "usage:", sys.argv[0], "out.db in1.db [in2.db ...]"

fnames = sys.argv[2:]
curs = [ sqlite3.connect(x).cursor() for x in fnames ]

tags = [ f[f.rfind('/')+1:].rstrip('.db').split('_') for f in fnames ]
ntags = len(tags[0])
if not all(len(x)==ntags for x in tags):
    print "file names must split into equal numbers of tags"

def all_same(xs, f=lambda x: x):
    x0 = None
    for x in xs:
        if x0 is None: x0 = f(x)
        elif f(x) != x0: return False
    return True

# print tags
for i in filter(lambda i: all_same(t[i] for t in tags), xrange(ntags)):
    for t in tags:
        t.pop(i)
print tags
ntags = len(tags[0])

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

for tab in tabs:
    if tab=='hist': continue
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

def test(s):
    print s
    return s

out_db = sqlite3.connect(sys.argv[1])
out_cur = out_db.cursor()
for i,f in enumerate(fnames):
    out_cur.execute(test('ATTACH DATABASE "{}" AS db{}'.format(f,i)))
out_cur.execute(test(
    'CREATE TABLE hist(' +
    ','.join('\n tag{} TEXT'.format(i) for i in xrange(ntags)) +
    ',' +
    ','.join('\n {} {}'.format(*col) for col in cols) +
    '\n)'
))
for i,f in enumerate(fnames):
    out_cur.execute(test(
        'INSERT INTO hist SELECT {}, * FROM db{}.hist'.format(
        ', '.join('"{}"'.format(tag) for tag in tags[i]), i)
    ))
    out_db.commit()

