#!/usr/bin/env python3

import sys, json, lzma

argc = len(sys.argv)
data = None

if argc<2:
    print("usage: {} hists.json[.xz] select.json".format(sys.argv[0]))
    sys.exit(0)

with lzma.open(sys.argv[1]) if sys.argv[1].endswith('.xz') \
else open(sys.argv[1]) as f:
    data = json.load(f)

hists      = data['histograms']
categories = data['annotation']['bins'][:-1]
weights    = data['annotation']['weights']

if argc<3:
    print('\033[34;1mHistograms:\033[0m')
    for hist_name in hists:
        print(hist_name)
    print('\n\033[34;1mCategories:\033[0m')
    for x in categories:
        print('{}: {}'.format(x[0],', '.join(x[1])))
    print('\n\033[34;1mWeights:\033[0m')
    for x in weights:
        print(x)
    sys.exit(0)

select = None
with open(sys.argv[2]) as f:
    select = json.load(f)

hist = hists[select["hist"]]
ci   = [ c[1].index(select[c[0]]) for c in categories ]
wi   = weights.index(select["weight"])

bins = [ ]
for b in hist["bins"]:
    for i in ci:
        if b is None: break
        b = b[i]
    bins.append([ b[0][wi], b[1] ])

print(json.dumps({ 'axes': hist['axes'], 'bins': bins },separators=(',',':')))

