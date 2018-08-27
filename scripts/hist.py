#!/usr/bin/env python3

import sys, json, lzma

argc = len(sys.argv)
data = None

with lzma.open(sys.argv[1]) if sys.argv[1].endswith('.xz') \
else open(sys.argv[1]) as f:
    data = json.load(f)

hists = data['histograms']
bins_fmt = data['annotation']['bins']

if argc==2:
    print('\033[34;1mHistograms:\033[0m')
    for hist_name in hists:
        print(hist_name)
    print('\n\033[34;1mBin values:\033[0m')
    for x in bins_fmt:
        print('{}: {}'.format(x[0],', '.join(x[1])))
    sys.exit(0)

ii = [ bins_fmt[i][1].index(x) for i,x in enumerate(sys.argv[3].split(':')) ] \
     if argc>3 else [ ]

hist = hists[sys.argv[2]]
bins = [ ]
for b in hist["bins"]:
    for i in ii:
        if b is None: break
        b = b[i]
    bins.append(b)

print(json.dumps({ 'axes': hist['axes'], 'bins': bins }, separators=(',',':')))

