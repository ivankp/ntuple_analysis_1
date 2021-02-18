#!/bin/bash

for ene in 13 100; do
for zoom in 1 2; do
  hed merged/H2jB_${ene}TeV_{eft,mtop}_all_antikt4.root \
    -o H2jB_${ene}TeV_zoom${zoom}.pdf \
    -e 'fl/.*_(eft|mtop)_.*/\1/' 'sd:weight2/all/all/no:' \
       "snt/sqrt_s34$/" "t+//  ${ene} TeV/" \
    -g 'rat 0.5 width; log y; leg' \
    --colors 602 98
done
done
