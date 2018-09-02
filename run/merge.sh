#!/bin/bash

../bin/merge -v -o H1j_B_hist.json.xz Hjets/*_H1j_B_hist.json.xz
../bin/merge -x -o H1j_B.json.xz H1j_B_hist.json.xz

