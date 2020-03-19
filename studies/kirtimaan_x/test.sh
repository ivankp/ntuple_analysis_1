#!/bin/bash
export LD_LIBRARY_PATH=/msu/data/t3work3/ivanp/gcc-9.2.0/lib64:/msu/data/t3work3/ivanp/gcc-9.2.0/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/hep/root-6.10.02/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/hep/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/gcc/lib64:/msu/data/t3work3/ivanp/gcc-7.2.0/gcc/lib:/msu/data/t3work3/ivanp/gcc-7.2.0/lib64:/msu/data/t3work3/ivanp/gcc-7.2.0/lib

/home/ivanp/work/ntuple_analysis/studies/kirtimaan_x/analysis - << CARD
{
  "input": [
    {
      "files": [
        "/msu/data/t3work4/luisonig/H2jets_ggf/EDNTuplesFiles/H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_100.root"
      ]
    }
  ],
  "analysis": {
    "jets": {
      "alg": [
        "antikt",
        0.4
      ],
      "njets_min": 2,
      "cuts": {
        "eta": 4.4,
        "pT": 30
      }
    },
    "binning": "/home/ivanp/work/ntuple_analysis/studies/kirtimaan_x/analysis.bins"
  },
  "output": "/home/ivanp/work/ntuple_analysis/studies/kirtimaan_x/out/H2jB_eft_test.root"
}
CARD
