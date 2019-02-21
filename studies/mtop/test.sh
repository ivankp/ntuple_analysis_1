#!/bin/bash

/home/ivanp/work/ntuple_analysis/studies/mtop/analysis - << CARD
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
      "njets_min": 1,
      "cuts": {
        "eta": 4.4,
        "pT": 30
      }
    },
    "binning": "/home/ivanp/work/ntuple_analysis/studies/mtop/mtop.bins"
  },
  "output": "test.root"
}
CARD
