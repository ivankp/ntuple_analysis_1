#!/bin/bash

/home/ivanp/work/ntuple_analysis/studies/scale/analysis - << CARD
{
  "input": [
    {
      "files": [
        "/msu/data/t3work4/luisonig/H1jets_ggf/NTuplesFiles/H1.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_100.root"
      ]
    }
  ],
  "entry_range": [0,100000],
  "reweighting": [
    {
      "pdf": "CT14nlo",
      "ren_fac": [
        [
          1,
          1
        ],
        [
          0.5,
          0.5
        ],
        [
          1,
          0.5
        ],
        [
          0.5,
          1
        ],
        [
          2,
          1
        ],
        [
          1,
          2
        ],
        [
          2,
          2
        ]
      ],
      "scale": "HT1",
      "pdf_var": false
    },
    {
      "pdf": "CT14nlo",
      "ren_fac": [
        [
          1,
          1
        ],
        [
          0.5,
          0.5
        ],
        [
          1,
          0.5
        ],
        [
          0.5,
          1
        ],
        [
          2,
          1
        ],
        [
          1,
          2
        ],
        [
          2,
          2
        ]
      ],
      "scale": "HT2",
      "pdf_var": false
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
    "binning": "/home/ivanp/work/ntuple_analysis/studies/scale/analysis.bins"
  },
  "output": "/home/ivanp/work/ntuple_analysis/studies/scale/test.db"
}
CARD
