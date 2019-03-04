#!/bin/bash

~/work/bh_analysis2/ang/root_hist_db/bin/make_db \
  H2jB_{eft,mtop}_antikt{4,6,8}.root \
  -o eft_mtop_H2jB.db \
  -l proc type jet weight isp photon_cuts higgs_y_cut \
     nsubjets flavor1 flavor2 var1 var2

